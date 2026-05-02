"""
Co-simulacion: acopla el mecanismo (ODE) con la sarta (PDE) y la bomba.

Esquema de paso compartido (Seccion 4.3 del documento):
1. Estado conocido: theta^n, omega^n, u_i^n, h_din^n
2. Calcular y_PR^n y TF^n desde la cinematica
3. Avanzar la PDE de la sarta un paso
4. Obtener F_PR^n del solver
5. Actualizar bomba (F_pump, q)
6. Calcular T_motor y T_neto
7. Avanzar theta y omega con RK4
"""

import numpy as np
from aib_simulator.geometry import (
    compute_kinematics, polished_rod_position,
    compute_kinematics_mark2,
)
from aib_simulator.rod_string import WaveSolver
from aib_simulator.pump import PumpModel
from aib_simulator.motor import torque_linear, torque_klauss
from aib_simulator.mechanism import (
    equivalent_inertia, counterweight_torque, compute_acceleration
)


class Simulator:
    """Simulador acoplado: motor + mecanismo + sarta + bomba.

    Parameters
    ----------
    params : dict — todos los parametros del sistema (de params.py)
    motor_model : str — 'linear' o 'klauss'
    """

    def __init__(self, params, motor_model='linear'):
        p = params
        self.p = p
        self.motor_model = motor_model
        self.dt = p['dt']

        # Geometria
        self.A = p['A_beam']
        self.I = p['I_geom']
        self.H = p['H_geom']
        self.R = p['R_crank']
        self.unit_class = p.get('unit_class', 'class1')

        # Pre-computar cinematica (tabla de alta resolucion)
        self._theta_table = np.linspace(0, 2 * np.pi, 7200)

        if self.unit_class == 'mark2':
            self.P_rocker = p['P_rocker']
            self.C_pitman = p['C_pitman']
            y_PR_arr, TF_arr, _, _ = compute_kinematics_mark2(
                self._theta_table, self.A, self.P_rocker, self.C_pitman,
                self.I, self.H, self.R
            )
        else:
            self.C = p['C_beam']
            self.P = p['P_pitman']
            y_PR_arr, TF_arr, _, _ = compute_kinematics(
                self._theta_table, self.A, self.C, self.I, self.H, self.P, self.R
            )
        self._y_PR_table = y_PR_arr
        self._TF_table = TF_arr
        self._y_PR_eq = (np.max(y_PR_arr) + np.min(y_PR_arr)) / 2

        # Solver de onda
        g_eff = p['g'] * (1 - p['rho_f'] / p['rod_rho'])
        w_fl = p['rod_w'] * (1 - p['rho_f'] / p['rod_rho'])
        self.solver = WaveSolver(
            N=p['N_nodes'], L=p['rod_L'], E=p['rod_E'], rho=p['rod_rho'],
            A_r=p['rod_A'], c_damp=p['rod_c'], dt=p['dt'], g_eff=g_eff
        )
        x = self.solver.x
        u_s = w_fl / (p['rod_E'] * p['rod_A']) * (p['rod_L'] * x - x**2 / 2)
        self._u_static = u_s

        # Bomba
        self.pump = PumpModel(
            A_p=p['pump_A'], P_wh=p['P_wh'], P_ann=p['P_ann'],
            L_bomba=p['L_bomba'], rho_f=p['rho_f'], h_din=p['h_din_0'], g=p['g']
        )

        # Motor params
        self.omega_0 = p['omega_0']
        self.N_total = p['N_total']
        self.eta_trans = p['eta_trans']

    def _get_y_PR(self, theta):
        return np.interp(theta % (2 * np.pi), self._theta_table, self._y_PR_table)

    def _get_TF(self, theta):
        return np.interp(theta % (2 * np.pi), self._theta_table, self._TF_table)

    def _get_motor_torque(self, omega_m):
        p = self.p
        if self.motor_model == 'klauss':
            return torque_klauss(omega_m, p['omega_0'], p['T_max'], p['s_max'])
        else:
            return torque_linear(omega_m, p['omega_0'], p['k_m'], p['T_nom'], p['s_nom'])

    def initialize(self, theta_0=0.0, omega_0_crank=None):
        """Inicializa el estado del simulador."""
        self.theta = theta_0

        if omega_0_crank is None:
            # SPM nominal
            self.omega = 2 * np.pi * self.p['SPM_base'] / 60
        else:
            self.omega = omega_0_crank

        # Inicializar sarta en posicion coherente con theta
        y_0 = self._get_y_PR(theta_0)
        u_rigid = -(y_0 - self._y_PR_eq) * (1 - self.solver.x / self.p['rod_L'])
        self.solver.initialize(self._u_static + u_rigid)

    def step(self):
        """Avanza un paso de co-simulacion. Retorna dict con estado actual."""
        p = self.p
        theta = self.theta
        omega = self.omega
        dt = self.dt

        # 1. Cinematica actual
        y_PR = self._get_y_PR(theta)
        TF = self._get_TF(theta)

        # 2. BC superior para la PDE
        u_top = -(y_PR - self._y_PR_eq)

        # 3. Bomba
        v_plunger = self.solver.get_velocity_bottom()
        F_pump, q_inst = self.pump.update(v_plunger, dt)

        # 4. Avanzar PDE
        self.solver.step(u_top=u_top, F_bottom=F_pump)

        # 5. Fuerza en la PR
        F_PR = self.solver.get_force_top()

        # 6. Torque del motor
        omega_m = omega * p['N_total']
        T_motor = self._get_motor_torque(omega_m)
        T_motor_crank = T_motor * p['N_total'] * p['eta_trans']

        # 7. Inercia equivalente (simplificada — valor aproximado constante)
        J_ef = (p['J_c'] + p['M_cw'] * p['L_cw']**2
                + p['J_b'] * (TF * p['C_beam'] / p['A_beam'])**2 * 0.5
                + p['m_pit'] * p['R_crank']**2 * 0.5
                + p['m_PR'] * TF**2
                + p['J_m'] * p['N_total']**2)

        # dJ_ef/dtheta (derivada numerica)
        dtheta_num = 0.001
        TF_plus = self._get_TF(theta + dtheta_num)
        TF_minus = self._get_TF(theta - dtheta_num)
        J_plus = (p['J_c'] + p['M_cw'] * p['L_cw']**2
                  + p['J_b'] * (TF_plus * p['C_beam'] / p['A_beam'])**2 * 0.5
                  + p['m_pit'] * p['R_crank']**2 * 0.5
                  + p['m_PR'] * TF_plus**2
                  + p['J_m'] * p['N_total']**2)
        J_minus = (p['J_c'] + p['M_cw'] * p['L_cw']**2
                   + p['J_b'] * (TF_minus * p['C_beam'] / p['A_beam'])**2 * 0.5
                   + p['m_pit'] * p['R_crank']**2 * 0.5
                   + p['m_PR'] * TF_minus**2
                   + p['J_m'] * p['N_total']**2)
        dJ_dtheta = (J_plus - J_minus) / (2 * dtheta_num)

        # 8. Aceleracion angular
        T_cw = counterweight_torque(theta, p['M_cw'], p['L_cw'], p['g'], p['tau_cw'])
        T_load = F_PR * TF
        T_friction = p['c_friction'] * omega
        T_inertia_corr = 0.5 * dJ_dtheta * omega**2

        T_neto = T_motor_crank - T_load - T_cw - T_friction - T_inertia_corr
        alpha = T_neto / J_ef

        # 9. Integrar omega y theta (Euler semi-implicito)
        omega_new = omega + alpha * dt
        theta_new = theta + omega_new * dt

        # Actualizar estado
        self.theta = theta_new
        self.omega = omega_new

        # Torque en el reductor (lado lento)
        T_gearbox = T_load + T_cw + T_friction

        return {
            'theta': theta, 'omega': omega, 'alpha': alpha,
            'y_PR': y_PR, 'F_PR': F_PR, 'TF': TF,
            'T_motor': T_motor, 'T_motor_crank': T_motor_crank,
            'T_gearbox': T_gearbox, 'T_cw': T_cw, 'T_neto': T_neto,
            'J_ef': J_ef,
            'F_pump': F_pump, 'q_inst': q_inst, 'v_plunger': v_plunger,
            'omega_m': omega_m,
        }
