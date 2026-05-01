"""
Modelo de bomba de fondo con valvulas TV y SV.

Implementa la logica de conmutacion de valvulas con histeresis
para evitar chattering, calcula presiones, caudal y carga sobre la sarta.

Referencia: AIB_modelo_simulacion.md, Seccion 3.6
"""

import numpy as np


class PumpModel:
    """Modelo de bomba de pistón con traveling valve (TV) y standing valve (SV).

    Estados:
    - UPSTROKE:   TV cerrada, SV abierta → F_pump = (P_d - P_s) * A_p
    - DOWNSTROKE: TV abierta, SV cerrada → F_pump = 0

    Parameters
    ----------
    A_p : float     — area del piston [m^2]
    P_wh : float    — presion en boca de pozo [Pa]
    P_ann : float   — presion en el anular [Pa]
    L_bomba : float — profundidad de la bomba [m]
    rho_f : float   — densidad del fluido [kg/m^3]
    h_din : float   — nivel dinamico inicial [m]
    g : float       — gravedad [m/s^2]
    v_threshold : float — umbral de velocidad para histeresis [m/s]
    """

    def __init__(self, A_p, P_wh, P_ann, L_bomba, rho_f, h_din, g=9.81,
                 v_threshold=0.001, v_smooth=0.01):
        self.A_p = A_p
        self.P_wh = P_wh
        self.P_ann = P_ann
        self.L_bomba = L_bomba
        self.rho_f = rho_f
        self.g = g
        self.v_threshold = v_threshold
        self.v_smooth = v_smooth  # ancho de transicion suave [m/s]

        # Estado inicial: downstroke (TV abierta)
        self.is_upstroke = False

        # Nivel dinamico (puede ser actualizado externamente)
        self.h_din = h_din

        # Acumuladores para caudal
        self.volume_cycle = 0.0
        self.volume_total = 0.0

    @property
    def P_d(self):
        """Presion de descarga [Pa]."""
        return self.P_wh + self.rho_f * self.g * self.L_bomba

    @property
    def P_s(self):
        """Presion de succion [Pa]."""
        return self.P_ann + self.rho_f * self.g * (self.L_bomba - self.h_din)

    @property
    def F_fluid(self):
        """Carga del fluido sobre el piston [N]."""
        return (self.P_d - self.P_s) * self.A_p

    def update(self, v_plunger, dt):
        """Actualiza estado de valvulas y calcula fuerza y caudal.

        Parameters
        ----------
        v_plunger : float — velocidad del embolo [m/s]
            Convencion: negativo = embolo sube (upstroke en coord. u)
        dt : float — paso de tiempo [s]

        Returns
        -------
        F_pump : float — fuerza en el fondo de la sarta [N] (traccion positiva)
        q_inst : float — caudal instantaneo [m^3/s]
        """
        # Transicion suave via sigmoid
        x = np.clip(-v_plunger / self.v_smooth, -20.0, 20.0)
        sigma = 1.0 / (1.0 + np.exp(-x))

        F_pump = self.F_fluid * sigma
        q_inst = self.A_p * max(0.0, -v_plunger) * sigma

        # Estado discreto (para diagnostico)
        self.is_upstroke = v_plunger < -self.v_threshold

        # Acumular volumen
        self.volume_cycle += q_inst * dt
        self.volume_total += q_inst * dt

        return F_pump, q_inst

    def reset_cycle_volume(self):
        """Resetea el acumulador de volumen del ciclo actual."""
        vol = self.volume_cycle
        self.volume_cycle = 0.0
        return vol
