"""
Cargador de datos de pozo desde JSON → dict compatible con el Simulator.

Uso:
    params = load_well('M-912D-365-168.json', well_key='well_C147')
    sim = Simulator(params, motor_model='linear')
"""

import json
import numpy as np


def load_well(json_path, well_key='well_C147', stroke_index=1):
    """Carga un JSON de unidad de bombeo y retorna dict de parametros.

    Parameters
    ----------
    json_path : str — ruta al JSON de la unidad
    well_key : str — clave del pozo dentro del JSON
    stroke_index : int — configuracion de carrera (1, 2 o 3)

    Returns
    -------
    dict con todos los parametros necesarios para el Simulator
    """
    with open(json_path, 'r') as f:
        data = json.load(f)

    geom = data['geometry_inches']
    well = data[well_key]
    strokes = data['strokes']
    ratings = data['ratings']

    # Geometria Mark II (pulgadas → metros)
    IN2M = 0.0254
    A_beam = geom['A'] * IN2M
    C_pitman = geom['C'] * IN2M
    P_rocker = geom['P'] * IN2M
    I_geom = geom['I'] * IN2M
    H_abs = geom['H'] * IN2M
    G_height = geom['G'] * IN2M
    H_offset = H_abs - G_height

    stroke_key = f'stroke_{stroke_index}'
    R_crank = strokes[stroke_key]['R_inches'] * IN2M

    # Sarta de varillas
    FT2M = 0.3048
    LBF2N = 4.4482
    rod_sections = []
    for sec in well['rod_string']:
        d_m = sec['diameter_in'] * IN2M
        L_m = sec['count'] * sec['length_ft'] * FT2M
        A_m2 = np.pi / 4 * d_m**2
        w_Nm = sec['weight_lb_ft'] * LBF2N / FT2M
        rod_sections.append({
            'name': f"{sec['diameter_in']}\" {sec['material']}-grade",
            'd': d_m, 'A': A_m2, 'L': L_m,
            'E': 200.0e9, 'rho': 7850.0, 'w': w_Nm,
            'count': sec['count'],
        })

    totals = well['rod_string_totals']
    rod_L = totals['total_length_m']
    W_rod = totals['weight_in_air_lbs'] * LBF2N
    W_rod_fl = totals['weight_in_fluid_lbs'] * LBF2N

    # Bomba
    pump_d = well['pump']['diameter_in'] * IN2M
    pump_A = np.pi / 4 * pump_d**2
    L_bomba = well['pump']['depth_m']

    # Motor
    motor = well['motor']

    # Velocidad de onda
    rod_E = 200.0e9
    rod_rho = 7850.0
    rod_a = np.sqrt(rod_E / rod_rho)

    # Parametros numericos
    N_nodes = max(100, int(rod_L / 15))  # ~15m por nodo
    dx = rod_L / N_nodes
    dt = 0.98 * dx / rod_a  # alpha ≈ 0.98
    rod_nu = 0.10
    rod_c = np.pi * rod_a * rod_nu / (2 * rod_L)

    # Seccion inferior para solver uniforme
    last_sec = rod_sections[-1]

    params = {
        # Geometria Mark II
        'A_beam': A_beam,
        'P_rocker': P_rocker,
        'C_pitman': C_pitman,
        'C_beam': P_rocker,    # compatibilidad: C_beam = rocker en Mark II
        'P_pitman': C_pitman,  # compatibilidad: P_pitman = pitman en Mark II
        'I_geom': I_geom,
        'H_geom': H_offset,   # para la cinematica usa H_offset
        'H_abs': H_abs,
        'G_height': G_height,
        'H_offset': H_offset,
        'R_crank': R_crank,
        'S_stroke': well.get('measured_stroke_in', strokes[stroke_key]['stroke_inches']) * IN2M,

        # Tipo de unidad
        'unit_class': 'mark2',
        'unit_model': data['model'],

        # Motor
        'motor_HP': motor['power_hp'],
        'motor_RPM': motor['rpm'],

        # Sarta
        'rod_tapered': rod_sections,
        'rod_L': rod_L,
        'rod_E': rod_E,
        'rod_rho': rod_rho,
        'rod_a': rod_a,
        'rod_d': last_sec['d'],
        'rod_A': last_sec['A'],
        'rod_w': last_sec['w'],
        'rod_nu': rod_nu,
        'rod_c': rod_c,
        'W_rod': W_rod,
        'W_rod_fl': W_rod_fl,

        # Bomba
        'pump_d': pump_d,
        'pump_A': pump_A,
        'pump_eta_v': 0.80,
        'L_bomba': L_bomba,

        # Ratings
        'max_torque_in_lbs': ratings['max_torque_in_lbs'],
        'PPRL_rated_lbs': ratings['PPRL_lbs'],
        'catalog_TF': ratings['torque_factor'],
        'catalog_phase_angle': ratings['phase_angle_deg'],

        # Numerico
        'N_nodes': N_nodes,
        'dx': dx,
        'dt': dt,
        'alpha_CFL': rod_a * dt / dx,

        # Constantes
        'g': 9.81,
    }

    return params


def add_estimates(params, overrides=None):
    """Agrega parametros estimados (motor, inercias, contrapeso, yacimiento).

    Estos son los valores que se ajustan con el optimizador.

    Parameters
    ----------
    params : dict — de load_well()
    overrides : dict — valores para sobreescribir los defaults

    Returns
    -------
    dict completo listo para el Simulator
    """
    g = params['g']

    # Motor (50 Hz Argentina, 4 polos)
    f_red = 50.0
    n_polos = 4
    omega_0 = 2 * np.pi * f_red / (n_polos / 2)
    P_nom = params['motor_HP'] * 745.7
    s_nom = 0.08
    omega_nom = omega_0 * (1 - s_nom)
    T_nom = P_nom / omega_nom
    k_m = T_nom / (omega_0 * s_nom)

    defaults = {
        # Motor
        'f_red': f_red,
        'omega_0': omega_0,
        'P_nom': P_nom,
        's_nom': s_nom,
        'omega_nom': omega_nom,
        'T_nom': T_nom,
        'k_m': k_m,
        'T_max': 2.5 * T_nom,
        's_max': 0.25,
        'J_m': 1.5,

        # Transmision
        'N_pulley': 3.5,
        'N_gearbox': 40.0,
        'N_total': 140.0,
        'eta_trans': 0.90,

        # SPM
        'SPM_nominal': 6.0,
        'SPM_base': 6.0,

        # Inercias
        'J_c': 800.0,
        'J_b': 3000.0,
        'm_pit': 500.0,
        'm_PR': 700.0,

        # Contrapeso
        'M_cw': 5000.0,
        'L_cw': 1.20,
        'tau_cw': np.pi / 2,

        # Friccion
        'c_friction': 5.0,

        # Pozo
        'L_pozo': params['L_bomba'] + 60,
        'D_csg': 0.13970,
        'D_tbg_o': 0.07300,
        'D_tbg_i': 0.06200,
        'P_wh': 0.5e6,
        'P_ann': 0.2e6,
        'A_ann': np.pi / 4 * (0.13970**2 - 0.07300**2),

        # Fluido
        'WC': 0.90,
        'rho_f': 1008.0,
        'mu_f': 5.0e-3,
        'rho_oil': 947.0,
        'rho_water': 1015.0,

        # Yacimiento
        'P_r': 8.0e6,
        'q_max': 30.0 / 86400,
        'h_din_0': 900.0,
    }

    # Aplicar defaults
    for k, v in defaults.items():
        if k not in params:
            params[k] = v

    # Aplicar overrides
    if overrides:
        params.update(overrides)

    # Valores derivados
    params['F_fluid'] = params['rho_f'] * g * params['L_bomba'] * params['pump_A']
    params['PPRL_est'] = params['W_rod_fl'] + params['F_fluid']
    params['MPRL_est'] = params['W_rod_fl']

    return params
