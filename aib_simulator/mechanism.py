"""
Dinamica del mecanismo — EOM Lagrange para el ciguenal.

J_ef(theta) * ddtheta + 0.5 * dJ_ef/dtheta * dtheta^2 = T_neto

Donde T_neto = T_motor*N*eta - T_cw - F_PR*TF - T_friction

Referencia: AIB_modelo_simulacion.md, Seccion 3.4
"""

import numpy as np


def equivalent_inertia(theta, J_c, J_b, m_pit, m_PR, A, C, P, R, I_geom, H,
                       M_cw, L_cw, J_m, N_total):
    """Inercia equivalente referida al ciguenal.

    Incluye contribuciones de:
    - Ciguenal + manivelas (J_c)
    - Balancin (J_b) proyectado via TF
    - Biela (m_pit * R^2 aproximado)
    - Cabezal + varilla pulida (m_PR * TF^2)
    - Contrapesos (M_cw * L_cw^2)
    - Motor (J_m * N^2) — inercia reflejada

    Parameters
    ----------
    theta : float — angulo del ciguenal [rad]
    (otros: parametros geometricos e inerciales)

    Returns
    -------
    J_ef : float — inercia equivalente [kg*m^2]
    """
    from aib_simulator.geometry import compute_kinematics

    # TF en este angulo
    theta_arr = np.array([theta - 0.001, theta, theta + 0.001])
    _, TF_arr, _, _ = compute_kinematics(theta_arr, A, C, I_geom, H, P, R)
    TF = TF_arr[1]

    # Contribuciones
    J_crank = J_c + M_cw * L_cw ** 2
    J_beam_projected = J_b * (TF * C / A) ** 2  # proyeccion del balancin
    J_pitman = m_pit * R ** 2 * 0.5  # aproximacion
    J_PR_projected = m_PR * TF ** 2
    J_motor_reflected = J_m * N_total ** 2

    J_ef = J_crank + J_beam_projected + J_pitman + J_PR_projected + J_motor_reflected
    return J_ef


def counterweight_torque(theta, M_cw, L_cw, g, tau_cw):
    """Torque del contrapeso sobre el ciguenal.

    T_cw = M_cw * g * L_cw * sin(theta + tau_cw)

    Parameters
    ----------
    theta : float — angulo del ciguenal [rad]
    M_cw : float — masa del contrapeso [kg]
    L_cw : float — brazo efectivo [m]
    g : float — gravedad [m/s^2]
    tau_cw : float — fase del contrapeso [rad] (180 deg = opuesto a la manivela)

    Returns
    -------
    T_cw : float — torque del contrapeso [N*m]
    """
    return M_cw * g * L_cw * np.sin(theta + tau_cw)


def compute_acceleration(theta, omega, F_PR, TF, T_motor_crank,
                         J_ef, dJ_ef_dtheta, M_cw, L_cw, g_val, tau_cw,
                         c_friction):
    """Calcula la aceleracion angular del ciguenal.

    J_ef * alpha = T_neto - 0.5 * dJ_ef/dtheta * omega^2

    Parameters
    ----------
    theta : float — angulo [rad]
    omega : float — velocidad angular [rad/s]
    F_PR : float — fuerza en la varilla pulida [N]
    TF : float — factor de torque [m/rad]
    T_motor_crank : float — torque del motor referido al ciguenal [N*m]
    J_ef : float — inercia equivalente [kg*m^2]
    dJ_ef_dtheta : float — derivada de J_ef respecto a theta
    M_cw, L_cw, g_val, tau_cw : contrapeso
    c_friction : float — coeficiente de friccion viscosa [N*m*s/rad]

    Returns
    -------
    alpha : float — aceleracion angular [rad/s^2]
    """
    T_cw = counterweight_torque(theta, M_cw, L_cw, g_val, tau_cw)
    T_load = F_PR * TF  # torque de la carga sobre el ciguenal
    T_friction = c_friction * omega
    T_inertia_correction = 0.5 * dJ_ef_dtheta * omega ** 2

    T_neto = T_motor_crank - T_load - T_cw - T_friction - T_inertia_correction

    alpha = T_neto / J_ef
    return alpha
