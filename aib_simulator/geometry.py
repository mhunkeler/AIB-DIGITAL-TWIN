"""
Cinematica del balancin — Clase I convencional (API Spec 11E)

Resuelve el cuadrilatero articulado ciguenal-biela-balancin para obtener
la posicion, velocidad, aceleracion y factor de torque de la varilla pulida
en funcion del angulo del ciguenal theta.

Referencia: AIB_modelo_simulacion.md, Seccion 3.2
"""

import numpy as np


def crank_pin(theta, I, H, R):
    """Posicion del perno de la manivela (crank pin) dado theta.

    Sistema de coordenadas: saddle bearing en el origen,
    centro del ciguenal en (I, -H).

    Parameters
    ----------
    theta : float or ndarray
        Angulo del ciguenal [rad], medido desde eje X positivo.
    I : float  — offset horizontal SB → ciguenal [m]
    H : float  — offset vertical SB → ciguenal [m]
    R : float  — radio de la manivela [m]

    Returns
    -------
    x_cp, y_cp : float or ndarray
    """
    x_cp = I + R * np.cos(theta)
    y_cp = -H + R * np.sin(theta)
    return x_cp, y_cp


def polished_rod_position(theta, A, C, I, H, P, R):
    """Posicion vertical de la varilla pulida y_PR para un angulo theta.

    Parameters
    ----------
    theta : float or ndarray  — angulo del ciguenal [rad]
    A : float  — brazo delantero del balancin [m]
    C : float  — brazo trasero del balancin [m]
    I : float  — offset horizontal SB → ciguenal [m]
    H : float  — offset vertical SB → ciguenal [m]
    P : float  — longitud de la biela (pitman) [m]
    R : float  — radio de la manivela [m]

    Returns
    -------
    y_PR : float or ndarray  — posicion vertical de la varilla pulida [m]
    """
    x_cp, y_cp = crank_pin(theta, I, H, R)

    # Distancia SB → crank pin
    J = np.sqrt(x_cp**2 + y_cp**2)

    # Angulo en el SB (ley de cosenos en el triangulo SB-EB-CP)
    cos_psi = (C**2 + J**2 - P**2) / (2 * C * J)
    cos_psi = np.clip(cos_psi, -1.0, 1.0)
    psi_b = np.arccos(cos_psi)

    # Angulo absoluto SB → CP
    phi = np.arctan2(y_cp, x_cp)

    # Posicion del equalizer bearing
    y_E = C * np.sin(phi + psi_b)

    # Posicion de la varilla pulida (Clase I: horsehead en lado opuesto)
    y_PR = -(A / C) * y_E

    return y_PR


def compute_kinematics(theta_array, A, C, I, H, P, R):
    """Calcula posicion, factor de torque y datos del equalizer para un array de theta.

    Returns
    -------
    y_PR : ndarray  — posicion de la varilla pulida [m]
    TF : ndarray    — factor de torque dy_PR/dtheta [m/rad]
    y_E : ndarray   — posicion vertical del equalizer bearing [m]
    x_E : ndarray   — posicion horizontal del equalizer bearing [m]
    """
    x_cp, y_cp = crank_pin(theta_array, I, H, R)
    J = np.sqrt(x_cp**2 + y_cp**2)

    cos_psi = (C**2 + J**2 - P**2) / (2 * C * J)
    cos_psi = np.clip(cos_psi, -1.0, 1.0)
    psi_b = np.arccos(cos_psi)

    phi = np.arctan2(y_cp, x_cp)
    angle = phi + psi_b

    x_E = C * np.cos(angle)
    y_E = C * np.sin(angle)
    y_PR = -(A / C) * y_E

    # Factor de torque: derivada numerica dy_PR/dtheta
    dtheta = theta_array[1] - theta_array[0]
    TF = np.gradient(y_PR, dtheta)

    return y_PR, TF, y_E, x_E


def full_kinematics(theta_array, omega, A, C, I, H, P, R):
    """Cinematica completa: posicion, velocidad, aceleracion, TF.

    Parameters
    ----------
    theta_array : ndarray  — angulos del ciguenal [rad]
    omega : float          — velocidad angular del ciguenal [rad/s]
    A, C, I, H, P, R : float  — geometria del mecanismo [m]

    Returns
    -------
    dict con claves: y_PR, v_PR, a_PR, TF, x_E, y_E
    """
    y_PR, TF, y_E, x_E = compute_kinematics(theta_array, A, C, I, H, P, R)

    dtheta = theta_array[1] - theta_array[0]

    # Velocidad: v_PR = (dy_PR/dtheta) * omega = TF * omega
    v_PR = TF * omega

    # Aceleracion: a_PR = (d2y_PR/dtheta2) * omega^2 + (dy_PR/dtheta) * alpha
    # Con omega constante (alpha = 0): a_PR = d2y_PR/dtheta2 * omega^2
    d2y_dtheta2 = np.gradient(TF, dtheta)
    a_PR = d2y_dtheta2 * omega**2

    return {
        'y_PR': y_PR,
        'v_PR': v_PR,
        'a_PR': a_PR,
        'TF': TF,
        'd2y_dtheta2': d2y_dtheta2,
        'x_E': x_E,
        'y_E': y_E,
    }


def mechanism_positions(theta, A, C, I, H, P, R):
    """Posiciones de todos los puntos del mecanismo para visualizacion.

    Returns
    -------
    dict con las coordenadas de cada punto clave del mecanismo.
    """
    # Centro del ciguenal
    x_crank_center = I
    y_crank_center = -H

    # Crank pin
    x_cp, y_cp = crank_pin(theta, I, H, R)

    # Saddle bearing (origen)
    x_sb, y_sb = 0.0, 0.0

    # Equalizer bearing
    J = np.sqrt(x_cp**2 + y_cp**2)
    cos_psi = np.clip((C**2 + J**2 - P**2) / (2 * C * J), -1.0, 1.0)
    psi_b = np.arccos(cos_psi)
    phi = np.arctan2(y_cp, x_cp)
    angle = phi + psi_b

    x_eb = C * np.cos(angle)
    y_eb = C * np.sin(angle)

    # Horsehead (extension del balancin al lado opuesto del equalizer)
    # Direccion SB → EB, luego extender al lado opuesto por distancia A
    beam_angle = np.arctan2(y_eb, x_eb)
    x_hh = -A * np.cos(beam_angle)
    y_hh = -A * np.sin(beam_angle)

    # Posicion de la varilla pulida (punto mas bajo del horsehead)
    y_PR = -(A / C) * y_eb

    return {
        'crank_center': (x_crank_center, y_crank_center),
        'crank_pin': (x_cp, y_cp),
        'saddle_bearing': (x_sb, y_sb),
        'equalizer': (x_eb, y_eb),
        'horsehead': (x_hh, y_hh),
        'polished_rod': (x_hh, y_PR),
    }
