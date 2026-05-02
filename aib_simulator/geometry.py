"""
Cinematica del balancin — Clase I y Mark II (API Spec 11E)

Clase I (convencional): pitman conecta DETRAS del SB, horsehead en lado opuesto.
  y_PR = -(A/C) * y_E

Mark II (Clase III / Unitorque): pitman conecta DELANTE del SB via equalizer,
  horsehead en el MISMO lado que el equalizer.
  y_PR = +(A/P_rocker) * y_E

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


# =========================================================================
# MARK II (Clase III / Unitorque) — cinematica
# =========================================================================

def polished_rod_position_mark2(theta, A, P_rocker, C_pitman, I, H_offset, R):
    """Posicion de la varilla pulida para unidad Mark II (Clase III).

    En Mark II el cuadrilatero es: ground → crank(R) → pitman(C) → rocker(P).
    El horsehead esta en el MISMO lado que el equalizer, a distancia A del SB.

    Parameters
    ----------
    theta : float or ndarray — angulo del ciguenal [rad]
    A : float         — SB → horsehead [m]
    P_rocker : float  — SB → equalizer bearing (rocker del 4-bar) [m]
    C_pitman : float  — longitud de la biela [m]
    I : float         — offset horizontal SB → ciguenal [m]
    H_offset : float  — offset vertical SB → ciguenal (= H_abs - G) [m]
    R : float         — radio de la manivela [m]
    """
    x_cp = I + R * np.cos(theta)
    y_cp = -H_offset + R * np.sin(theta)

    J = np.sqrt(x_cp**2 + y_cp**2)

    # Ley de cosenos: triangulo SB(origen) - EB - CP
    # Lados: P_rocker (SB→EB), C_pitman (EB→CP), J (SB→CP)
    cos_psi = (P_rocker**2 + J**2 - C_pitman**2) / (2 * P_rocker * J)
    cos_psi = np.clip(cos_psi, -1.0, 1.0)
    psi = np.arccos(cos_psi)

    phi = np.arctan2(y_cp, x_cp)

    # Equalizer bearing: a distancia P_rocker del SB
    y_E = P_rocker * np.sin(phi + psi)

    # Mark II: el ratio de amplificacion es C_pitman/P_rocker (no A/P_rocker).
    # Esto se verifica contra las 3 carreras del catalogo con error < 1%.
    # El ratio A/P sobreestima la carrera en ~23%.
    y_PR = (C_pitman / P_rocker) * y_E

    return y_PR


def compute_kinematics_mark2(theta_array, A, P_rocker, C_pitman, I, H_offset, R):
    """Cinematica Mark II: posicion, TF, equalizer."""
    x_cp = I + R * np.cos(theta_array)
    y_cp = -H_offset + R * np.sin(theta_array)
    J = np.sqrt(x_cp**2 + y_cp**2)

    cos_psi = (P_rocker**2 + J**2 - C_pitman**2) / (2 * P_rocker * J)
    cos_psi = np.clip(cos_psi, -1.0, 1.0)
    psi = np.arccos(cos_psi)
    phi = np.arctan2(y_cp, x_cp)
    angle = phi + psi

    x_E = P_rocker * np.cos(angle)
    y_E = P_rocker * np.sin(angle)
    y_PR = (C_pitman / P_rocker) * y_E

    dtheta = theta_array[1] - theta_array[0]
    TF = np.gradient(y_PR, dtheta)

    return y_PR, TF, y_E, x_E


def full_kinematics_mark2(theta_array, omega, A, P_rocker, C_pitman, I, H_offset, R):
    """Cinematica completa Mark II: posicion, velocidad, aceleracion, TF."""
    y_PR, TF, y_E, x_E = compute_kinematics_mark2(
        theta_array, A, P_rocker, C_pitman, I, H_offset, R
    )
    dtheta = theta_array[1] - theta_array[0]
    v_PR = TF * omega
    d2y_dtheta2 = np.gradient(TF, dtheta)
    a_PR = d2y_dtheta2 * omega**2

    return {
        'y_PR': y_PR, 'v_PR': v_PR, 'a_PR': a_PR, 'TF': TF,
        'd2y_dtheta2': d2y_dtheta2, 'x_E': x_E, 'y_E': y_E,
    }


def mechanism_positions_mark2(theta, A, P_rocker, C_pitman, I, H_offset, R, G=0.0):
    """Posiciones del mecanismo Mark II para visualizacion.

    Parameters
    ----------
    G : float — altura del ciguenal sobre la base [m] (para dibujar el suelo)
    """
    x_crank_center = I
    y_crank_center = -H_offset

    x_cp = I + R * np.cos(theta)
    y_cp = -H_offset + R * np.sin(theta)

    x_sb, y_sb = 0.0, 0.0

    J = np.sqrt(x_cp**2 + y_cp**2)
    cos_psi = np.clip((P_rocker**2 + J**2 - C_pitman**2) / (2 * P_rocker * J), -1.0, 1.0)
    psi = np.arccos(cos_psi)
    phi = np.arctan2(y_cp, x_cp)
    angle = phi + psi

    # Equalizer bearing (a distancia P del SB, mismo lado que horsehead)
    x_eb = P_rocker * np.cos(angle)
    y_eb = P_rocker * np.sin(angle)

    # Horsehead (extension en la MISMA direccion, a distancia A del SB)
    x_hh = A * np.cos(angle)
    y_hh = A * np.sin(angle)

    y_PR = (C_pitman / P_rocker) * y_eb

    return {
        'crank_center': (x_crank_center, y_crank_center),
        'crank_pin': (x_cp, y_cp),
        'saddle_bearing': (x_sb, y_sb),
        'equalizer': (x_eb, y_eb),
        'horsehead': (x_hh, y_hh),
        'polished_rod': (x_hh, y_PR),
        'base_height': -(H_offset + G),  # nivel del suelo
    }
