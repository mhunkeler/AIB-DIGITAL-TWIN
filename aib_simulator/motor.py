"""
Motor de induccion — modelos de torque-velocidad.

Modelo lineal: T_m = T_0 - k_m * (omega_m - omega_0)
Modelo Klauss: T_m = 2*T_max / (s/s_max + s_max/s)

Referencia: AIB_modelo_simulacion.md, Seccion 3.3
"""

import numpy as np


def torque_linear(omega_m, omega_0, k_m, T_nom, s_nom):
    """Torque del motor — modelo lineal.

    T_m(omega_m) = T_0 - k_m * omega_m
    donde T_0 = T_nom + k_m * omega_nom

    Parameters
    ----------
    omega_m : float — velocidad angular del motor [rad/s]
    omega_0 : float — velocidad sincronica [rad/s]
    k_m : float — pendiente [N*m*s/rad]
    T_nom : float — torque nominal [N*m]
    s_nom : float — slip nominal

    Returns
    -------
    T_m : float — torque del motor [N*m]
    """
    omega_nom = omega_0 * (1 - s_nom)
    T_0 = T_nom + k_m * omega_nom
    T_m = T_0 - k_m * omega_m
    return T_m  # permite torque negativo (frenado regenerativo)


def torque_klauss(omega_m, omega_0, T_max, s_max):
    """Torque del motor — modelo Klauss.

    T_m(s) = 2*T_max / (s/s_max + s_max/s)

    Parameters
    ----------
    omega_m : float — velocidad angular del motor [rad/s]
    omega_0 : float — velocidad sincronica [rad/s]
    T_max : float — torque pull-out [N*m]
    s_max : float — slip al pull-out

    Returns
    -------
    T_m : float — torque del motor [N*m]
    """
    s = (omega_0 - omega_m) / omega_0
    if abs(s) < 0.001:
        s = 0.001 if s >= 0 else -0.001
    T_m = 2 * T_max / (s / s_max + s_max / s)
    return T_m  # permite torque negativo (frenado regenerativo super-sincrono)
