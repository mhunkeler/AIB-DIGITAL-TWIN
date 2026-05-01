"""
Modelo de yacimiento — IPR Vogel + balance de nivel dinamico.

El caudal de entrada al pozo se calcula con la correlacion de Vogel
para pozos saturados. El nivel dinamico evoluciona segun el balance
volumetrico en el espacio anular.

Referencia: AIB_modelo_simulacion.md, Seccion 3.7
"""

import numpy as np


def vogel_ipr(P_wf, P_r, q_max):
    """Caudal de entrada del yacimiento (IPR Vogel).

    q_in/q_max = 1 - 0.2*(P_wf/P_r) - 0.8*(P_wf/P_r)^2

    Parameters
    ----------
    P_wf : float — presion de fondo fluyente [Pa]
    P_r : float  — presion estatica del yacimiento [Pa]
    q_max : float — AOFP (caudal absoluto de pozo abierto) [m^3/s]

    Returns
    -------
    q_in : float — caudal de entrada [m^3/s]
    """
    if P_wf <= 0:
        return q_max
    if P_wf >= P_r:
        return 0.0
    ratio = P_wf / P_r
    q_in = q_max * (1 - 0.2 * ratio - 0.8 * ratio ** 2)
    return max(0.0, q_in)


def bottom_hole_pressure(h_din, P_ann, rho_f, g, L_bomba):
    """Presion de fondo fluyente a partir del nivel dinamico.

    P_wf = P_ann + rho_f * g * (L_bomba - h_din)

    Parameters
    ----------
    h_din : float  — nivel dinamico (profundidad desde superficie) [m]
    P_ann : float  — presion en el anular [Pa]
    rho_f : float  — densidad del fluido [kg/m^3]
    g : float      — gravedad [m/s^2]
    L_bomba : float — profundidad de la bomba [m]

    Returns
    -------
    P_wf : float — presion de fondo fluyente [Pa]
    """
    submergence = L_bomba - h_din  # columna de fluido sobre la bomba
    if submergence < 0:
        submergence = 0  # bomba por encima del nivel (seca)
    return P_ann + rho_f * g * submergence


class ReservoirModel:
    """Modelo de yacimiento con balance de nivel dinamico.

    dh_din/dt = (q_pump - q_in) / A_ann

    Parameters
    ----------
    P_r : float    — presion del yacimiento [Pa]
    q_max : float  — AOFP [m^3/s]
    P_ann : float  — presion anular [Pa]
    rho_f : float  — densidad del fluido [kg/m^3]
    L_bomba : float — profundidad de la bomba [m]
    A_ann : float  — area del anular [m^2]
    h_din_0 : float — nivel dinamico inicial [m]
    g : float      — gravedad [m/s^2]
    """

    def __init__(self, P_r, q_max, P_ann, rho_f, L_bomba, A_ann, h_din_0,
                 g=9.81):
        self.P_r = P_r
        self.q_max = q_max
        self.P_ann = P_ann
        self.rho_f = rho_f
        self.L_bomba = L_bomba
        self.A_ann = A_ann
        self.g = g
        self.h_din = h_din_0

    @property
    def P_wf(self):
        """Presion de fondo fluyente actual [Pa]."""
        return bottom_hole_pressure(
            self.h_din, self.P_ann, self.rho_f, self.g, self.L_bomba
        )

    @property
    def q_in(self):
        """Caudal de entrada del yacimiento [m^3/s]."""
        return vogel_ipr(self.P_wf, self.P_r, self.q_max)

    def update(self, q_pump_avg, dt_slow):
        """Actualiza el nivel dinamico (paso lento).

        Parameters
        ----------
        q_pump_avg : float — caudal promedio de la bomba en el periodo [m^3/s]
        dt_slow : float — paso de tiempo lento [s]

        Returns
        -------
        dict con h_din, P_wf, q_in, q_pump_avg
        """
        q_in = self.q_in

        # Balance volumetrico: si la bomba extrae mas de lo que entra,
        # el nivel baja (h_din aumenta = mas profundo)
        dh_dt = (q_pump_avg - q_in) / self.A_ann
        self.h_din += dh_dt * dt_slow

        # Limites fisicos
        self.h_din = np.clip(self.h_din, 0.0, self.L_bomba)

        return {
            'h_din': self.h_din,
            'P_wf': self.P_wf,
            'q_in': self.q_in,
            'q_pump_avg': q_pump_avg,
        }
