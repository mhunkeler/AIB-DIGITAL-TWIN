"""
Sarta de varillas — Solver de la ecuacion de onda de Gibbs

Esquema leapfrog (diferencias finitas explicitas, orden 2) para:

    d2u/dt2 = a^2 * d2u/dx2 - c * du/dt + g_eff

con condiciones de borde:
  - Superior (Dirichlet): u(0,t) = u_top(t)
  - Inferior (Neumann):   E*A_r * du/dx|_L = F_bottom(t)

El termino g_eff (gravedad efectiva) permite trabajar con desplazamientos
totales, incluyendo el estiramiento por peso propio de la sarta.

Referencia: AIB_modelo_simulacion.md, Secciones 3.5 y 4.1
"""

import numpy as np


class WaveSolver:
    """Solver 1D de la ecuacion de onda con amortiguamiento (sarta uniforme).

    Parameters
    ----------
    N : int       — numero de nodos (0..N, total N+1 nodos)
    L : float     — longitud de la sarta [m]
    E : float     — modulo de Young [Pa]
    rho : float   — densidad [kg/m^3]
    A_r : float   — area transversal [m^2]
    c_damp : float — coeficiente de amortiguamiento [1/s]
    dt : float    — paso de tiempo [s]
    """

    def __init__(self, N, L, E, rho, A_r, c_damp, dt, g_eff=0.0):
        self.N = N
        self.L = L
        self.E = E
        self.rho = rho
        self.A_r = A_r
        self.c_damp = c_damp
        self.dt = dt
        self.g_eff = g_eff  # gravedad efectiva (con buoyancy) [m/s^2]

        self.dx = L / N
        self.a = np.sqrt(E / rho)  # velocidad de onda [m/s]
        self.alpha = self.a * dt / self.dx  # numero de Courant

        if self.alpha > 1.0:
            raise ValueError(
                f"CFL violada: alpha = {self.alpha:.4f} > 1.0. "
                f"Reducir dt o aumentar N."
            )

        # Arrays de estado: u_curr (n), u_prev (n-1), u_prev2 (n-2)
        self.u_curr = np.zeros(N + 1)
        self.u_prev = np.zeros(N + 1)
        self._u_bottom_prev2 = 0.0  # u[N] de 2 pasos atras (para vel. central)

        # Pre-calcular coeficientes
        self._alpha2 = self.alpha ** 2
        self._c_dt = c_damp * dt
        self._g_dt2 = g_eff * dt ** 2  # termino de gravedad

    def initialize(self, u0, v0=None):
        """Condicion inicial: u(x, 0) = u0, du/dt(x, 0) = v0.

        Parameters
        ----------
        u0 : ndarray (N+1,)  — desplazamiento inicial
        v0 : ndarray (N+1,) or None  — velocidad inicial (default: 0)
        """
        self.u_curr = u0.copy()
        if v0 is not None:
            self.u_prev = u0 - self.dt * v0
        else:
            self.u_prev = u0.copy()

    def step(self, u_top, F_bottom):
        """Avanza un paso de tiempo.

        Parameters
        ----------
        u_top : float     — desplazamiento impuesto en x=0 (BC Dirichlet)
        F_bottom : float  — fuerza aplicada en x=L (BC Neumann, positiva = traccion)

        Returns
        -------
        u_new : ndarray (N+1,) — desplazamiento en el nuevo paso
        """
        N = self.N
        a2 = self._alpha2
        c_dt = self._c_dt
        u = self.u_curr
        u_prev = self.u_prev

        u_new = np.empty(N + 1)

        # BC superior (Dirichlet)
        u_new[0] = u_top

        # Nodos interiores (leapfrog con amortiguamiento + gravedad)
        g_dt2 = self._g_dt2
        i = np.arange(1, N)
        u_new[i] = (
            2 * u[i] - u_prev[i]
            + a2 * (u[i + 1] - 2 * u[i] + u[i - 1])
            - c_dt * (u[i] - u_prev[i])
            + g_dt2
        )

        # BC inferior (Neumann) con nodo fantasma
        # du/dx|_L = F_bottom / (E * A_r)  (positivo = traccion)
        u_ghost = u[N - 1] + 2 * self.dx * F_bottom / (self.E * self.A_r)
        u_new[N] = (
            2 * u[N] - u_prev[N]
            + a2 * (u_ghost - 2 * u[N] + u[N - 1])
            - c_dt * (u[N] - u_prev[N])
            + g_dt2
        )

        # Actualizar estado
        self._u_bottom_prev2 = self.u_prev[N]  # guardar u[N] de 2 pasos atras
        self.u_prev = u.copy()
        self.u_curr = u_new

        return u_new

    def get_force_at(self, i):
        """Fuerza axial en el nodo i: F = E*A_r * du/dx."""
        if i == 0:
            # Forward difference O(h^2) con 3 puntos: (-3u0 + 4u1 - u2) / (2dx)
            dudx = (-3*self.u_curr[0] + 4*self.u_curr[1] - self.u_curr[2]) / (2 * self.dx)
        elif i == self.N:
            # Backward difference O(h^2) con 3 puntos
            dudx = (3*self.u_curr[self.N] - 4*self.u_curr[self.N-1] + self.u_curr[self.N-2]) / (2 * self.dx)
        else:
            # Central difference
            dudx = (self.u_curr[i + 1] - self.u_curr[i - 1]) / (2 * self.dx)
        return self.E * self.A_r * dudx

    def get_force_top(self):
        """Fuerza axial en x=0 (varilla pulida)."""
        return self.get_force_at(0)

    def get_velocity_bottom(self):
        """Velocidad en x=L (embolo de la bomba), diferencia backward."""
        return (self.u_curr[self.N] - self.u_prev[self.N]) / self.dt

    @property
    def x(self):
        """Posiciones de los nodos."""
        return np.linspace(0, self.L, self.N + 1)


class TaperedWaveSolver:
    """Solver de onda para sarta combinada (multi-tapered).

    Cada tramo tiene sus propias propiedades (E, rho, A_r, L).
    En las interfaces se impone continuidad de desplazamiento y fuerza.

    Parameters
    ----------
    sections : list of dict
        Cada dict tiene claves: 'L', 'E', 'rho', 'A', 'n_nodes'
    c_damp : float — coeficiente de amortiguamiento global [1/s]
    dt : float — paso de tiempo [s]
    """

    def __init__(self, sections, c_damp, dt):
        self.sections = sections
        self.c_damp = c_damp
        self.dt = dt
        self.n_sections = len(sections)

        # Construir malla global
        self._build_mesh()

        # Verificar CFL en todos los tramos
        for k, sec in enumerate(sections):
            a_k = np.sqrt(sec['E'] / sec['rho'])
            dx_k = sec['L'] / sec['n_nodes']
            alpha_k = a_k * dt / dx_k
            if alpha_k > 1.0:
                raise ValueError(
                    f"CFL violada en tramo {k}: alpha = {alpha_k:.4f}"
                )

        # Arrays de estado
        self.u_curr = np.zeros(self.n_total)
        self.u_prev = np.zeros(self.n_total)

    def _build_mesh(self):
        """Construye la malla global y mapeo de indices."""
        self.n_total = 0
        self.node_x = []
        self.node_a2 = []      # (a*dt/dx)^2 por nodo
        self.node_c_dt = []    # c*dt por nodo
        self.node_EA = []      # E*A por nodo
        self.node_dx = []      # dx por nodo
        self.interfaces = []   # indices globales de las interfaces

        x_offset = 0.0
        global_idx = 0

        for k, sec in enumerate(self.sections):
            n_k = sec['n_nodes']
            dx_k = sec['L'] / n_k
            a_k = np.sqrt(sec['E'] / sec['rho'])
            alpha_k = a_k * self.dt / dx_k
            EA_k = sec['E'] * sec['A']

            # Nodos de este tramo
            if k == 0:
                start = 0
            else:
                start = 1  # el primer nodo del tramo k coincide con el ultimo de k-1
                self.interfaces.append(global_idx - 1)

            for j in range(start, n_k + 1):
                self.node_x.append(x_offset + j * dx_k)
                self.node_a2.append(alpha_k ** 2)
                self.node_c_dt.append(self.c_damp * self.dt)
                self.node_EA.append(EA_k)
                self.node_dx.append(dx_k)
                global_idx += 1

            x_offset += sec['L']

        self.n_total = global_idx
        self.node_x = np.array(self.node_x)
        self.node_a2 = np.array(self.node_a2)
        self.node_c_dt = np.array(self.node_c_dt)
        self.node_EA = np.array(self.node_EA)
        self.node_dx = np.array(self.node_dx)
        self.L = x_offset

        # Para la BC inferior
        last_sec = self.sections[-1]
        self._bottom_E = last_sec['E']
        self._bottom_A = last_sec['A']
        self._bottom_dx = last_sec['L'] / last_sec['n_nodes']

    def initialize(self, u0, v0=None):
        """Condicion inicial."""
        self.u_curr = u0.copy()
        if v0 is not None:
            self.u_prev = u0 - self.dt * v0
        else:
            self.u_prev = u0.copy()

    def step(self, u_top, F_bottom):
        """Avanza un paso de tiempo."""
        N = self.n_total - 1
        u = self.u_curr
        u_prev = self.u_prev
        a2 = self.node_a2
        c_dt = self.node_c_dt
        EA = self.node_EA
        dx = self.node_dx

        u_new = np.empty(self.n_total)

        # BC superior
        u_new[0] = u_top

        # Nodos interiores (excluyendo interfaces)
        interface_set = set(self.interfaces)
        for i in range(1, N):
            if i in interface_set:
                # Interfaz: ecuacion especial con continuidad de fuerza
                # E_L * A_L * (u[i] - u[i-1]) / dx_L = E_R * A_R * (u[i+1] - u[i]) / dx_R
                # Usar promedio ponderado
                EA_L = self.node_EA[i - 1] if i > 0 else EA[i]
                EA_R = self.node_EA[i + 1] if i < N else EA[i]
                dx_L = dx[i]
                dx_R = dx[i + 1] if i < N else dx[i]

                # Avance con continuidad de fuerza:
                # Fuerza izquierda: EA_L * (u[i] - u[i-1]) / dx_L
                # Fuerza derecha:   EA_R * (u[i+1] - u[i]) / dx_R
                # Aceleracion del nodo masa promedio
                rho_avg = 0.5 * (self.sections[0]['rho'] + self.sections[-1]['rho'])
                A_avg = 0.5 * (EA_L / self.sections[0]['E'] + EA_R / self.sections[-1]['E'])
                mass_node = rho_avg * A_avg * 0.5 * (dx_L + dx_R)

                F_L = EA_L * (u[i - 1] - u[i]) / dx_L
                F_R = EA_R * (u[i + 1] - u[i]) / dx_R
                accel = (F_L + F_R) / mass_node

                u_new[i] = (
                    2 * u[i] - u_prev[i]
                    + self.dt ** 2 * accel
                    - c_dt[i] * (u[i] - u_prev[i])
                )
            else:
                u_new[i] = (
                    2 * u[i] - u_prev[i]
                    + a2[i] * (u[i + 1] - 2 * u[i] + u[i - 1])
                    - c_dt[i] * (u[i] - u_prev[i])
                )

        # BC inferior (Neumann)
        u_ghost = u[N - 1] - 2 * self._bottom_dx * F_bottom / (self._bottom_E * self._bottom_A)
        u_new[N] = (
            2 * u[N] - u_prev[N]
            + a2[N] * (u_ghost - 2 * u[N] + u[N - 1])
            - c_dt[N] * (u[N] - u_prev[N])
        )

        self.u_prev = u.copy()
        self.u_curr = u_new
        return u_new

    def get_force_top(self):
        """Fuerza axial en x=0."""
        dudx = (self.u_curr[1] - self.u_curr[0]) / self.node_dx[0]
        return self.node_EA[0] * dudx / self.sections[0]['E'] * self.sections[0]['E']

    def get_velocity_bottom(self):
        """Velocidad en el ultimo nodo."""
        return (self.u_curr[-1] - self.u_prev[-1]) / self.dt

    @property
    def x(self):
        return self.node_x
