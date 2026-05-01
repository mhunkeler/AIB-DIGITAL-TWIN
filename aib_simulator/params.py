"""
Parametros del caso de prueba — Seccion 6 de AIB_modelo_simulacion.md

Todos los valores en unidades SI (m, kg, s, N, Pa, rad).
Caso: pozo tipico Cuenca del Golfo San Jorge, unidad C-228D-200-86.
"""

import numpy as np

# ---------------------------------------------------------------------------
# 6.1 Pozo
# ---------------------------------------------------------------------------
L_bomba = 1100.0        # Profundidad de la bomba [m]
L_pozo = 1200.0         # Profundidad total del pozo [m]
D_csg = 0.13970         # Diametro interior del casing 5-1/2" [m]
D_tbg_o = 0.07300       # Diametro exterior del tubing 2-7/8" [m]
D_tbg_i = 0.06200       # Diametro interior del tubing [m]
P_wh = 0.5e6            # Presion en boca de pozo (manometrica) [Pa]
P_ann = 0.2e6           # Presion en anular (manometrica) [Pa]

# Area anular (entre casing y tubing exterior)
A_ann = np.pi / 4 * (D_csg**2 - D_tbg_o**2)  # [m^2]

# ---------------------------------------------------------------------------
# 6.2 Unidad de bombeo — Clase I convencional, C-228D-200-86
# ---------------------------------------------------------------------------
# Geometria del mecanismo (API Spec 11E)
A_beam = 4.000          # Brazo delantero del balancin (SB → horsehead) [m]
C_beam = 2.500          # Brazo trasero del balancin (SB → equalizer) [m]
I_geom = 2.000          # Offset horizontal SB → ciguenal [m]
H_geom = 4.500          # Altura vertical SB → ciguenal [m]
P_pitman = 3.500        # Longitud de la biela (pitman) [m]
R_crank = 0.900         # Radio de la manivela [m]
S_stroke = 2.88         # Carrera teorica resultante [m] (aprox 2R*A/C)

SPM_nominal = 8.0       # SPM nominal de operacion [ciclos/min]
SPM_base = 6.0          # SPM caso base (evita fluid pound) [ciclos/min]

# Inercias
J_c = 200.0             # Ciguenal + manivelas (sin contrapesos) [kg*m^2]
J_b = 800.0             # Balancin alrededor del SB [kg*m^2]
m_pit = 250.0           # Masa total de la biela [kg]
m_PR = 350.0            # Masa cabezal + varilla pulida + adaptadores [kg]

# Contrabalanceo (sobre la manivela)
M_cw = 2500.0           # Masa del contrapeso total [kg]
L_cw = 0.85             # Brazo efectivo del contrapeso [m]
tau_cw = np.pi          # Fase respecto a la manivela: 180 deg [rad]

# ---------------------------------------------------------------------------
# 6.3 Motor + transmision (NEMA D, 50 Hz, Argentina)
# ---------------------------------------------------------------------------
P_nom = 30.0e3          # Potencia nominal [W] (40 HP)
f_red = 50.0            # Frecuencia de red [Hz]
n_polos = 4             # Numero de polos
omega_0 = 2 * np.pi * f_red / (n_polos / 2)  # Vel. sincronica [rad/s] = 157.08
s_nom = 0.12            # Slip nominal (NEMA D)
omega_nom = omega_0 * (1 - s_nom)             # Vel. nominal [rad/s] = 138.23
T_nom = P_nom / omega_nom                     # Torque nominal [N*m] = ~217
k_m = T_nom / (omega_0 - omega_nom)           # Pendiente lineal [N*m*s/rad] = ~11.5

# Modelo Klauss
T_max = 2.5 * T_nom     # Torque pull-out [N*m] = ~540
s_max = 0.40            # Slip al pull-out

# Inercia del rotor
J_m = 0.45              # [kg*m^2]

# Transmision
N_pulley = 5.0          # Relacion poleas (motor → reductor rapido)
N_gearbox = 30.0        # Relacion caja reductora
N_total = N_pulley * N_gearbox  # Relacion total = 150
eta_trans = 0.92        # Eficiencia de la transmision

# Friccion del mecanismo (modelo viscoso simple)
c_friction = 5.0        # Coeficiente de friccion viscosa [N*m*s/rad] (estimado)

# ---------------------------------------------------------------------------
# 6.4 Sarta de varillas
# ---------------------------------------------------------------------------

# --- Configuracion base: sarta uniforme 3/4" D-grade ---
rod_L = 1100.0          # Longitud total [m]
rod_d = 0.01905         # Diametro nominal 3/4" [m]
rod_A = 285.0e-6        # Area transversal [m^2]
rod_w = 23.5            # Peso lineal en aire (con couplings) [N/m]
rod_E = 200.0e9         # Modulo de Young [Pa]
rod_rho = 7850.0        # Densidad del acero [kg/m^3]
rod_a = 5050.0          # Velocidad de onda [m/s]
rod_nu = 0.10           # Factor de amortiguamiento Gibbs (adimensional)
rod_c = np.pi * rod_a * rod_nu / (2 * rod_L)  # Coef. amortiguamiento [1/s] = 0.722

# Peso total de la sarta
W_rod = rod_w * rod_L   # Peso en aire [N] = 25850

# --- Configuracion avanzada: sarta combinada 7/8" + 3/4" ---
rod_tapered = [
    {
        'name': 'Superior 7/8"',
        'd': 0.02223,          # 7/8" [m]
        'A': 388.0e-6,         # Area transversal [m^2]
        'L': 550.0,            # Longitud [m]
        'E': 200.0e9,          # Modulo de Young [Pa]
        'rho': 7850.0,         # Densidad [kg/m^3]
        'w': 30.7,             # Peso lineal [N/m] (estimado para 7/8")
    },
    {
        'name': 'Inferior 3/4"',
        'd': 0.01905,          # 3/4" [m]
        'A': 285.0e-6,         # Area transversal [m^2]
        'L': 550.0,            # Longitud [m]
        'E': 200.0e9,          # Modulo de Young [Pa]
        'rho': 7850.0,         # Densidad [kg/m^3]
        'w': 23.5,             # Peso lineal [N/m]
    },
]

# ---------------------------------------------------------------------------
# 6.5 Bomba
# ---------------------------------------------------------------------------
pump_d = 0.04445         # Diametro del piston 1-3/4" [m]
pump_A = 1552.0e-6       # Area del piston [m^2]
pump_eta_v = 0.80        # Eficiencia volumetrica esperada

# ---------------------------------------------------------------------------
# 6.6 Fluido y yacimiento
# ---------------------------------------------------------------------------
# Fluido producido (mezcla agua + petroleo, alto WC)
WC = 0.90               # Corte de agua [fraccion]
rho_oil = 947.0          # Densidad del petroleo 18 API [kg/m^3]
rho_water = 1015.0       # Densidad del agua de formacion [kg/m^3]
rho_f = 1008.0           # Densidad efectiva mezcla [kg/m^3]
mu_f = 5.0e-3            # Viscosidad efectiva a 60 C [Pa*s] (5 cP)

# Yacimiento (IPR Vogel)
P_r = 7.0e6             # Presion estatica del yacimiento [Pa]
q_max = 35.0 / 86400    # AOFP [m^3/s] (35 m^3/d)
h_din_0 = 800.0         # Nivel dinamico inicial [m]

# ---------------------------------------------------------------------------
# 6.7 Verificacion de coherencia (valores precalculados del documento)
# ---------------------------------------------------------------------------
g = 9.81                # Aceleracion gravitacional [m/s^2]

# Peso flotado de la sarta
W_rod_fl = W_rod * (1 - rho_f / rod_rho)  # = 25850 * 0.872 ≈ 22541 N

# Carga del fluido sobre el piston
F_fluid = rho_f * g * L_bomba * pump_A  # ≈ 16880 N

# PPRL y MPRL estimadas (sin efectos dinamicos)
PPRL_est = W_rod_fl + F_fluid  # ≈ 39420 N ≈ 8.86 klb
MPRL_est = W_rod_fl            # ≈ 22541 N ≈ 5.07 klb

# Caudal teorico a 6 SPM
q_th_6spm = pump_A * S_stroke * SPM_base / 60  # [m^3/s]

# ---------------------------------------------------------------------------
# 4.4 Parametros numericos
# ---------------------------------------------------------------------------
N_nodes = 100            # Nodos de discretizacion de la sarta
dx = rod_L / N_nodes     # Espaciado [m] = 11.0
dt = 2.15e-3             # Paso de tiempo [s] (elegido para alpha ~ 0.99)
alpha_CFL = rod_a * dt / dx  # Numero de Courant ≈ 0.987 (< 1, minima dispersion)
