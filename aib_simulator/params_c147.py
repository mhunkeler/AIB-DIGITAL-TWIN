"""
Parametros del pozo C147 — Mark II M-912D-365-168 (CAPSA, Diadema)

Datos firmes: geometria API 11E, sarta de varillas, bomba, motor HP/RPM.
Datos estimados: inercias, transmision, contrapesos, amortiguamiento, yacimiento.
Los estimados se ajustan con el notebook 08_ajuste_c147.ipynb.

Fuente de datos firmes: M-912D-365-168.json, wells_atlas.json
"""

import numpy as np

# ---------------------------------------------------------------------------
# Geometria Mark II (M-912D-365-168) — DATOS FIRMES (del JSON)
# ---------------------------------------------------------------------------
# En Mark II: C = pitman, P = rocker (SB→equalizer), ratio = C/P
A_beam = 334.0 * 0.0254        # 8.484 m (SB → horsehead)
P_rocker = 193.5 * 0.0254      # 4.915 m (SB → equalizer bearing)
C_pitman = 270.0 * 0.0254      # 6.858 m (longitud de la biela)
I_geom = 202.56 * 0.0254       # 5.145 m (offset horizontal SB → ciguenal)
H_abs = 295.13 * 0.0254        # 7.496 m (altura SB sobre base)
G_height = 112.13 * 0.0254     # 2.848 m (altura ciguenal sobre base)
H_offset = H_abs - G_height    # 4.648 m (offset vertical SB → ciguenal)
R_crank = 63.56 * 0.0254       # 1.614 m (stroke 1 = 168")

S_stroke = 165.6 * 0.0254      # 4.206 m (carrera medida en campo)

# Ratings del catalogo
max_torque = 912000             # in-lbs (capacidad del reductor)
PPRL_rated = 36500              # lbs (carga maxima de diseno)
structural_unbalance = -5385    # lbs
catalog_TF = 75.207             # in/rad (torque factor API)
catalog_phase_angle = 19        # grados

# ---------------------------------------------------------------------------
# Motor — DATOS FIRMES (HP y RPM del JSON)
# ---------------------------------------------------------------------------
motor_HP = 100
motor_RPM = 1200
motor_pulley_mm = 220

# Motor electrico (50 Hz Argentina, 6 polos para 1200 RPM sincronica)
# 1200 RPM sincronica → 6 polos a 60 Hz O 8 polos a 60 Hz...
# En Argentina (50 Hz): sync = 120*f/polos
# Para ~1200 RPM: 120*50/5 = 1200 (no existe 5 polos)
# 6 polos: 120*50/6 = 1000 RPM sync → motor de 1000 RPM
# 4 polos: 120*50/4 = 1500 RPM sync
# El motor de 1200 RPM nominal probablemente es 4 polos con slip alto
# o es un motor de 60 Hz (importado). Usamos 4 polos 50 Hz.
f_red = 50.0                    # Hz (Argentina)
n_polos = 4
omega_0 = 2 * np.pi * f_red / (n_polos / 2)  # 157.08 rad/s (1500 RPM sync)

P_nom = motor_HP * 745.7       # 74,570 W
s_nom = 0.08                   # slip nominal (100 HP NEMA B/C)
omega_nom = omega_0 * (1 - s_nom)  # 144.51 rad/s (1380 RPM)
T_nom = P_nom / omega_nom      # 516 N*m
k_m = T_nom / (omega_0 * s_nom)  # 41.1 N*m*s/rad

# Klauss
T_max = 2.5 * T_nom            # 1290 N*m
s_max = 0.25                   # pull-out slip

J_m = 1.5                      # kg*m² (rotor 100 HP, estimado)

# ---------------------------------------------------------------------------
# Transmision — ESTIMADO (a ajustar)
# ---------------------------------------------------------------------------
# Motor 1500 sync → ciguenal ~5-8 SPM
# Si SPM = 6: N_total = 1380/6/2pi*60 ≈ 1380/(0.628) = 2197... no
# N_total = RPM_motor / RPM_ciguenal = 1380 / (SPM/60*2pi...
# RPM_ciguenal = SPM (son lo mismo en unidades de rev/min)
# N_total = 1380 / 6 = 230 (para 6 SPM)
# N_total = 1380 / 8 = 172.5 (para 8 SPM)
N_pulley = 3.5                  # relacion de poleas (220mm motor / ~770mm driven)
N_gearbox = 40.0                # reductor doble reduccion
N_total = N_pulley * N_gearbox  # = 140 → SPM ≈ 1380/140 = 9.9
eta_trans = 0.90                # eficiencia de la transmision

SPM_nominal = 6.0               # SPM nominal (del informe de campo)
SPM_base = 6.0

# ---------------------------------------------------------------------------
# Sarta de varillas — DATOS FIRMES (del JSON)
# ---------------------------------------------------------------------------
rod_E = 200.0e9                 # Modulo de Young acero [Pa]
rod_rho = 7850.0                # Densidad del acero [kg/m^3]
rod_a = np.sqrt(rod_E / rod_rho)  # 5050 m/s

# Sarta combinada: 3 tramos D-grade
rod_tapered = [
    {
        'name': '1" D-grade (superior)',
        'd': 1.000 * 0.0254,       # 25.4 mm
        'A': np.pi / 4 * (1.000 * 0.0254)**2,  # 506.7 mm²
        'L': 91 * 30 * 0.3048,     # 832.1 m (2730 ft)
        'E': rod_E,
        'rho': rod_rho,
        'w': 2.904 * 4.4482 / 0.3048,  # 42.4 N/m
        'count': 91,
    },
    {
        'name': '7/8" D-grade (medio)',
        'd': 0.875 * 0.0254,       # 22.2 mm
        'A': np.pi / 4 * (0.875 * 0.0254)**2,  # 387.8 mm²
        'L': 111 * 30 * 0.3048,    # 1014.7 m (3330 ft)
        'E': rod_E,
        'rho': rod_rho,
        'w': 2.224 * 4.4482 / 0.3048,  # 32.4 N/m
        'count': 111,
    },
    {
        'name': '3/4" D-grade (inferior)',
        'd': 0.750 * 0.0254,       # 19.05 mm
        'A': np.pi / 4 * (0.750 * 0.0254)**2,  # 285.0 mm²
        'L': 96 * 30 * 0.3048,     # 877.8 m (2880 ft)
        'E': rod_E,
        'rho': rod_rho,
        'w': 1.634 * 4.4482 / 0.3048,  # 23.8 N/m
        'count': 96,
    },
]

rod_L = 8940 * 0.3048           # 2724.9 m total
rod_nu = 0.10                   # factor de amortiguamiento Gibbs
rod_c = np.pi * rod_a * rod_nu / (2 * rod_L)  # 0.291 1/s

# Para el solver uniforme: area promedio ponderada por longitud
_total_L = sum(s['L'] for s in rod_tapered)
rod_A = sum(s['A'] * s['L'] for s in rod_tapered) / _total_L  # 391 mm²
rod_w = sum(s['w'] * s['L'] for s in rod_tapered) / _total_L
rod_d = np.sqrt(4 * rod_A / np.pi)

# Pesos de la sarta (del JSON)
W_rod = 20040 * 4.4482          # 89,142 N en aire
W_rod_fl = 17514 * 4.4482       # 77,907 N en fluido

# ---------------------------------------------------------------------------
# Bomba — DATOS FIRMES (del JSON)
# ---------------------------------------------------------------------------
pump_d = 1.75 * 0.0254          # 0.04445 m
pump_A = np.pi / 4 * pump_d**2  # 1552 mm²
pump_eta_v = 0.80               # eficiencia volumetrica

# ---------------------------------------------------------------------------
# Pozo — DATOS FIRMES + ESTIMADOS
# ---------------------------------------------------------------------------
L_bomba = 2740.0                # m (profundidad de la bomba, del JSON)
L_pozo = 2800.0                 # m (profundidad total, estimado)
D_csg = 0.13970                 # m (casing 5-1/2", estimado)
D_tbg_o = 0.07300               # m (tubing 2-7/8" OD, estimado)
D_tbg_i = 0.06200               # m (tubing ID, estimado)
P_wh = 0.5e6                   # Pa (presion boca de pozo, estimado)
P_ann = 0.2e6                  # Pa (presion anular, estimado)

A_ann = np.pi / 4 * (D_csg**2 - D_tbg_o**2)

# ---------------------------------------------------------------------------
# Fluido — ESTIMADO (tipico Diadema, alto corte de agua)
# ---------------------------------------------------------------------------
WC = 0.90
rho_f = 1008.0                  # kg/m³
mu_f = 5.0e-3                   # Pa*s (5 cP)

# ---------------------------------------------------------------------------
# Yacimiento — ESTIMADO (tipico Golfo San Jorge)
# ---------------------------------------------------------------------------
P_r = 8.0e6                    # Pa (presion estatica, estimado)
q_max = 30.0 / 86400           # m³/s (30 m3/d AOFP, estimado)
h_din_0 = 900.0                # m (nivel dinamico inicial, estimado)

# ---------------------------------------------------------------------------
# Inercias — ESTIMADO (escalar desde peso de la unidad 40,000 lbs)
# ---------------------------------------------------------------------------
J_c = 800.0                    # kg*m² (ciguenal + manivelas, sin CW)
J_b = 3000.0                   # kg*m² (balancin alrededor del SB)
m_pit = 500.0                  # kg (masa de la biela)
m_PR = 700.0                   # kg (cabezal + varilla pulida + adaptadores)

# ---------------------------------------------------------------------------
# Contrapeso — ESTIMADO (4 OORO blocks)
# ---------------------------------------------------------------------------
M_cw = 5000.0                  # kg (4 bloques ~1250 kg c/u)
L_cw = 1.20                    # m (brazo efectivo)
tau_cw = np.pi / 2             # 90 deg (a optimizar)

# ---------------------------------------------------------------------------
# Friccion — ESTIMADO
# ---------------------------------------------------------------------------
c_friction = 5.0                # N*m*s/rad

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------
g = 9.81

# ---------------------------------------------------------------------------
# Parametros numericos
# ---------------------------------------------------------------------------
N_nodes = 200                   # nodos (sarta mas larga → mas nodos)
dx = rod_L / N_nodes            # 13.6 m
dt = 2.65e-3                    # paso de tiempo (alpha ~0.98 para dx=13.6)
alpha_CFL = rod_a * dt / dx    # ~0.98

# Verificacion CFL
assert alpha_CFL < 1.0, f"CFL violada: alpha = {alpha_CFL:.3f}"

# ---------------------------------------------------------------------------
# Valores de referencia del campo (para validacion)
# ---------------------------------------------------------------------------
field_F_min_lbs = 7580.53       # fuerza minima medida [lbs]
field_F_max_lbs = 27014.31      # fuerza maxima medida [lbs]
field_F_min = field_F_min_lbs * 4.4482  # 33,717 N
field_F_max = field_F_max_lbs * 4.4482  # 120,170 N
field_stroke_in = 165.6         # carrera medida [in]
