# Modelo matemático y simulación dinámica de un Aparato Individual de Bombeo (AIB)

**Documento técnico del proyecto — versión 1.0**

---

## Índice

1. [Objetivo y alcance](#1-objetivo-y-alcance)
2. [Descripción física del sistema](#2-descripción-física-del-sistema)
3. [Marco teórico — formulación matemática](#3-marco-teórico--formulación-matemática)
4. [Métodos numéricos](#4-métodos-numéricos)
5. [Plan de implementación por fases](#5-plan-de-implementación-por-fases)
6. [Parámetros del caso de prueba](#6-parámetros-del-caso-de-prueba)
7. [Arquitectura del código](#7-arquitectura-del-código)
8. [Referencias](#8-referencias)

---

## 1. Objetivo y alcance

### 1.1 Objetivo general

Desarrollar un simulador en Python que reproduzca el comportamiento dinámico de un AIB (sucker-rod pumping unit) operando un pozo de petróleo, partiendo del **torque del motor eléctrico** como entrada primaria y resolviendo la propagación de fuerzas, posiciones y velocidades a lo largo de toda la cadena mecánica e hidráulica hasta el yacimiento.

### 1.2 Objetivos específicos

- Simular el comportamiento dinámico **forward** (fuerza-en, movimiento-out) del sistema, sin imponer la cinemática del balancín.
- Capturar las fluctuaciones de velocidad angular dentro del ciclo (típicamente ±5–15% del SPM nominal).
- Resolver la propagación de ondas elásticas en la sarta de varillas con amortiguamiento.
- Modelar el comportamiento de las válvulas TV y SV de la bomba de fondo, incluyendo carga de fluido y eventualmente patologías (golpe de fluido, llenado parcial, interferencia por gas).
- Acoplar la dinámica rápida del bombeo con la dinámica lenta del yacimiento (IPR Vogel + balance de nivel dinámico).
- Calcular las salidas de interés operativo: cartas dinamométricas de superficie y de fondo, torque neto en la caja reductora, potencia consumida por el motor, caudal producido y eficiencia global del sistema.

### 1.3 Alcance y simplificaciones del modelo base

Se incluyen explícitamente:

- Geometría tipo Clase I (convencional API), parametrizada por las dimensiones $A, C, I, H, P, R$ de la API Spec 11E.
- Sarta uniforme o multi-tapered (varios tramos con distinto diámetro).
- Motor de inducción tipo NEMA D modelado por curva torque–velocidad linealizada (extensible a modelo Klauss completo).
- Caja reductora rígida con relación $N$ y eficiencia $\eta$.
- Contrabalanceo en la manivela (rotativo).
- Bomba de pistón con válvulas TV (traveling valve) y SV (standing valve).
- IPR Vogel para flujo del yacimiento.
- Nivel dinámico variable en el espacio anular.

Se simplifican o excluyen del modelo base (extensibles en versiones posteriores):

- Inercia y dinámica detallada de la biela (se trata como cuerpo rígido con CG conocido).
- Pandeo lateral de la sarta (se considera solamente la dinámica axial 1D).
- Pérdidas por fricción varilla-tubing (se incorporan vía coeficiente de amortiguamiento global $\nu$).
- Compresibilidad de la columna de fluido en el tubing (modelo cuasi-estático).
- Gas libre dentro de la bomba (modelo monofásico).
- Pozos desviados (se asume pozo vertical).
- Transitorios eléctricos del motor (se asume tensión y frecuencia constantes).

---

## 2. Descripción física del sistema

### 2.1 Cadena de transmisión

El flujo de potencia y la propagación de variables ocurren en la siguiente secuencia:

```
Motor eléctrico → poleas + caja reductora → cigüeñal/manivela → biela (pitman)
    → balancín (walking beam) → cabezal (horse head) → varilla pulida (PR)
    → sarta de varillas (rod string) → émbolo (plunger) → bomba (TV + SV)
    → columna de fluido en el tubing → boca de pozo
```

Y en paralelo, la **realimentación física** (no de información sino de fuerzas):

```
Yacimiento (P_r, q_max) → entrada al espacio anular (annulus)
    → nivel dinámico h_din → presión de succión P_s en la bomba
    → estado de las válvulas → carga F_pump sobre la sarta
    → carga F_PR en la varilla pulida → torque sobre el cigüeñal
    → carga reflejada al motor
```

### 2.2 Variables observables del sistema

| Símbolo | Descripción | Unidad |
|---------|-------------|--------|
| $\theta$ | Ángulo del cigüeñal | rad |
| $\omega = \dot\theta$ | Velocidad angular del cigüeñal | rad/s |
| $y_{PR}$ | Posición vertical de la varilla pulida | m |
| $u(x,t)$ | Desplazamiento axial de la sarta a profundidad $x$ | m |
| $F_{PR}$ | Fuerza axial en la varilla pulida (carga del cabezal) | N |
| $F(x,t)$ | Fuerza axial en la sarta a profundidad $x$ | N |
| $T_m$ | Torque del motor | N·m |
| $T_g$ | Torque neto en la caja reductora (lado lento) | N·m |
| $h_{din}$ | Nivel dinámico de fluido en el anular | m |
| $P_{wf}$ | Presión de fondo fluyente | Pa |
| $q_{pump}$ | Caudal volumétrico instantáneo de la bomba | m³/s |

---

## 3. Marco teórico — formulación matemática

### 3.1 Análisis de grados de libertad

El sistema tiene grados de libertad de naturaleza heterogénea:

| Subsistema | DOFs aparentes | DOFs independientes | Tipo |
|------------|----------------|---------------------|------|
| Motor + reductor + cigüeñal + manivela + biela + balancín + PR | 7 | **1** ($\theta$) | mecánico rígido (mecanismo de 1 DOF) |
| Sarta de varillas (continuum) | $\infty$ | $N{-}1$ (discretizada en $N$ nodos) | continuum elástico (PDE) |
| Bomba | — | estados discretos (válvulas) | conmutado |
| Nivel dinámico | 1 | 1 | dinámica lenta |

La cadena rígida motor→PR colapsa a un único grado de libertad porque el cuadrilátero articulado cigüeñal-biela-balancín tiene 1 DOF, y las conexiones rígidas anteriores (motor, reductor, transmisión por poleas) no agregan DOFs cinemáticos independientes.

### 3.2 Cinemática del balancín (Clase I, API Spec 11E)

Sistema de coordenadas: saddle bearing (pivote del balancín) en el origen, centro del cigüeñal en $(I, -H)$. Para un ángulo de manivela $\theta$ (medido desde el eje X positivo):

$$
\begin{aligned}
\text{Crank pin:} \quad & (x_{cp}, y_{cp}) = (I + R\cos\theta,\; -H + R\sin\theta) \\
\text{Distancia SB→CP:} \quad & J(\theta) = \sqrt{x_{cp}^2 + y_{cp}^2} \\
\text{Ángulo en el SB:} \quad & \psi_b(\theta) = \arccos\!\left(\frac{C^2 + J^2 - P^2}{2\,C\,J}\right) \\
\text{Ángulo absoluto SB→CP:} \quad & \varphi(\theta) = \mathrm{atan2}(y_{cp}, x_{cp}) \\
\text{Equalizer bearing:} \quad & (x_E, y_E) = C\,(\cos(\varphi+\psi_b),\; \sin(\varphi+\psi_b)) \\
\text{Posición vertical de PR:} \quad & y_{PR}(\theta) = -\frac{A}{C}\,y_E
\end{aligned}
$$

A partir de $y_{PR}(\theta)$ se obtienen, por derivación numérica:

- **Velocidad de la varilla pulida:** $v_{PR}(\theta) = \dfrac{dy_{PR}}{d\theta}\,\dot\theta$
- **Aceleración de la varilla pulida:** $a_{PR}(\theta) = \dfrac{d^2 y_{PR}}{d\theta^2}\,\dot\theta^{\,2} + \dfrac{dy_{PR}}{d\theta}\,\ddot\theta$
- **Factor de torque:** $TF(\theta) = \dfrac{dy_{PR}}{d\theta}$ (m/rad)

El factor de torque es la pieza clave del acoplamiento mecánico: traduce fuerza vertical en la PR a torque sobre el cigüeñal.

### 3.3 Dinámica del mecanismo rígido (Lagrange)

Con $\theta$ como única coordenada generalizada, la **ecuación de Euler-Lagrange** se reduce a:

$$
\boxed{\;J_{ef}(\theta)\,\ddot\theta \;+\; \tfrac{1}{2}\,\frac{dJ_{ef}}{d\theta}\,\dot\theta^{\,2} \;=\; T_{neto}(\theta,\dot\theta,t)\;}
$$

#### Inercia equivalente proyectada

$$
J_{ef}(\theta) = N^2 J_m \;+\; J_{red} \;+\; J_c \;+\; J_{biela,ef}(\theta) \;+\; J_b\!\left(\!\frac{d\psi}{d\theta}\!\right)^{\!2} \;+\; m_{PR}\!\left(\!\frac{dy_{PR}}{d\theta}\!\right)^{\!2}
$$

donde:
- $N^2 J_m$: inercia del rotor del motor reflejada al lado lento del reductor
- $J_{red}$: inercia de los engranajes y eje del reductor (lado lento)
- $J_c$: inercia del cigüeñal y manivelas (incluido el contrapeso) alrededor del eje del cigüeñal
- $J_{biela,ef}(\theta)$: inercia efectiva de la biela (combinación de rotación y traslación de su CG)
- $J_b$: inercia del balancín alrededor del saddle bearing
- $\psi(\theta)$: ángulo del balancín, geométricamente determinado
- $m_{PR}$: masa del cabezal + varilla pulida + porción superior de la sarta (parte rígida en superficie)

#### Torque neto sobre el cigüeñal

$$
T_{neto} = N\,\eta\,T_m(\omega_m) \;-\; T_{cw}(\theta) \;-\; F_{PR}(t)\,TF(\theta) \;-\; T_{loss}(\dot\theta)
$$

con:

- $T_m(\omega_m)$: torque del motor (modelo eléctrico, ver §3.4)
- $\omega_m = N\,\dot\theta$: velocidad del eje del motor
- $T_{cw}(\theta) = M_{cw}\,g\,L_{cw}\,\sin(\theta - \tau_{cw})$: torque gravitacional del contrapeso
- $F_{PR}\,TF$: carga reflejada de la varilla pulida (el lazo de realimentación con la sarta)
- $T_{loss}$: pérdidas por fricción (modelo viscoso simple: $T_{loss} = c_f\,\dot\theta$ en primera aproximación)

### 3.4 Modelo del motor de inducción (NEMA D)

#### Modelo lineal (primera aproximación)

$$
T_m(\omega_m) = T_0 - k_m\,(\omega_m - \omega_0)
$$

con $\omega_0$ = velocidad sincrónica, $T_0$ = torque a $\omega = 0$ extrapolado, $k_m$ = pendiente de la región estable. Se calibra con dos puntos: $(\omega_{nom}, T_{nom})$ y $(\omega_0, 0)$.

#### Modelo Klauss (extensión)

$$
T_m(s) = \frac{2\,T_{\max}}{s/s_{\max} + s_{\max}/s}, \qquad s = \frac{\omega_0 - \omega_m}{\omega_0}
$$

donde $s$ es el deslizamiento (slip), $s_{\max}$ el slip al torque máximo. Captura mejor el comportamiento alejado del punto nominal y el arranque.

### 3.5 Sarta de varillas — ecuación de onda de Gibbs

#### Ecuación de gobierno (sarta uniforme)

$$
\frac{\partial^2 u}{\partial t^2} = a^2\,\frac{\partial^2 u}{\partial x^2} - c\,\frac{\partial u}{\partial t}
$$

con:
- $u(x,t)$: desplazamiento axial respecto al equilibrio estático no estirado
- $a = \sqrt{E/\rho}$: velocidad de propagación de onda longitudinal en el material (≈ 5050 m/s para acero)
- $c$: coeficiente de amortiguamiento (1/s), relacionado con el factor adimensional $\nu$ vía $c = \pi\,a\,\nu / (2L)$
- $\nu$: típicamente 0.05–0.20 (Gibbs–Neely)

#### Sarta multi-tapered

Para sarta con $K$ tramos de propiedades $(E_k, \rho_k, A_k, L_k)$, cada tramo cumple su propia ecuación de onda con $a_k = \sqrt{E_k/\rho_k}$. En cada interfaz entre tramos $k$ y $k{+}1$ deben cumplirse:

- Continuidad de desplazamiento: $u_k(x_k^-, t) = u_{k+1}(x_k^+, t)$
- Continuidad de fuerza axial: $E_k A_k\,\partial_x u_k = E_{k+1} A_{k+1}\,\partial_x u_{k+1}$

#### Condiciones de borde

- **Borde superior** (Dirichlet, viene de la cinemática del balancín):
  $$u(0, t) = y_{PR}(\theta(t)) - y_{PR,eq}$$
  con $y_{PR,eq}$ una posición de referencia (típicamente el punto medio del recorrido).

- **Borde inferior** (Neumann no lineal, viene de la bomba):
  $$E\,A_r\,\frac{\partial u}{\partial x}\bigg|_{x=L} = -F_{pump}(t)$$

#### Fuerza axial a cualquier profundidad

$$F(x,t) = E\,A_r(x)\,\frac{\partial u}{\partial x}$$

En particular, la **carga en la varilla pulida** que cierra el lazo con el mecanismo:

$$F_{PR}(t) = E\,A_r(0)\,\frac{\partial u}{\partial x}\bigg|_{x=0} + W_{r,susp}$$

donde $W_{r,susp}$ es el peso suspendido de la sarta en condiciones estáticas (referencia respecto a la cual se mide $u$).

### 3.6 Modelo de la bomba

#### Estados de las válvulas

La bomba alterna entre dos modos según la dirección del movimiento del émbolo:

| Modo | Velocidad émbolo | TV | SV | Carga sobre la sarta |
|------|------------------|----|----|----------------------|
| Carrera ascendente | $\dot u(L,t) > 0$ | cerrada | abierta | $F_{pump} = (P_d - P_s)\,A_p$ |
| Carrera descendente | $\dot u(L,t) < 0$ | abierta | cerrada | $F_{pump} = 0$ |

donde:
- $P_d$: presión de descarga (peso de la columna de fluido en el tubing)
- $P_s$: presión de succión (sumergencia de la bomba)
- $A_p$: área transversal del pistón

#### Presiones (modelo cuasi-estático)

$$P_d = P_{wh} + \rho_f\,g\,L_{bomba}$$
$$P_s = P_{ann} + \rho_f\,g\,(L_{bomba} - h_{din})$$

con $P_{wh}$ = presión en boca de pozo (manifold), $P_{ann}$ = presión en el anular (típicamente atmosférica + columna de gas), $L_{bomba}$ = profundidad de instalación de la bomba.

#### Caudal volumétrico instantáneo

$$q_{pump}(t) = A_p \cdot \max(0, \dot u(L,t)) \cdot \mathbb{1}_{\text{TV cerrada}}$$

(solo se desplaza fluido cuando la válvula viajera está cerrada y el émbolo sube)

#### Eficiencia volumétrica

El caudal real producido se reduce por:
- Llenado parcial del barril (si la bomba está sobredimensionada respecto al inflow del yacimiento)
- Escurrimiento entre pistón y barril (slippage)
- Compresibilidad de gas libre (en modelos avanzados)

Factor de eficiencia volumétrica típica: $\eta_v \in [0.70, 0.90]$.

### 3.7 Yacimiento — IPR Vogel y balance de nivel

#### Inflow Performance Relationship (Vogel, pozos saturados)

$$\frac{q_{in}}{q_{\max}} = 1 - 0.2\,\frac{P_{wf}}{P_r} - 0.8\!\left(\frac{P_{wf}}{P_r}\right)^{\!2}$$

con $q_{\max} = AOFP$ = caudal absoluto de pozo abierto, $P_r$ = presión estática del yacimiento.

#### Balance volumétrico del anular

$$A_{ann}\,\frac{dh_{din}}{dt} = q_{in}(P_{wf}) - q_{pump,prom}$$

donde $A_{ann}$ es el área anular (entre tubing y casing), $q_{pump,prom}$ es el caudal promedio de la bomba en el último ciclo. La constante de tiempo de esta ecuación es del orden de horas, mucho mayor que el ciclo de bombeo (segundos).

#### Presión de fondo fluyente

$$P_{wf} = P_{ann} + \rho_f\,g\,(L_{bomba} - h_{din})$$

---

## 4. Métodos numéricos

### 4.1 Ecuación de onda — diferencias finitas explícitas (leapfrog)

Discretización: malla espacial $x_i = i\,\Delta x$ con $i = 0, 1, \ldots, N$ y $\Delta x = L/N$. Esquema centrado en tiempo y espacio (orden 2):

$$
u_i^{n+1} = 2 u_i^n - u_i^{n-1} + \alpha^2\,(u_{i+1}^n - 2 u_i^n + u_{i-1}^n) - c\,\Delta t\,(u_i^n - u_i^{n-1})
$$

con $\alpha = a\,\Delta t / \Delta x$ (número de Courant). **Restricción CFL:** $\alpha \le 1$ (estricta para esquema explícito).

#### Condición de borde inferior con punto fantasma

La condición de Neumann $E A_r\,\partial_x u\big|_L = -F_{pump}$ se implementa con un nodo fantasma $u_{N+1}$:

$$
\frac{u_{N+1}^n - u_{N-1}^n}{2\,\Delta x} = -\frac{F_{pump}^n}{E A_r}
\;\Longrightarrow\;
u_{N+1}^n = u_{N-1}^n - \frac{2\,\Delta x\,F_{pump}^n}{E A_r}
$$

y luego se aplica la fórmula central en $i = N$ usando $u_{N+1}^n$.

#### Interfaces de sarta combinada

En el nodo $i^*$ donde cambia la sección, se imponen las condiciones de continuidad mediante una ecuación adicional que reemplaza la fórmula estándar en ese nodo. La derivación detallada se hará en Fase 2.

### 4.2 ODE del mecanismo — Runge-Kutta 4

Sistema de orden 1 equivalente (con $\omega = \dot\theta$):

$$
\begin{cases}
\dot\theta = \omega \\
\dot\omega = \dfrac{T_{neto}(\theta,\omega,t) - \tfrac{1}{2}(dJ_{ef}/d\theta)\,\omega^2}{J_{ef}(\theta)}
\end{cases}
$$

Integrado con RK4 de paso fijo $\Delta t$ (el mismo paso que la sarta para co-simulación).

### 4.3 Esquema de co-simulación con paso compartido

El acoplamiento ODE-PDE se resuelve con paso compartido y desfase de un tick (explícito en el acoplamiento):

```
Para cada paso n → n+1:
  1. Estado conocido: θ^n, ω^n, u_i^n, u_i^(n-1), nivel h_din^n, estado válvulas
  2. Calcular cinemática: y_PR(θ^n), TF(θ^n), J_ef(θ^n), dJ_ef/dθ
  3. Imponer BC superior: u_0^n = y_PR(θ^n)
  4. Calcular F_pump^n (lógica de válvulas según signo de v_émbolo)
  5. Avanzar PDE de la sarta un paso → u_i^(n+1)
  6. Calcular F_PR^n a partir de u_i^n
  7. Avanzar ODE del mecanismo con RK4 → θ^(n+1), ω^(n+1)
  8. Cada N_slow pasos: actualizar h_din con balance volumétrico
```

### 4.4 Estabilidad y elección de timestep

**Criterio dominante**: CFL de la sarta. Para $L = 1100$ m, $a = 5050$ m/s, $N = 100$ nodos:

$$\Delta x = 11\text{ m}, \qquad \Delta t_{CFL} = \Delta x / a = 2.18\text{ ms}$$

Se adopta $\Delta t = 1.5$ ms para tener margen ($\alpha \approx 0.69$).

**Frecuencia de bombeo**: 8 SPM = 0.133 Hz → período = 7.5 s → **5000 pasos por ciclo**. Para correr 3–5 ciclos transitorios + 5 ciclos de régimen permanente: ~50000 pasos. Tiempo de ejecución estimado en Python+NumPy vectorizado: 2–10 segundos.

### 4.5 Régimen permanente y convergencia

Para detectar régimen permanente se monitorea la **diferencia ciclo a ciclo** de variables clave:

$$\epsilon_k = \max_t |F_{PR}^{(k)}(t) - F_{PR}^{(k-1)}(t)| / \max_t |F_{PR}^{(k)}(t)|$$

Se considera convergido cuando $\epsilon_k < 10^{-3}$ en dos ciclos consecutivos.

---

## 5. Plan de implementación por fases

Cada fase produce un módulo verificable independientemente. La validación incluye criterios cuantitativos (rangos esperados) y cualitativos (forma de las curvas).

### Fase 1 — Cinemática del balancín ✓ (completada)

**Output**: $y_{PR}(\theta), v_{PR}(\theta), a_{PR}(\theta), TF(\theta)$ a $\omega$ constante.

**Validaciones**:

| Criterio | Esperado | Estado |
|----------|----------|--------|
| $v_{PR} = 0$ en BDC y TDC | sí | ✓ |
| Carrera $S \approx 2R\,(A/C)$ | aprox., depende de geometría | ✓ |
| Asimetría upstroke/downstroke | sí (~5–15° en Clase I) | ✓ |
| $|a_{PR}|_{\max}$ a 8 SPM | 0.1–0.5 g | ✓ (0.13 g) |

### Fase 2 — Dinámica de la sarta (PDE)

**Sub-fase 2a**: sarta uniforme, sin acoplamiento al balancín.

Implementar el integrador leapfrog con BCs simples (impuestas analíticamente).

**Validaciones**:

1. **Onda viajera sin amortiguamiento**: condición inicial $u(x,0) = \sin(\pi x/L)$, $\dot u(x,0) = 0$, BCs Dirichlet homogéneos. La solución debe oscilar manteniendo amplitud con período $T = 2L/a$. Tolerancia: error < 1% tras 10 períodos.

2. **Decaimiento exponencial**: mismo caso pero con $c > 0$. La amplitud debe decaer como $e^{-c t/2}$.

3. **Estiramiento estático**: BC superior fija $u(0,t) = 0$, BC inferior $F = -W_{rod}$ (peso de la sarta), sin amortiguamiento, condición inicial relajada. El estado estacionario debe ser $u(x) = -W_{rod}\,x\,(2L - x)/(2EA_r)$ aprox., y el estiramiento total $\Delta L = W_{rod}\,L/(EA_r)$. Para acero D-grade y $L = 1100$ m, esperar $\Delta L \approx 0.5$–1.5 m.

4. **Frecuencia natural fundamental**: la sarta libre-libre tiene $f_1 = a/(2L)$. Para nuestros datos: $f_1 \approx 2.3$ Hz. Verificar mediante FFT del decaimiento libre.

**Sub-fase 2b**: sarta combinada (multi-tapered), implementar interfaces.

**Validación**: para sarta de dos tramos iguales del mismo material, el resultado debe coincidir con la sarta uniforme equivalente.

### Fase 3 — Acoplamiento balancín → sarta

Imponer movimiento del balancín como BC superior de la sarta, con bomba ideal ($F_{pump}$ alternante simple según el signo de $\dot y_{PR}$).

**Validaciones**:

1. **Carga estática esperada**: en el upstroke, la carga en la PR debe ser aproximadamente
   $$F_{PR,up} \approx W_{r,fl} + (P_d - P_s)\,A_p$$
   donde $W_{r,fl}$ es el peso flotado de la sarta. Para nuestro caso: ~38–55 kN.

2. **Carga dinámica en el downstroke**: $F_{PR,down} \approx W_{r,fl}$ (la bomba ya no aporta). Para nuestro caso: ~22–25 kN.

3. **Forma de la carta dinamométrica de superficie**: $F_{PR}$ vs $y_{PR}$ debe ser un paralelogramo con esquinas redondeadas (ondas de tensión y compresión), girando en sentido antihorario (BDC abajo-derecha, TDC arriba-izquierda).

4. **PPRL/MPRL**: peak y minimum polished rod load. Verificar que $PPRL > MPRL$ y que la diferencia coincide aproximadamente con la carga del fluido más efectos dinámicos.

### Fase 4 — Modelo de la bomba con válvulas

Implementar la lógica de cambio de estado de TV/SV, incluyendo histéresis para evitar chattering. Calcular caudal instantáneo y carrera efectiva del émbolo.

**Validaciones**:

1. **Caudal teórico vs real**: el caudal nominal (sin pérdidas) es
   $$q_{th} = A_p \cdot S_p \cdot \text{SPM}/60$$
   donde $S_p$ es la carrera del émbolo. Para nuestro caso: ~22 m³/d. El caudal real debe ser 70–90% de este valor (eficiencia volumétrica).

2. **Carta dinamométrica de fondo**: $F_{plunger}$ vs posición del émbolo. Debe ser un rectángulo (idealizado) o paralelogramo con esquinas suavizadas. Las dos cargas características son: carga de fluido ($F_o$) durante upstroke, ~0 durante downstroke.

3. **Fluid pound check**: si la bomba está sobredimensionada (caudal superior al inflow), debe aparecer una caída brusca de carga al inicio del downstroke (golpe de fluido). Útil para validar el modelo en condiciones patológicas.

### Fase 5 — Dinámica del motor + mecanismo (forward dynamic)

Reemplazar la cinemática impuesta a $\omega$ constante por la integración de la EOM del mecanismo. El motor ahora aporta el torque de entrada según su curva, y la velocidad angular es resultado.

**Validaciones**:

1. **SPM promedio**: debe coincidir con el SPM de diseño (función del punto de operación motor + reductor + carga). Tolerancia: ±5%.

2. **Variación de $\omega$ dentro del ciclo**: típicamente ±5–15% del valor medio. Mayor variación indica contrabalanceo desbalanceado o motor sub-dimensionado.

3. **Curva de torque en la caja reductora**: $T_g(\theta)$ debe ser oscilante con dos picos por ciclo (upstroke y downstroke). El pico en upstroke típicamente es mayor en una unidad mal contrabalanceada.

4. **Balance energético**: trabajo neto del motor por ciclo = trabajo de elevación del fluido + pérdidas. Tolerancia: < 5%.

5. **Comparación cualitativa con Fase 3**: con un motor "fuerte" (k_m muy grande), debe recuperarse el comportamiento de cinemática impuesta.

### Fase 6 — Yacimiento + régimen permanente

Acoplamiento lento entre la dinámica de bombeo (escala de segundos) y la dinámica del nivel dinámico (escala de horas).

**Validaciones**:

1. **Convergencia del nivel dinámico**: para condiciones constantes, $h_{din}$ debe estabilizarse en un valor donde $q_{in}(P_{wf}) = q_{pump}$.

2. **Curva IPR-Outflow Performance**: el punto de operación del sistema corresponde a la intersección de la curva IPR (Vogel) con la curva del sistema de bombeo. Verificar gráficamente.

3. **Respuesta a perturbación**: si se aumenta el SPM repentinamente, $h_{din}$ debe descender (mayor extracción), $P_{wf}$ disminuye y $q_{in}$ aumenta hasta encontrar nuevo equilibrio.

### Fase 7 — Postprocesamiento

Generar todos los outputs operativos del simulador:

- **Cartas dinamométricas**: superficie y fondo (con ejes en las unidades estándar de la industria si se desea).
- **Torque en la caja reductora**: vs $\theta$, identificar pico, comparar con la capacidad nominal de la unidad.
- **Potencia consumida**: $P(t) = T_m\,\omega_m$, integrar para potencia promedio.
- **Eficiencia global**: $\eta_{sys} = P_{hidráulica útil} / P_{eléctrica entregada}$.
- **Reporte resumen**: PPRL, MPRL, torque pico, caudal, eficiencias, fluid pound flag.

---

## 6. Parámetros del caso de prueba

### 6.1 Pozo (Cuenca del Golfo San Jorge, caso típico)

| Parámetro | Símbolo | Valor | Unidad |
|-----------|---------|-------|--------|
| Profundidad de la bomba | $L_{bomba}$ | 1100 | m |
| Profundidad total del pozo | $L_{pozo}$ | 1200 | m |
| Diámetro interior del casing | $D_{csg}$ | 139.7 (5½") | mm |
| Diámetro exterior del tubing | $D_{tbg,o}$ | 73.0 (2⅞") | mm |
| Diámetro interior del tubing | $D_{tbg,i}$ | 62.0 | mm |
| Presión en boca de pozo | $P_{wh}$ | 0.5 (manométrica) | MPa |
| Presión en anular (gas) | $P_{ann}$ | 0.2 (manométrica) | MPa |

### 6.2 Unidad de bombeo (convencional API Clase I)

Geometría representativa de unidad mediana tipo **C-228D-200-86** operada a stroke de 144" (3.66 m):

| Parámetro | Símbolo | Valor | Unidad | Equivalente imperial |
|-----------|---------|-------|--------|----------------------|
| Brazo delantero del balancín | $A$ | 4.000 | m | 157.5 in |
| Brazo trasero del balancín | $C$ | 2.500 | m | 98.4 in |
| Offset horizontal del cigüeñal | $I$ | 2.000 | m | 78.7 in |
| Altura vertical SB → cigüeñal | $H$ | 4.500 | m | 177.2 in |
| Longitud de la biela (pitman) | $P$ | 3.500 | m | 137.8 in |
| Radio de la manivela | $R$ | 0.900 | m | 35.4 in |
| Carrera teórica resultante | $S$ | ~2.88 | m | ~113 in |
| SPM nominal de operación | — | 8 | ciclos/min | — |

**Inercias** (estimadas para una unidad mediana; refinar con catálogo del fabricante):

| Componente | Símbolo | Valor | Unidad |
|------------|---------|-------|--------|
| Cigüeñal + manivelas (sin contrapesos) | $J_c$ | 200 | kg·m² |
| Balancín (alrededor del SB) | $J_b$ | 800 | kg·m² |
| Biela (masa total) | $m_{pit}$ | 250 | kg |
| Cabezal + varilla pulida + adaptadores | $m_{PR}$ | 350 | kg |

**Contrabalanceo** (sobre la manivela):

| Parámetro | Símbolo | Valor | Unidad |
|-----------|---------|-------|--------|
| Masa del contrapeso (total entre dos manivelas) | $M_{cw}$ | 2500 | kg |
| Brazo efectivo | $L_{cw}$ | 0.85 | m |
| Fase respecto a la manivela | $\tau_{cw}$ | 180° | grados |

(Para que el contrabalanceo "asista" durante el upstroke, el contrapeso baja cuando la PR sube, lo que requiere $\tau_{cw} = 180°$ en la convención adoptada.)

### 6.3 Motor + transmisión (típico Argentina, 50 Hz)

| Parámetro | Símbolo | Valor | Unidad | Notas |
|-----------|---------|-------|--------|-------|
| Potencia nominal | $P_{nom}$ | 30 | kW (40 HP) | NEMA D |
| Frecuencia de red | $f$ | 50 | Hz | Argentina |
| Número de polos | — | 4 | — | — |
| Velocidad sincrónica | $\omega_0$ | 157.08 | rad/s (1500 rpm) | $4\pi f/P$ |
| Slip nominal | $s_{nom}$ | 0.12 | — | NEMA D |
| Velocidad nominal | $\omega_{nom}$ | 138.23 | rad/s (1320 rpm) | $\omega_0(1-s_{nom})$ |
| Torque nominal | $T_{nom}$ | 217 | N·m | $P_{nom}/\omega_{nom}$ |
| Pendiente lineal | $k_m$ | 11.5 | N·m·s/rad | $T_{nom}/(\omega_0-\omega_{nom})$ |
| Torque pull-out (Klauss) | $T_{\max}$ | 540 | N·m | 2.5× $T_{nom}$ |
| Slip al pull-out | $s_{\max}$ | 0.40 | — | típico NEMA D |
| Inercia del rotor | $J_m$ | 0.45 | kg·m² | depende del fabricante |

**Transmisión**:

| Parámetro | Símbolo | Valor |
|-----------|---------|-------|
| Relación poleas (motor → reductor rápido) | $N_p$ | 5.0 |
| Relación caja reductora | $N_g$ | 30.0 |
| Relación total | $N$ | 150 |
| Eficiencia de la transmisión | $\eta$ | 0.92 |

Verificación: $\omega_{crank} = \omega_{nom}/N = 138.23/150 = 0.922$ rad/s = 8.8 SPM. ✓ (cerca del SPM nominal de operación)

### 6.4 Sarta de varillas

**Configuración base** (para Fase 2): sarta uniforme 3/4" D-grade.

| Parámetro | Símbolo | Valor | Unidad |
|-----------|---------|-------|--------|
| Longitud total | $L$ | 1100 | m |
| Diámetro nominal de la varilla | $d_r$ | 19.05 (3/4") | mm |
| Área transversal | $A_r$ | 285 | mm² |
| Peso lineal en aire (con couplings) | $w_r$ | 23.5 | N/m (2.40 kg/m) |
| Módulo de Young | $E$ | 200 | GPa |
| Densidad del acero | $\rho_s$ | 7850 | kg/m³ |
| Velocidad de onda | $a$ | 5050 | m/s |
| Factor de amortiguamiento Gibbs | $\nu$ | 0.10 | adimensional |
| Coeficiente de amortiguamiento | $c$ | 0.722 | 1/s |

**Configuración avanzada** (para Fase 2b): sarta combinada 7/8" + 3/4".

| Tramo | Diámetro | Longitud |
|-------|----------|----------|
| Superior | 7/8" (22.23 mm) | 550 m |
| Inferior | 3/4" (19.05 mm) | 550 m |

### 6.5 Bomba

| Parámetro | Símbolo | Valor | Unidad |
|-----------|---------|-------|--------|
| Diámetro del pistón | $d_p$ | 44.45 (1¾") | mm |
| Área del pistón | $A_p$ | 1552 | mm² |
| Tipo | — | Tubing pump | — |
| Eficiencia volumétrica esperada | $\eta_v$ | 0.80 | — |

### 6.6 Fluido y yacimiento

**Fluido producido** (mezcla agua + petróleo, alto WC típico Golfo San Jorge):

| Parámetro | Símbolo | Valor | Unidad |
|-----------|---------|-------|--------|
| Corte de agua | WC | 90 | % |
| Densidad del petróleo (18 °API) | $\rho_o$ | 947 | kg/m³ |
| Densidad del agua de formación | $\rho_w$ | 1015 | kg/m³ |
| Densidad efectiva mezcla | $\rho_f$ | 1008 | kg/m³ |
| Viscosidad efectiva (60°C) | $\mu_f$ | 5 | cP |

**Yacimiento (IPR Vogel)**:

| Parámetro | Símbolo | Valor | Unidad |
|-----------|---------|-------|--------|
| Presión estática | $P_r$ | 7.0 | MPa |
| AOFP (caudal a $P_{wf}$=0) | $q_{\max}$ | 35 | m³/d |
| Nivel dinámico inicial | $h_{din,0}$ | 800 | m |

### 6.7 Verificación de coherencia preliminar

**Cargas estáticas estimadas** (sin dinámica):

- Peso de la sarta en aire: $W_r = w_r \cdot L = 23.5 \cdot 1100 = 25.85$ kN
- Peso flotado: $W_{r,fl} = W_r \cdot (1 - \rho_f/\rho_s) = 25.85 \cdot 0.872 = 22.54$ kN
- Carga del fluido sobre el pistón (tubing lleno):
  $$F_o = \rho_f\,g\,L_{bomba}\,A_p = 1008 \cdot 9.81 \cdot 1100 \cdot 1.552\times 10^{-3} = 16.88\text{ kN}$$
- **PPRL estimada** (sin efectos dinámicos): $W_{r,fl} + F_o = 39.4$ kN ≈ 8.86 klb
- **MPRL estimada** (sin efectos dinámicos): $W_{r,fl} = 22.5$ kN ≈ 5.06 klb

Capacidad estructural de la unidad C-228D-200-86: 86 klb. **Uso al 10–12%**. Esto indica que la unidad está sobredimensionada, lo cual es habitual en pozos de Golfo San Jorge donde se priorizan unidades robustas y de bajo mantenimiento.

**Caudal teórico estimado**:

$$q_{th} = A_p \cdot S \cdot \text{SPM}/60 = 1.552\times 10^{-3} \cdot 2.88 \cdot 8/60 = 5.96\times 10^{-4}\text{ m³/s} = 51.5\text{ m³/d}$$

Con $\eta_v = 0.80$: $q_{real} \approx 41$ m³/d. Como $q_{\max,res} = 35$ m³/d, la bomba **podría estar sobredimensionada**, lo cual produciría llenado parcial y golpe de fluido. Esto es un escenario interesante para validar el modelo en condiciones patológicas.

Para evitar este efecto en el caso base, se puede:

- Reducir el SPM a 6 → $q_{th} \approx 38.6$ m³/d, mejor ajuste.
- O reducir la carrera operativa.
- O aumentar la productividad del yacimiento ($q_{\max}$).

**Decisión para el caso base**: operar a **6 SPM** para evitar fluid pound en la primera corrida. Reservar el caso de 8 SPM para análisis de patología.

---

## 7. Arquitectura del código

### 7.1 Estructura de módulos

```
aib_simulator/
├── core/
│   ├── geometry.py          # Cinemática del balancín (Fase 1)
│   ├── rod_string.py        # PDE de la sarta (Fase 2)
│   ├── pump.py              # Lógica de válvulas y carga (Fase 4)
│   ├── motor.py             # Curva torque-velocidad (Fase 5)
│   ├── reservoir.py         # IPR y nivel dinámico (Fase 6)
│   └── mechanism.py         # EOM del mecanismo + co-simulación (Fase 5)
├── solver/
│   ├── wave_solver.py       # Esquema leapfrog FD
│   ├── ode_solver.py        # RK4 para la EOM
│   └── coupling.py          # Orquestador de la co-simulación
├── postprocess/
│   ├── dynamometer.py       # Cartas de superficie y fondo
│   ├── torque_analysis.py   # Curva de torque, contrabalanceo
│   └── reports.py           # Resumen operativo
├── data/
│   ├── case_GSJ_typical.yaml  # Parámetros del caso de prueba
│   └── motor_catalog.yaml     # Curvas de motores comerciales
├── tests/
│   └── test_*.py             # Tests unitarios por módulo
└── main.py                   # Entry point
```

### 7.2 Convenciones

- **Unidades internas**: SI estricto (m, kg, s, N, Pa, rad). Conversión a unidades imperiales solo en presentación.
- **Variables de estado**: arrays NumPy de doble precisión.
- **Configuración**: archivos YAML con todos los parámetros del caso, separados por subsistema.
- **Logging**: en cada paso se almacenan las variables clave en un buffer circular para postprocesamiento; al final, se exportan a HDF5 o pickle.
- **Tests**: cada módulo tiene tests unitarios con casos analíticos cerrados (los descritos en las validaciones).

### 7.3 Dependencias

Mínimas (Fases 1–4):

- `numpy >= 1.24`
- `scipy >= 1.10` (para integradores de referencia y FFT)
- `matplotlib >= 3.7` (visualización)
- `pyyaml` (configuración)

Aceleración (opcional):

- `numba >= 0.58` (JIT para los lazos críticos del solver de la PDE)

Análisis interactivo:

- `jupyter`
- `pandas` (manejo de datos de salida)

---

## 8. Referencias

### Fundamentos teóricos

1. **Gibbs, S. G. (1963).** *Predicting the behavior of sucker-rod pumping systems.* JPT 15(7): 769–778. [Paper fundacional de la ecuación de onda aplicada a sucker-rod pumps.]
2. **Gibbs, S. G., Neely, A. B. (1966).** *Computer diagnosis of downhole conditions in sucker-rod pumping wells.* JPT 18(1): 91–98. [Inversión de la carta de superficie a la de fondo.]
3. **Everitt, T. A., Jennings, J. W. (1992).** *An improved finite-difference calculation of downhole dynamometer cards for sucker-rod pumps.* SPE Production Engineering 7(1): 121–127.
4. **Doty, D. R., Schmidt, Z. (1983).** *An improved model for sucker-rod pumping.* SPEJ 23(1): 33–41. [Acoplamiento sarta-bomba con dinámica de fluido.]

### Manuales y referencias de ingeniería

5. **Takács, G. (2015).** *Sucker-Rod Pumping Manual* (2nd ed.). PennWell. [Referencia moderna integral, incluye geometría API y ejemplos de cálculo.]
6. **API Spec 11E.** *Specification for Pumping Units.* American Petroleum Institute. [Geometrías estandarizadas Clase I, II, III.]
7. **API RP 11L.** *Recommended Practice for Design Calculations for Sucker-Rod Pumping Systems.* [Método clásico de diseño preliminar.]
8. **Lufkin Oilfield Manufacturing.** *Pumping Unit Performance Catalog.* [Datos de inercia y geometría de unidades comerciales.]

### IPR y comportamiento de yacimiento

9. **Vogel, J. V. (1968).** *Inflow performance relationships for solution-gas drive wells.* JPT 20(1): 83–92.
10. **Brown, K. E. (1980).** *The Technology of Artificial Lift Methods, Vol. 4: Production Optimization of Oil and Gas Wells.* PennWell.

### Métodos numéricos

11. **LeVeque, R. J. (2007).** *Finite Difference Methods for Ordinary and Partial Differential Equations.* SIAM. [Esquemas para la ecuación de onda.]
12. **Press, W. H. et al. (2007).** *Numerical Recipes* (3rd ed.). Cambridge University Press. [RK4, FFT, métodos de raíces.]

### Recursos para el contexto argentino

13. **Cuenca del Golfo San Jorge** — datos públicos de producción y características de pozos disponibles en publicaciones de la Secretaría de Energía y de IAPG (Instituto Argentino del Petróleo y del Gas).
14. **MEC S.A.** (Catriel, Río Negro) — fabricante nacional de unidades AIB tipo C convencional.
15. **AMELCO Formación Profesional** (Neuquén) — material didáctico sobre componentes de AIB.

---

**Fin del documento — versión 1.0**

*Próximas actualizaciones: detalles de implementación de la sarta combinada (Fase 2b), modelo Klauss completo del motor con transitorios eléctricos (Fase 5 extendida), y modelo de gas libre en la bomba (Fase 4 extendida).*
