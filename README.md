# AIB Digital Twin — Lufkin Mark II M-912D-365-168

## Modelo de Referencia

El **Lufkin Mark II M-912D-365-168** es una de las unidades de bombeo mecanico mas
utilizadas en los yacimientos de la Cuenca del Golfo San Jorge (Chubut/Santa Cruz, Argentina),
operados por empresas como CAPSA, YPF, PAE, Tecpetrol y Sinopec.

Es el mismo modelo instalado en el pozo de prueba **C147 (CAPSA, Diadema)**.

## Designacion API 11E

```
M - 912D - 365 - 168
|    |      |     |
|    |      |     └─ Carrera maxima: 168 pulgadas
|    |      └─ PPRL (Peak Polished Rod Load): 36,500 lbs
|    └─ Torque maximo reductor: 912,000 in-lbs
└─ Mark II (geometria Unitorque)
```

## Especificaciones Principales

| Parametro | Valor |
|-----------|-------|
| **Tipo** | Mark II (Unitorque) |
| **Fabricante** | Lufkin Industries (ahora GE Oil & Gas / Baker Hughes) |
| **Torque maximo** | 912,000 in-lbs |
| **PPRL** | 36,500 lbs |
| **CCE** | 912,000 |
| **Carrera maxima** | 168 in (4.27 m) |
| **Carreras disponibles** | 168 / 149 / 130 in |
| **Structural Unbalance (SU)** | -5,385 lbs |
| **Torque Factor (TF)** | 75.207 |
| **Phase Angle** | 19 degrees |

## Geometria del Mecanismo (API 11E)

```
                    Walking Beam
         ┌─────────────────────────────┐
         │            A                │
    ─────┼──────────┬──────────────────┤
         │          │  Saddle Bearing   │ Horsehead
         │          │                  │
         │     P    │                  │
         │          │                  │
         │     Pitman                  
         │          │                  
    ─────┼──── C ───┤ Equalizer        
         │          │ Bearing          
         │          │                  
         │     Crank│                  
         │     R    │                  
    ═════╪══════════╧═══════════ Crankshaft
         │
         │ Base
```

### Dimensiones del Mecanismo

| Parametro | Valor | Descripcion |
|-----------|-------|-------------|
| **A** | 334.0 in (8.48 m) | Distancia del saddle bearing al horsehead (walking beam front) |
| **C** | 270.0 in (6.86 m) | Longitud del pitman (connecting rod) |
| **I** | 202.56 in (5.14 m) | Distancia horizontal del crankshaft al saddle bearing |
| **P** | 193.5 in (4.91 m) | Distancia del equalizer bearing al saddle bearing |
| **H** | 295.13 in (7.50 m) | Altura del saddle bearing sobre la base |
| **G** | 112.13 in (2.85 m) | Altura del crankshaft sobre la base |

### Radios de Manivela y Carreras

| Carrera | Radio (R) | Stroke |
|---------|-----------|--------|
| Stroke 1 (max) | R1 = 63.56 in | 168 in |
| Stroke 2 (med) | R2 = 56.56 in | 149 in |
| Stroke 3 (min) | R3 = 49.56 in | 130 in |

## Dimensiones Generales Estimadas

| Parametro | Valor |
|-----------|-------|
| **Altura total** | ~25 ft (7.6 m) — desde base hasta horsehead en posicion superior |
| **Longitud total** | ~30 ft (9.1 m) — desde contrapeso trasero hasta pozo |
| **Ancho** | ~6 ft (1.8 m) |
| **Peso total** | ~35,000 - 45,000 lbs (16 - 20 ton) — sin contrapesos |
| **Peso con contrapesos** | ~50,000 - 65,000 lbs (23 - 30 ton) |

## Motor

| Parametro | Valor (tipico C147) |
|-----------|-------------------|
| Tipo | Electrico trifasico |
| Potencia | 30 - 100 HP |
| RPM | 1200 |
| Polea motor | 220 mm |
| Reductor | Doble reduccion |

## Contrapesos

El Mark II usa contrapesos rotatorios montados en el brazo de la manivela:

| Tipo | Cantidad | Peso unitario |
|------|----------|---------------|
| OORO (tipico C147) | 4 | ~1,500 - 2,000 lbs c/u |

Los contrapesos se ajustan en distancia al eje para balancear la carga del pozo.

## Sarta de Varillas (C147)

| Tramo | Diametro | Material | Cantidad | Longitud | Peso/pie |
|-------|----------|----------|----------|----------|----------|
| 1 (superior) | 1" | D | 91 varillas | 30 ft c/u | 2.904 lb/ft |
| 2 (medio) | 7/8" | D | 111 varillas | 30 ft c/u | 2.224 lb/ft |
| 3 (inferior) | 3/4" | D | 96 varillas | 30 ft c/u | 1.634 lb/ft |

**Longitud total**: 2,730 + 3,330 + 2,880 = **8,940 ft (2,725 m)**
**Peso en aire**: 7,928 + 7,406 + 4,706 = **20,040 lbs**
**Peso en fluido** (buoyancy ~0.87): **~17,514 lbs**

## Bomba de Profundidad

| Parametro | Valor |
|-----------|-------|
| Diametro piston | 1.75 in |
| Profundidad | 2,740 m |
| Tipo | Tubing pump |

## Velocidad de Operacion

| Parametro | Rango tipico |
|-----------|-------------|
| GPM (golpes por minuto) | 2 - 6 |
| Periodo del ciclo | 10 - 30 segundos |
| Velocidad del motor | 1200 RPM |
| Relacion de reduccion | ~30:1 a 50:1 |

## Aplicacion en Argentina

El Lufkin Mark II es ampliamente utilizado en:
- **Cuenca del Golfo San Jorge** (Chubut/Santa Cruz): CAPSA, YPF, PAE, Sinopec
- **Cuenca Neuquina** (Neuquen/Mendoza): YPF, Pluspetrol, Tecpetrol
- **Cuenca Cuyana** (Mendoza): YPF

Condiciones tipicas:
- Profundidades: 800 - 3,000 m
- Temperaturas ambiente: -15 C a +45 C
- Produccion: 5 - 50 m3/dia
- Corte de agua: 30% - 95%

## Fuentes

- [Lufkin Mark II Unitorque Pumping Unit](https://www.thehistorycenteronline.com/uploads/resources/Lufkin_Mark_II_Unitorque_Pumping_Unit.pdf)
- [Lufkin Pumping Units Catalog](https://www.scribd.com/doc/275118013/Catalog-Pumping-Units-Lufkin-pdf)
- [Mark II Installation Manual](https://www.yumpu.com/en/document/view/11297459/mark-ii-pumping-units-installation-manual-lufkin-industries)
- [API SPEC 11E Specifications for Pumping Units](https://standards.globalspec.com/std/14568280/spec-11e)
- Datos internos del proyecto Well Optimus (structural_data_raw.py, config C147)
