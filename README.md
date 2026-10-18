# Procesamiento y Optimización de Imágenes Satelitales

**Optimización y Heurística | Proyecto Universitario | MATLAB**

![MATLAB](https://img.shields.io/badge/MATLAB-Procesamiento%20de%20Imágenes-orange)
![Sentinel-2](https://img.shields.io/badge/Datos-Sentinel--2-blue)
![Estado](https://img.shields.io/badge/Estado-En%20desarrollo-yellow)
![Universidad](https://img.shields.io/badge/Universidad-ULPGC-green)

## Descripción del proyecto

Este repositorio alberga un proyecto educativo que se enfoca en el análisis y mejora de imágenes satelitales utilizando métodos de optimización y heurística, empleando imágenes multiespectrales del satélite Sentinel-2 junto con la plataforma de programación MATLAB.

El área de estudio elegida se ubica en Telde, Gran Canaria (Canarias, España), e incluye zonas urbanas, áreas verdes, costas y cuerpos de agua.

El proyecto integra técnicas de teledetección, análisis digital de imágenes y estrategias de optimización para investigar cómo diversas configuraciones de parámetros afectan la identificación y clasificación de diferentes superficies.

La meta fundamental es crear una metodología que sea reproducible para el procesamiento de imágenes satelitales, formular un problema de optimización y evaluar varias estrategias de solución a través de experimentos computacionales.

---

## Objetivos

Los objetivos principales de esta iniciativa son:

1. Obtener y preparar imágenes multiespectrales del satélite Sentinel-2.
2. Visualizar las imágenes satelitales mediante composiciones RGB.
3. Identificar píxeles válidos y descartar datos no confiables.
4. Calcular índices espectrales para el análisis de vegetación, masas de agua y superficies edificadas.
5. Definir un problema de optimización vinculado a la clasificación de superficies.
6. Implementar y comparar diferentes algoritmos de optimización heurística.
7. Evaluar la calidad de las soluciones, la eficiencia en el uso de recursos computacionales y la solidez de los métodos.
8. Analizar los hallazgos y extraer conclusiones fundamentadas.

---

## Zona de estudio y datos utilizados

**Localización:** Telde, Gran Canaria, Islas Canarias.

**Satélite:** Sentinel-2

**Sensor:** Multispectral Instrument (MSI)

**Fecha de adquisición:** 28 de septiembre de 2026

**Fuente de datos:** Copernicus Data Space Ecosystem

**Entorno de desarrollo:** MATLAB

La zona seleccionada presenta diferentes tipos de cobertura terrestre, lo que permite estudiar el comportamiento de los índices espectrales en distintos entornos.

### Bandas espectrales utilizadas

| Banda | Descripción | Resolución original |
|---|---|---|
| B02 | Azul | 10 m |
| B03 | Verde | 10 m |
| B04 | Rojo | 10 m |
| B08 | Infrarrojo cercano (NIR) | 10 m |
| B11 | Infrarrojo de onda corta (SWIR) | 20 m |

Las bandas exportadas se trabajan sobre una cuadrícula común de 10 metros. Esto no implica que la banda B11 adquiera mayor detalle espacial que su resolución original.

---

## Metodología

El desarrollo del proyecto se organiza en varias etapas.

### 1. Adquisición y preprocesamiento

Las imágenes satelitales son obtenidas del navegador de datos de Copernicus.

Luego, las bandas seleccionadas se importan a MATLAB para llevar a cabo:

- Lectura de archivos TIFF.
- Verificación de dimensiones.
- Detección de píxeles sin datos válidos.
- Creación de una máscara de validez.
- Preparación de las bandas para su análisis.

### 2. Visualización RGB

Se crea una composición en color natural utilizando las siguientes bandas:

- Canal rojo: B04.
- Canal verde: B03.
- Canal azul: B02.

Esta representación facilita la interpretación visual del área de estudio y la identificación de elementos geográficos

### 3. Cálculo de índices espectrales

Se calculan tres índices espectrales.

**NDVI — Índice de Vegetación de Diferencia Normalizada**

NDVI = (B08 - B04) / (B08 + B04)

Permite analizar la presencia y el estado relativo de la vegetación.

**NDWI — Índice de Agua de Diferencia Normalizada**

NDWI = (B03 - B08) / (B03 + B08)

Se utiliza para resaltar superficies de agua.

**NDBI — Índice de Áreas Construidas de Diferencia Normalizada**

NDBI = (B11 - B08) / (B11 + B08)

Se emplea para analizar superficies construidas, aunque otros tipos de cobertura pueden presentar valores similares.

### 4. Formulación del problema de optimización

La etapa siguiente implica la definición de un problema de optimización vinculado a los umbrales de clasificación que se derivan de los índices espectrales.

Se explorarán diversas estrategias para identificar cuál configuración produce mejores resultados.

Los algoritmos candidatos incluyen:

- **Optimización por Enjambre de Partículas (PSO)**: técnica de optimización.
- **Recocido Simulado (SA)**: método de simulación.
- **Evolución Diferencial (DE)**: aproximación evolutiva.

La elección final de los algoritmos y la función objetivo se llevará a cabo durante el desarrollo de esta fase.

### 5. Evaluación experimental

Se planificarán experimentos para contrastar los métodos seleccionados a través de:

- El valor de la función objetivo.
- La rapidez de convergencia.
- El tiempo necesario para la ejecución.
- La estabilidad entre ejecuciones independientes.
- La sensibilidad a los parámetros.
- La calidad de la clasificación, en caso de contar con datos de referencia apropiados.

Se incluirán métodos de referencia para proporcionar un contexto a los resultados.

---

## Estructura del repositorio

```text
Satellite-imagery-processing/
│
├── data/
│   ├── B02.tiff
│   ├── B03.tiff
│   ├── B04.tiff
│   ├── B08.tiff
│   └── B11.tiff
│
├── src/
│   ├── cargar_visualizar.m
│   └── [algoritmos de optimización pendientes]
│
├── results/
│   ├── rgb.png
│   ├── valid_mask.png
│   ├── ndvi.png
│   ├── ndwi.png
│   └── ndbi.png
│
├── docs/
│   └── [memoria y documentación pendientes]
│
├── main.m
└── README.md
```

*Esta estructura representa la organización prevista. Algunos archivos todavía están pendientes de incorporación.*

---

## Estado actual del proyecto

### Fase 1. Procesamiento inicial

- [x] Establecimiento del repositorio en GitHub.
- [x] Elección de la región de estudio.
- [x] Adquisición de imágenes de Sentinel-2.
- [x] Carga de las bandas en MATLAB.
- [x] Creación de la composición RGB.
- [x] Generación de la máscara de píxeles válidos.
- [x] Cálculo del NDVI.
- [x] Cálculo del NDWI.
- [x] Cálculo del NDBI.

### Fase 2. Optimización

- [ ] Definición matemática del problema.
- [ ] Elaboración de la función objetivo.
- [ ] Selección de algoritmos.
- [ ] Implementación en MATLAB.
- [ ] Realización de pruebas de funcionamiento.
      
### Fase 3. Experimentación

- [ ] Diseño de experimentos.
- [ ] Ejecución de los algoritmos.
- [ ] Comparación de resultados.
- [ ] Análisis estadístico.
- [ ] Estudio de la convergencia.
- [ ] Análisis de sensibilidad.

### Fase 4. Documentación

- [ ] Interpretación de resultados.
- [ ] Redacción de la memoria científica.
- [ ] Elaboración de conclusiones.
- [ ] Preparación de la presentación final.

---

## Herramientas y tecnologías

**MATLAB:** procesamiento de imágenes, operaciones matriciales, implementación de algoritmos y análisis experimental.

**Sentinel-2:** adquisición de información multiespectral de la superficie terrestre.

**Copernicus Data Space Ecosystem:** plataforma utilizada para obtener las imágenes satelitales.

**GitHub:** almacenamiento del código, control de versiones y colaboración entre integrantes del grupo.

---

## Resultados esperados

Para el término del proyecto se aspira a contar con:

1. Un protocolo replicable para el procesamiento de imágenes multiespectrales.
2. Visualizaciones gráficas de la área de investigación.
3. Cartografía de índices espectrales.
4. Una formulación matemática del desafío de optimización.
5. Implementaciones operativas de los algoritmos elegidos.
6. Resultados experimentales que sean comparables.
7. Un análisis sobre la eficiencia, calidad y restricciones de cada técnica.
8. Un informe científico que registre el trabajo llevado a cabo.

---

## Fuente de datos

Las imágenes utilizadas proceden de [Copernicus Data Space Ecosystem](https://dataspace.copernicus.eu/), plataforma de acceso a los datos del programa europeo Copernicus.

---

## Autores

Proyecto desarrollado en grupo para la asignatura de **Optimización y Heurística**.

**Grado en Ciencia e Ingeniería de Datos**

**Universidad de Las Palmas de Gran Canaria (ULPGC)**

Curso académico 2026–2027.

---

## Estado del proyecto

**En desarrollo — Octubre de 2026**

Actualmente se ha completado la etapa inicial de carga, visualización y cálculo de índices espectrales.

Las siguientes etapas se centrarán en formular el problema de optimización, implementar los algoritmos seleccionados y desarrollar el estudio experimental comparativo.
