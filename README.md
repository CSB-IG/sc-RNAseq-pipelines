# scRNA-seq Hackathon 

## ** Objetivo**
Desarrollar pipelines óptimos de análisis de single-cell RNA-seq en R, Python y un enfoque híbrido (R + Python), asegurando modularidad, documentación clara y benchmarks de rendimiento.

---

## ** Estructura del Repositorio**
```
├── data/          # Datos de prueba
├── scripts/       # Código de cada pipeline
│   ├── R/
│   ├── Python/
│   ├── Hybrid/
├── docs/          # Documentación y tutoriales
├── results/       # Resultados y benchmarks
├── benchmarks/    # Evaluaciones de rendimiento
├── .github/       # Configuración de issues y workflows
├── README.md      # Descripción del proyecto
├── LICENSE        # Licencia (MIT, Apache, etc.)
└── CONTRIBUTING.md # Reglas de contribución
```

---

## ** ¿Cómo contribuir?**
1. **Clona el repositorio**
   ```bash
   git clone https://github.com/CSB-IG/sc-RNAseq-pipelines.git
   ```
2. **Crea una rama para tu contribución**
   ```bash
   git checkout -b feature-nueva
   ```
3. **Sube tus cambios y haz un pull request (PR)**
   ```bash
   git add .
   git commit -m "Descripción clara del cambio"
   git push origin feature-nueva
   ```
4. **Espera revisión y aprobación** 

---

## **📌 Milestones e Issues**
### 🔹 **Fases del Hackatón**
✅ **Set-up inicial**: Configuración del repositorio y datos de prueba.
✅ **Implementación del pipeline**: Desarrollo en R, Python y versión híbrida.
✅ **Pruebas y optimización**: Evaluación de performance y validaciones.
✅ **Documentación y tutorial**: Creación de tutoriales en Jupyter/RMarkdown.
✅ **Benchmarks y paralelización**: Medición de tiempos y consumo de recursos.

###  **Principales Issues**
- Implementar Preprocesamiento (QC y Alineamiento) 
- Implementar Filtrado y Normalización 
- Implementar Clustering y Reducción de Dimensionalidad 
- Implementar Identificación de Tipos Celulares 
- Implementar Comparación Casos vs Controles 
- Implementar RNA Velocity y Pseudotiempo 
- Implementar Interacciones Célula-Célula 
- Realizar Benchmarks de Tiempos y Memoria 
- Identificar Pasos Paralelizables 
- Crear Documentación y Tutorial 

---

## **📖 Documentación y Tutorial**
Cada pipeline tendrá:
- Documentación clara en `docs/`
- Notebooks de ejemplo (`tutorial.ipynb` en Python, `tutorial.Rmd` en R)
- Explicaciones paso a paso para ejecutar cada script

---

## **📊 Evaluación de Rendimiento**
Cada implementación será evaluada en:
- **Tiempo de ejecución** (segundos)
- **Consumo de memoria RAM** (MB)
- **Uso de almacenamiento en disco** (MB)

---

## ** Licencia**
Este proyecto se distribuye bajo la licencia MIT. Share and Enjoy. 

