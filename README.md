
# STERapp Stereology Pipeline  
Modular and reproducible pipeline for ovarian data analysis and stereological fecundity estimation in fish

This repository contains the scripts and documentation associated with a complete workflow for analyzing ovarian data using stereological principles. The methodology integrates histological measurements, image analysis with STERapp, statistical modeling in R, and final calculation of potential and real fecundity for each specimen.

The goal is to provide a **reproducible, traceable, and automated** tool that minimizes human error in labor‑intensive tasks such as oocyte classification, area and density calculations, or ovary volume estimation.

---

## 📂 Repository Structure

```
STERapp-stereology-pipeline/
├── src/
│   └── V1_stereology_template_for_STERapp.R
├── documentation/
│   ├── Script_Documentation.docx
│   └── Script_Functionality_Overview.docx
└── README.md
```

---

## 🧬 Pipeline Overview

The main R script automates the entire analytical process:

### **1. Loading packages and setting project paths**
Initializes the working directories (`input_data/` and `outputs/`) and checks that all required R packages are installed.

### **2. Data import and cleaning**
The script imports three key datasets:
- `diameters_spp.csv`: oocyte and nucleus diameters  
- `ovaries_data.csv`: fresh and fixed ovary weights  
- `sterapp_data.csv`: stereology data exported from STERapp  

It also corrects common formatting issues such as:
- extra delimiters  
- decimal separators  
- inconsistent labels  
- trailing characters  

### **3. Oocyte diameter class creation**
Automatic generation of diameter classes and calculation of:
- frequency per class  
- mean nucleus diameter  
- linear model predicting nucleus diameter from oocyte diameter  

### **4. Frequency correction**
Model-based correction factors are applied to adjust oocyte frequency distributions.

### **5. Stereology data processing**
Includes:
- file name normalization  
- extraction of fish ID and image ID  
- developmental stage filtering and recoding  
- calculation of total area, stage area, object counts, and densities per fish  

### **6. Ovary volume calculation**
Based on:
- linear modeling of fixed vs. fresh ovary weight  
- Scherle coefficient (estimated or default)  
- optional histological shrinkage correction  

### **7. Stereological parameter calculation**
The script computes:
- β coefficient (axis ratio)  
- M1 and M3 moments of the diameter distribution  
- K coefficient  

### **8. Fecundity estimation**
Two metrics are produced:
- **Potential fecundity:** sum of all oocytes estimated  
- **Real fecundity:** viable oocytes minus atretic ones  
Results are exported to `fecundidad_spp.xlsx`.

### **9. Visualization**
A bar plot comparing potential vs. real fecundity per fish is created and saved as `fecundidad_plot.png`.

---

## 📥 Required Input Files

All files must be placed inside the folder `input_data/`:

| File | Required columns | Purpose |
|------|------------------|---------|
| `diameters_spp.csv` | `fish ID`, `diam_ooc`, `diam_nuc` | Oocyte sizing and diameter modeling |
| `ovaries_data.csv` | `fish_ID`, `ovary_fresh_weight`, `ovary_fixed_weight` | Ovary biometrics and volume estimation |
| `sterapp_data.csv` | STERapp full export | Image metadata, areas, developmental stages, counts |

---

## ▶️ Running the Script

1. Place all CSV input files into `input_data/`
2. Open R or RStudio
3. Run:

```r
source("src/V1_stereology_template_for_STERapp.R")
```

Outputs will appear in the `outputs/` directory:

- `diameters_calc.csv`  
- `stages_filtrados.csv`  
- `df_sterapp.csv`  
- `fecundidad_spp.xlsx`  
- `fecundidad_plot.png`  

---

## 📈 Outputs

The pipeline generates:

- automatic estimation of potential and real fecundity  
- comparison plots per fish  
- cleaned and standardized intermediate datasets  
- reproducible workflow with informative processing messages  

---

## 🧪 Requirements

- R ≥ 4.0  
- Required packages:  
  `dplyr`, `tidyr`, `readr`, `ggplot2`, `writexl`, `here`, `tools`, `scales`

---

## 📄 Additional Documentation

Located in the `documentation/` folder:

- **Script_Documentation.docx**  
  Full methodological explanation, stereological formulas, data flow.

- **Script_Functionality_Overview.docx**  
  Conceptual overview of the automation, required data, and preprocessing considerations.

---

## ✨ Author

**Sonia Rábade Uberos**  
Design of the pipeline, workflow standardization, and implementation of the stereological methodology.

---

## 📜 License

You may add a license (MIT, GPL‑3, etc.) depending on how you wish to share the project.
