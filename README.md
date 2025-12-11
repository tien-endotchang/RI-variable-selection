# Variable Selection Using Relative Importance Ranking

This repository contains the codes and data for the paper:

**Tien-En Chang & Argon Chen (2025). Variable Selection Using Relative Importance Ranking.**  

We propose a class of filter-based variable selection methods built on **Relative Importance (RI)** measures, including General Dominance (GD), Comprehensive Relative Importance (CRI), and a *new*, *computationally efficient* variant, **CRI.Z**. Our methods are robust to multicollinearity and competitive with state-of-the-art approaches like the lasso and relaxed lasso.

---

## Installation & Setup

### Prerequisites
Make sure you have the following installed:
*  **R version:** 4.4.1 (2024-06-14 ucrt)
*  **Dependencies:** `glmnet`, `relaimpo`, `care`, `ggplot2`, `ggh4x` etc.

    You can install all required packages with:
    ```r
    install.packages(c("glmnet", "relaimpo", "care", "ggplot2", "ggh4x"))
    ```
* Make sure to set your working directory `~/RI-variable-selection/`.

### Setup
1. Clone or download this repository.
2. Open R/RStudio.
3. **Important:** Set your working directory to the project root:
   ```r
   setwd("~/RI-variable-selection/")
   ```

---

## Reproduction Instructions

All analysis scripts are located in the `main/` directory. Results are automatically saved to the `results/` directory.

### 1. Part 1: Variable Ranking (Paper Section 4.2)
These scripts simulate the ranking capabilities of the proposed methods across low, medium, and high problem dimension settings.

| Simulation Setting | Script Path | Output Data (`.rds`) | Output Figure (`.pdf`) |
| :--- | :--- | :--- | :--- |
| **low** | `main/part1_ranking/sim.lo.select.R` | `results/rds/part1/lo/` | Figs 1 & 2 |
| **medium** | `main/part1_ranking/sim.med.select.R` | `results/rds/part1/med/` | *(See Supp. Material)* |
| **high-50** | `main/part1_ranking/sim.hi50.select.R` | `results/rds/part1/hi50/` | *(See Supp. Material)* |
| **high-100** | `main/part1_ranking/sim.hi100.select.R` | `results/rds/part1/hi100/` | Figs 3 & 4 |

To run all ranking simulations:
```r
source("main/part1_ranking/sim.lo.select.R")
source("main/part1_ranking/sim.med.select.R")
source("main/part1_ranking/sim.hi50.select.R")
source("main/part1_ranking/sim.hi100.select.R")
```

### 2. Part 2: Modeling (Paper Section 4.3)
These scripts evaluate the predictive performance of the selected models.

| Simulation Setting | Script Path | Output Data (`.rds`) | Output Figure (`.pdf`) |
| :--- | :--- | :--- | :--- |
| **low** | `main/part2_modeling/sim.lo.R` | `results/rds/part2/lo/` | Figs 5 & 6 |
| **medium** | `main/part2_modeling/sim.med.R` | `results/rds/part2/med/` | *(See Supp. Material)* |
| **high-50** | `main/part2_modeling/sim.hi50.R` | `results/rds/part2/hi50/` | *(See Supp. Material)* |
| **high-100** | `main/part2_modeling/sim.hi100.R` | `results/rds/part2/hi100/` | Figs 7 & 8 |

To run all modeling simulations:
```r
source("main/part2_modeling/sim.lo.R")
source("main/part2_modeling/sim.med.R")
source("main/part2_modeling/sim.hi50.R")
source("main/part2_modeling/sim.hi100.R")
```

#### Runtime Analysis
To generate the summarized running time for Part 1 and Part 2 simulations:
```r
source("main/runtime.R")
```
**Output:** Tables 2 and 3 in the paper (saved to `results/tab/`).

### 3. Real-world Examples (Paper Section 5)
We utilize two benchmark gene expression datasets:
1.  **Leukemia:** (n,p)=(72,7129) from Golub et al. (Science, 1999)
2.  **Glioma:** (n,p)=(85,22283) from Freije et al. (Cancer Research, 2004)

*Original data obtained from the [Arizona State University Feature Selection Repository](http://featureselection.asu.edu/).*

Run the analysis:
```r
source("main/real-world/real-world.R")
```
**Output:** Tables 4 and 5 in the paper (saved to `results/tab/`).

---

## Repository Structure

```text
RI-variable-selection/
├── data/                    # Datasets (.csv)
│   ├── aml/                 # Leukemia dataset
│   └── gli_85/              # Glioma dataset
├── external/                # Code adapted from HTT (2020)
├── main/                    # Source code
│   ├── part1_ranking/       # Section 4.2 simulations
│   ├── part2_modeling/      # Section 4.3 simulations
│   ├── real-world/          # Section 5 real-world dataset examples
│   ├── fig.df_modified      # Generate ECDF (fig 9)
│   └── runtime.R            # Script to compute summary runtimes
├── R/                       # Helper functions                               
├── results/                 # Generated outputs
│   ├── fig/                 # Plots (.pdf)
│   ├── rds/                 # Raw simulation objects (.rds)
│   └── tab/                 # Tables and runtime logs
└── README.md
```

---