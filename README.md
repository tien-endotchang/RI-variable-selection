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

---

## Reproducing Results

The experiments are divided into three parts matching the sections of the manuscript. Run these scripts from the `main/` directory.

### 1. Variable Ranking (Section 4.2)
Simulations comparing the ranking accuracy of different RI measures.

| Script (`main/part1_ranking/`) | Output Location | Key Figures Generated |
| :--- | :--- | :--- |
| `sim.lo.select.R` | `results/rds/part1/lo` | Fig 1 (`sim.lo.S.sub`) |
| `sim.med.select.R` | `results/rds/part1/med` | Fig 2 (`sim.lo.Pr.sub`) |
| `sim.hi50.select.R` | `results/rds/part1/hi50` | (Supplementary Figures) |
| `sim.hi100.select.R` | `results/rds/part1/hi100` | Fig 3-4 (`sim.hi100.S/Pr.sub`) |

## Run Simulations

Run the following R scripts from the `main/` directory.

### 1. Part 1: Variable Ranking (Section 4.2 of paper)
```r
source("main/part1_ranking/sim.lo.select.R")
source("main/part1_ranking/sim.med.select.R")
source("main/part1_ranking/sim.hi50.select.R")
source("main/part1_ranking/sim.hi100.select.R")
```
These will generate .rds simulation files in `results/rds/part1/{lo,med,hi50,hi100}` respectively and .pdf figures `results/fig/part1`. Specifically, Figures 1-4 (`sim.lo.S.sub`, `sim.lo.Pr.sub`, `sim.hi100.S.sub` and `sim.hi100.Pr.sub`) are stored in `results/fig/part1`. The rest of figures are shown in Supplementary Material.

### 2. Part 2: Modeling (Section 4.3 of paper)
```r
source("main/part2_modeling/sim.lo.R")
source("main/part2_modeling/sim.med.R")
source("main/part2_modeling/sim.hi50.R")
source("main/part2_modeling/sim.hi100.R")
```
These will generate .rds simulation files in `results/rds/part2/{lo,med,hi50,hi100}` respectively and .pdf figures `results/fig/part2`. Specifically, Figures 5-8 (`sim.n100.p10.val.F.sub`, `sim.n100.p10.val.err.rel.sub`, `sim.n100.p1000.val.F.sub` and `sim.n100.p1000.val.err.rel.sub`) are stored in `results/fig/part2`. The rest of figures are shown in Supplementary Material.

### 3. Real-world Dataset Example (Section 5 of paper)
```r
source("main/real-world/real-world.R")
```
This will use real-world datasets in `data/{aml,gli_85}` to run variable selection experiments. These datasets are derived (from `.mat` to `.csv`) from the leukemia dataset of Golub et al. (Science, 1999) and the glioma dataset of Freije et al. (Cancer Research, 2004). The original data were obtained from the public scikit-feature feature selection repository (Arizona State University). This will generate tables and experiment files `results/tab/res.{aml,gli_85}.{csv,rds}, resulting in Tables 2 and 3. 

## Results
Running the above scripts will create a `results/` folder containing:
* `fig/`: all generated simulation plots (`.pdf`)
* `rds/`: `.RDS` files containing raw simulation ouputs
* `tab/`: table results (running time, real-world examples)

<<<<<<< HEAD
the summarized running time of part 1 and 2 can be computed by running `main/runtime.R`.
---
=======
---
>>>>>>> 2dfc1888bd435ac77f45e79ade9ba11b2cdf4fe7
