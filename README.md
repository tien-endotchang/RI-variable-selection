# Variable Selection Using Relative Importance Ranking

This repository contains the codes for the paper:

**Tien-En Chang & Argon Chen (2025). Variable Selection Using Relative Importance Ranking.**  [Link to paper coming soon]

We propose a class of filter-based variable selection methods built on **Relative Importance (RI)** measures, including General Dominance (GD), Comprehensive Relative Importance (CRI), and a new computationally efficient variant, **CRI.Z**. Our methods are robust to multicollinearity and competitive with state-of-the-art approaches like the lasso and relaxed lasso.

---

## Installation

Make sure you have the following installed:
*  R (>= 4.0.0)
*  R packages: `glmnet`, `relaimpo`, `care`, `ggplot2`, `ggh4x` etc.

    You can install all required packages with:
    ```r
    install.packages(c("glmnet", "relaimpo", "care", "ggplot2", "ggh4x"))
    ```
* Make sure to set your working directory `~/RI-variable-selection/`.

---

## Run Simulations

Run the following R scripts from the `main/` directory.

### 1. Part 1: Variable Ranking (Section 4.2 of paper)
```r
source("main/part1_ranking/sim.lo.select.R")
source("main/part1_ranking/sim.med.select.R")
source("main/part1_ranking/sim.hi50.select.R")
source("main/part1_ranking/sim.hi100.select.R")
```

### 2. Part 2: Modeling (Section 4.3 of paper)
```r
source("main/part2_modeling/sim.lo.R")
source("main/part2_modeling/sim.med.R")
source("main/part2_modeling/sim.hi50.R")
source("main/part2_modeling/sim.hi100.R")
```

## Results
Running the above scripts will create a `results/` folder containing:
* `fig/`: all generated simulation plots (`.pdf`)
* `res.`: `.RDS` files containing raw simulation ouputs

---