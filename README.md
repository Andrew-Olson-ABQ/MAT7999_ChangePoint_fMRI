# Graph-Based Change-Point Detection in fMRI Data

**Master's Project — Wayne State University, Department of Mathematics**
**Andrew Olson | MA in Applied Mathematics | May 2026**

---

## Overview

This project applies and extends graph-based change-point detection methods to high-dimensional functional MRI (fMRI) data. The primary goal is to identify structural changes in multivariate brain activity signals and compare how classical and nonparametric statistical methods detect those changes.

A key contribution of this work is a **novel temporally weighted graph construction** that incorporates temporal proximity into the distance function used to build the similarity graph. This modification reveals alternative, more stable change-point locations that are not detected by the standard approach.

---

## Data

Data were obtained from the [OpenNeuro](https://openneuro.org) public repository:

> **Dataset:** Targeting Dynamic Facial Processing Mechanisms in Superior Temporal Sulcus Using fMRI Neurofeedback  
> **DOI:** [10.18112/openneuro.ds004141.v1.0.5](https://doi.org/10.18112/openneuro.ds004141.v1.0.5)

- 20 participants (experimental and sham groups)
- 4D BOLD NIfTI scans: `64 × 64 × 33 × 300` per functional run
- After masking: `300 time points × 133,056 voxel signals`
- Dimensionality reduced to `300 × 30` via PCA prior to change-point analysis

---

## Methods

Five change-point detection methods are applied to the PCA-reduced multivariate time series:

| Method | Type |
|---|---|
| PELT | Classical, penalized cost minimization |
| Binary Segmentation | Classical, recursive splitting |
| Energy / Rank-Energy | Nonparametric, distributional |
| Kernel / MMD | Nonparametric, kernel-based |
| Graph-Based Edge Count | Nonparametric, graph-based (MST) |

### Temporally Weighted Graph (Novel Contribution)

The standard graph-based method constructs a minimum spanning tree using feature similarity alone. This project introduces a modified distance:

```
d_λ(i, j) = (1 − cor(Z_i, Z_j)) + λ · |i − j| / n
```

where `λ ≥ 0` controls the influence of temporal proximity. When `λ = 0` the method reduces to the original. As `λ` increases, the graph favors connections between temporally adjacent observations, revealing structurally different segmentations.

---

## Key Results

- **Standard graph-based method** detects a highly significant change point at **τ̂ = 267** (p ≈ 0)
- **Temporally weighted method** (λ ≥ 0.5) shifts the estimate to **τ̂ ≈ 154**, which also remains highly significant
- Multi-subject analysis (n = 20) shows the temporally weighted method produces **more concentrated and consistent** change-point estimates across subjects (paired t-test: t = 2.00, p = 0.047)

---

## Repository Contents

```
├── MAT7999_ChangePoint_fMRI.R          # Full R analysis script
├── MAT7999_ChangePoint_fMRI.pdf        # Formal project report (LaTeX/Overleaf)
├── MAT7999_ChangePoint_fMRI_Presentation.pptx  # Project presentation
└── README.md
```

---

## Requirements

All analysis was performed in **R**. Key packages used:

- `oro.nifti` — reading NIfTI fMRI files
- `changepoint` — PELT and binary segmentation
- `ecp` — energy-based change-point detection
- `ggplot2` — visualization

---

## How to Run

1. Download the fMRI dataset from [OpenNeuro ds004141](https://openneuro.org/datasets/ds004141)
2. Update the file path to the NIfTI `.nii.gz` file in the R script
3. Run `MAT7999_ChangePoint_fMRI.R` in order — the script handles masking, PCA, and all five detection methods

---

## Author

**Andrew Olson**  
MA, Applied Mathematics — Wayne State University  
BS, Applied Mathematics (Minors: Physics, Astronomy) — University of Arizona  
📧 hw0460@wayne.edu
