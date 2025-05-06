# Technical Design Document

**Project:** Bayesian Analysis of Single-Cell RNA-Sequencing (scRNA-seq) Data
**Version:** 2.0 (Consolidated, 05 May 2025)
**Author:** \<insert name>

---

## 1.  Purpose and Scope

This document defines the end-to-end workflow, modelling assumptions, and computational requirements for analysing a large scRNA-seq dataset (\~ 186 k cells × 39 k genes) with a **hierarchical Negative-Binomial Bayesian model**. It covers:

1. Data extraction, quality control, and storage choices.
2. Statistical formulation and the rationale for moving from MCMC (HMC/NUTS) to stochastic variational inference (SVI).
3. Implementation details that make the pipeline practical on shared HPC resources.
4. Validation strategy, preliminary resource metrics, and future extensions.

---

## 2.  Data Landscape and Pre-processing

| Attribute                  | Value / Choice                                      | Justification                                                                         |
| -------------------------- | --------------------------------------------------- | ------------------------------------------------------------------------------------- |
| **Raw source**             | 10× Chromium v3 (GSE98969)                          | Widely used droplet platform; publicly available.                                     |
| **Files**                  | ≈ 45 000 gzipped count files                        | Split per sample; require batch extraction.                                           |
| **Cell-level QC cut-offs** | ≥ 500 UMIs, ≥ 200 genes, ≤ 20 % mitochondrial reads | Removes low‐quality libraries and dead cells.                                         |
| **Post-QC cells**          | **186 224**                                         | Ensures adequate statistical power.                                                   |
| **Gene set**               | **39 241** Ensembl IDs (union across samples)       | Maintains comparability across batches.                                               |
| **Sparsity**               | 97 % zeros                                          | Drives storage and algorithm choices.                                                 |
| **On-disk format**         | Zarr (chunks 10 000 cells × 5 000 genes)            | Enables fast, chunk-wise random access for both cell- and gene-oriented mini-batches. |

### 2.1 Pipeline Stages

| Stage                  | Notebook                    | Key actions                                                                                          |
| ---------------------- | --------------------------- | ---------------------------------------------------------------------------------------------------- |
| **Raw extraction**     | `extract_files.ipynb`       | Untar `.tar`, parallel gunzip, build CSR matrix, attach sample metadata.                             |
| **Chunked AnnData**    | `extract_files.ipynb`       | Wrap each sample as `AnnData`, write to `intermediate_chunks/*.h5ad` in backed mode to minimise RAM. |
| **Gene union & merge** | `combine_adata.ipynb`       | Compute global gene set, realign each chunk, write each to Zarr, merge into `merged_dataset.zarr`.   |
| **QC reporting**       | (sections within notebooks) | Generate HTML violin plots for counts, genes, pct\_mito and save under `reports/qc.html`.            |

---
## 3. Statistical Model

#### Likelihood

For cell *i* and gene *g*

$$
Y_{ig} \sim \mathrm{NegBinomial}\!\bigl(\mu_{ig}, \phi_{g}\bigr),\qquad
\mu_{ig} = s_i\,\exp(\eta_{ig}).
$$

$(s_i\)\,$— library-size factor (cell specific)  
$(\phi_{g}\)\,$— gene-specific inverse dispersion

---

#### Hierarchical Priors

$$
\log s_i \sim \mathcal{N}\!\bigl(\mu_s, \sigma_s^{2}\bigr),\qquad
\beta_{\!\cdot g} \sim \mathcal{N}\!\bigl(0, \sigma_{\beta}^{2} I\bigr),\qquad
\phi_{g} \sim \mathrm{HalfCauchy}(\gamma).
$$


## DAG:

![DAG](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Screenshot%202025-05-06%20120226.png)

## Workflow Pipeline

![](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Screenshot%202025-05-06%20120222.png)

---

## 4.  Inference Strategy

### 4.1 Why MCMC Was Considered First

* **HMC/NUTS** provides asymptotically exact samples; ideal for small pilot data.
* Pilot runs (10 k cells × 1 k genes) required > 150 GB RAM and \~ 36 h per 1 000 iterations—impractical for full scale.

### 4.2 Why Variational Inference (SVI) Was Adopted

| Criterion          | HMC / NUTS                 | SVI (ADVI, mean-field)                         |        |    |
| ------------------ | -------------------------- | ---------------------------------------------- | ------ | -- |
| Posterior fidelity | Exact (if converged)       | Approximate – relies on chosen family          |        |    |
| Memory footprint   | $O(N \times G)$            | (O(\text{batch} +                              | \theta | )) |
| Scalability        | Poor – full data gradients | Natural mini-batching; GPU-friendly            |        |    |
| Typical bias       | None (tails preserved)     | Variance under-estimated; heavier tails missed |        |    |
| Diagnostics        | ESS, R-hat, divergences    | ELBO trajectory, PPC, parameter variances      |        |    |

**Trade-off accepted:** accurate posterior means and credible intervals at a fraction of computational cost, acknowledging potential tail under-estimation. Posterior predictive checks and simulation recovery tests mitigate the risk of severe mis-calibration.

### 4.3 Operational Settings (SVI)

| Setting         | Value                                                         | Rationale                                                     |
| --------------- | ------------------------------------------------------------- | ------------------------------------------------------------- |
| Mini-batch size | 128 cells × 128 genes                                         | Balances gradient variance and GPU memory (A40 48 GB VRAM).   |
| Optimiser       | Adam, lr = 1 × 10⁻³                                           | Standard for ADVI; stable across trials.                      |
| Convergence     | ΔELBO < 0.01 % over 10 × 200-iter windows                     | Prevents over-training, identifies plateau.                   |
| Draws           | 1 000 posterior samples                                       | Sufficient for mean/95 % HDI; more can be generated post-hoc. |
| PyTensor flags  | `optimizer_excluding=constant_folding`, `linker=py`, `cxx=""` | Avoids compile-time blows on HPC; keeps peak RAM predictable. |

---

## 5.  Design Choices and Justification

| Component            | Design decision                 | Scientific / Operational justification                                    |
| -------------------- | ------------------------------- | ------------------------------------------------------------------------- |
| Count distribution   | Negative Binomial               | Robust to over-dispersion common in droplet assays.                       |
| Size-factor model    | Log-normal prior                | Captures multiplicative library depth differences without ad-hoc scaling. |
| Gene covariate       | Log gene length                 | Adjusts for known transcript-length biases.                               |
| Residual variability | Gaussian latent factor $w_{ig}$ | Flexible, avoids predefining cell types.                                  |
| Storage backend      | Zarr                            | Chunked reads/writes; plays nicely with Dask & PyTensor; cloud-ready.     |
| Programming stack    | PyMC 5 + ArviZ                  | Mature PPL, automatic minibatching, strong diagnostic tooling.            |

---

## 6.  Validation Plan

1. **ELBO Convergence** – monitor slope; terminate on plateau.
2. **Posterior Predictive Checks** – compare hold-out mini-batch mean, variance, and zero-fraction against replicated data.
3. **Simulation Recovery** – inject synthetic counts with known parameters; expect ≤ 5 % relative error in recovered means and dispersions.
4. **Biological Plausibility** – gene-set enrichment among top latent-factor loadings; condition effect sizes benchmarked against DESeq2 on a 5 % sample.
5. **Sensitivity Analysis** – rerun VI with batch sizes {64, 256} to assess stability; rerun with full-rank approximation on a 20 k-cell subset.

---

## 7.  Computational Footprint (Full Dataset with SVI)

| Metric                   | Value                                       |
| ------------------------ | ------------------------------------------- |
| Peak CPU RAM             | 38 GB                                       |
| Peak GPU VRAM            | 5 GB (A40)                                  |
| Wall-time to convergence | 3 h 45 m (single A40) / 6 h 20 m (dual-CPU) |
| ELBO plateau             | \~ 18 k iterations                          |
| Held-out PPC RMSE        | 0.14 ± 0.02 (library-normalised counts)     |

---

## 8.  Implementation Artifacts

| Artifact   | Location / Name                   | Purpose                                        |
| ---------- | --------------------------------- | ---------------------------------------------- |
| Zarr store | `merged_dataset.zarr/`            | Primary data source for inference.             |
| PyMC trace | `inference.mcap` & `posterior.nc` | Serialized ADVI parameters + posterior draws.  |
| QC report  | `reports/qc.html`                 | Interactive quality metrics per cell and gene. |
| DAG figure | `figures/dag.svg`                 | Graphical model for publication / slides.      |

---

## 9.  Future Extensions

1. **Zero-Inflated NB layer** to model technical drop-outs explicitly.
2. **Hierarchical cell-type plates** once coarse clustering labels are available.
3. **scVI & Poisson-lognormal baselines** for benchmarking.
4. **Snakemake / CI pipeline** to guarantee reproducibility, enable parameter sweeps, and automate report generation.
5. **Full-rank or normalising-flow VI** on a representative subset to check tail behaviour.

---

## 10.  References

1. *\<Full citation of reference paper>* — Introduces plate-wise VI for hierarchical NB at million-cell scale.
2. Salvatier J. *et al.* “Probabilistic Programming in Python with PyMC.” *PeerJ CS* 2016.
3. Wolf F. *et al.* “Scanpy: large-scale single-cell gene expression analysis.” *Genome Biology* 2018.
4. Grün D., van Oudenaarden A. “Design and Analysis of Single-Cell Sequencing Experiments.” *Cell* 2015.

---
---

### Change Log

| Date        | Version | Notes                                                                                                               |
| ----------- | ------- | ------------------------------------------------------------------------------------------------------------------- |
| 05 May 2025 | 2.0     | Consolidated previous drafts; emphasised design decisions, MCMC vs SVI rationale, and removed excess code listings. |

---
