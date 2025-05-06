### `Bayesian Inference for Anndata.ipynb` – Detailed Documentation

---

#### 1  Purpose

This notebook is a **proof-of-concept implementation of stochastic variational inference (SVI) for a hierarchical Negative-Binomial model** on a **subset** of the merged microglia dataset. It demonstrates:

* How to stream counts from an `AnnData` object into PyMC minibatches
* How to build a non-centred NB model that corrects for sequencing depth and gene length
* How to fit the model with Automatic Differentiation Variational Inference (ADVI) over 300 000 iterations
* How to inspect convergence (ELBO) and obtain an approximate posterior trace for downstream exploratory analysis ([GitHub][1])

---

#### 2  Input Artefacts

| File                  | Description                                          | Produced by           |
| --------------------- | ---------------------------------------------------- | --------------------- |
| `merged_dataset.h5ad` | 186 224 cells × 39 241 genes – merged AnnData object | `combine_adata.ipynb` |

> Only a **tiny slice** (1 000 cells × 500 genes) is sampled for this demo to keep compilation and ADVI tractable on CPU.

---

#### 3  Software Requirements

| Package              | Version (tested)         | Comment                                                                    |
| -------------------- | ------------------------ | -------------------------------------------------------------------------- |
| **PyMC**             | 5.12                     | Probabilistic programming                                                  |
| **ArviZ**            | ≥ 0.17                   | Diagnostics                                                                |
| **Scanpy / AnnData** | ≥ 1.10                   | Data loading                                                               |
| **PyTensor**         | auto-installed with PyMC | Config overrides disable C-compilation for HPC clusters without a compiler |

All dependencies are already listed in `envs/pymc.yml`.

---

#### 4  Data Pre-processing Steps

1. **Load** `merged_dataset.h5ad` with Scanpy.

2. **Sub-select** 1 000 random cells and 500 random genes for faster iteration:

   ```python
   adata_sub = adata[:1000, :500].copy()
   counts = adata_sub.X.toarray().T       # shape → (n_genes, n_cells)
   ```

3. **Generate or impute** auxiliary covariates (if absent):

   * `gene_lengths` – uniform 1 000-3 000 bp
   * `seq_depth` – median of `Bases` per cell (falls back to 10 000 if missing)
     These are injected as offsets in the log-linear predictor.

---

#### 5  Model Formulation (minibatch version)

For each gene *g* and cell *i*:

$$
\begin{aligned}
\log(\mu_{ig}) &= \underbrace{\bigl(\mu_{\text{gene}} + \sigma_{\text{gene}} z_g\bigr)}_{\text{gene-level}} \;+\;  
                \underbrace{\bigl(\mu_{\text{cond}} + \sigma_{\text{cond}} z_i\bigr)}_{\text{cell-level}} \;+\;
                \log\!\frac{s}{s_\mathrm{ref}} \;+\;
                \log\!\frac{L_{g}}{L_\mathrm{ref}} \\[4pt]
y_{ig} &\sim \text{NB}\!\bigl(\mu_{ig}, \alpha\bigr)
\end{aligned}
$$

* **Gene effects** $z_g\sim\mathcal N(0,1)$
* **Cell (condition) effects** $z_i\sim\mathcal N(0,1)$
* **Dispersion** $\alpha\sim\text{Half-Cauchy}(2)$

A non-centred parameterisation for both gene- and cell-hierarchies improves optimisation.
The entire flattened count matrix is streamed through **`pm.Minibatch`** arrays (batch = 128) so that ELBO is computed on mini-subsets each gradient step ([GitHub][1]).

---

#### 6  Inference Settings

```python
with pm.Model() as model_vi:
    ...
    approx = pm.fit(
        n=300_000,                 # iterations
        method="advi",
        callbacks=[pm.callbacks.CheckParametersConvergence(tolerance=0.01)],
        progressbar=True
    )
    trace_vi = approx.sample(draws=1_000)
```

* **Optimizer** – default ADAM (PyMC)
* **Learning rate** – automatically tuned
* **Convergence check** – terminates early if the running SD of parameters stabilises within 0.01

---

#### 7  Notebook Workflow

| Stage                  | Cells | Description                                                             |
| ---------------------- | ----- | ----------------------------------------------------------------------- |
| Imports & env-flags    | 1-3   | Disable PyTensor C-compiler (use `linker=py`) for portability           |
| Read data / slice      | 4-6   | Load `AnnData`, take 1 000 × 500 subset, cast to dense `np.ndarray`     |
| Covariate generation   | 7-10  | Fill missing sequencing depth and gene lengths                          |
| Model build & ADVI fit | 11-12 | Define hierarchical NB + offsets, run 300 k iterations with minibatches |
| Posterior summary      | 13    | `az.summary(trace_vi)`                                                  |
| Convergence plot       | 14    | Plot ELBO history (`approx.hist`)                                       |

Total wall-clock ≈ 40–60 min on a modern laptop CPU (or <10 min on a single GPU if JAX backend is configured).

---

#### 8  Outputs

| File                   | Contents                                 |
| ---------------------- | ---------------------------------------- |
| *In-memory* `trace_vi` | 1 000 posterior draws for all parameters |
| *Optional*             | ELBO curve figure displayed inline       |

No artefacts are written to disk; export with:

```python
az.to_netcdf(trace_vi, "svi_nb_trace.nc")
```

---

#### 9  Interpretation & Next Steps

* Focus here is **methodological feasibility**: demonstrating that an NB model with tens of thousands of parameters can be fitted by SVI on scRNA-seq counts.
* For production‐scale runs:

  * **Increase** gene/cell subsets, potentially switch to **black-box VI on GPU** (JAX).
  * **Incorporate** real covariates (cell cycle, batch, treatment) instead of simulated ones.
  * **Compare** ADVI to NUTS on a smaller gene set to assess posterior quality.

---

#### 10  Troubleshooting

| Issue                            | Likely Cause                      | Remedy                                                       |
| -------------------------------- | --------------------------------- | ------------------------------------------------------------ |
| ELBO divergence / no convergence | Learning rate too high            | Set `pm.fit(..., obj_optimizer=pm.adam(learning_rate=1e-4))` |
| `MemoryError` when flattening    | Dense conversion too big          | Keep matrix in sparse COO and sample indices instead         |
| Slow on CPU                      | PyTensor using pure-python linker | Install C-compiler or enable JAX GPU backend                 |

---

[1]: https://raw.githubusercontent.com/HarshiniR4/SCRNAseq_Bayesian_Inference/main/src/Bayesian%20Inference%20for%20Anndata.ipynb "raw.githubusercontent.com"
