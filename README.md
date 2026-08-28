# SCINET

**Reconstruction of cell-type-specific interactomes at single-cell resolution.**

Reference interactomes such as PCNet or STRING are aggregates: one network
standing in for every cell type in the body. But a gene pair that interacts in a
T cell need not interact in a neuron, and the wiring that matters for a disease
is usually the wiring in one particular cell population.

SCINET takes a single-cell expression profile and a reference interactome and
returns a separate, weighted network per cell type or cell state. Rather than
inferring networks from scratch, which single-cell data are too sparse to
support, it rescores the edges the reference already proposes:

1. **Gene activity.** Expression is transformed into activity scores by a
   rank-based inverse normal transform, optionally corrected for each gene's
   baseline expression across the dataset.
2. **Subsampling.** Cells within a population are repeatedly subsampled, so each
   edge receives a *distribution* of scores rather than a single point estimate.
   That distribution is what makes the difference between populations testable.
3. **Edge rescoring.** Each edge of the reference network is scored from the
   activities of its two endpoints, giving one weighted network per subsample.

Applied to disease gene sets, the resulting networks show that
disease-associated perturbations are modular in a cell-type-specific way: the
genes implicated in a disorder cluster tightly in some cell types and not in
others.

## Requirements

| Dependency | Notes |
| --- | --- |
| R >= 4.0 | |
| A C++17 compiler | GCC 7+ or Clang 6+, with OpenMP for multithreading |
| Rcpp, RcppArmadillo | Installed automatically |
| Matrix, igraph, methods | Installed automatically |

On Linux you may also want a tuned BLAS:

```bash
sudo apt-get install libopenblas-dev
```

## Install

```r
install.packages("devtools")
devtools::install_github("shmohammadi86/SCINET")
```

Or from a clone:

```bash
git clone https://github.com/shmohammadi86/SCINET.git
R CMD INSTALL SCINET
```

SCINET is usable on its own with any imputed expression matrix and any reference
network. Pairing it with [ACTIONet](https://github.com/shmohammadi86/ACTIONet)
additionally gives you archetype-decoupled activity scores, which are cheaper to
compute over many cell states.

## Quick start

```r
library(SCINET)
library(igraph)

# `A` is an imputed gene-by-cell expression matrix with gene symbols as rownames.
# `G` is a reference interactome as an igraph object with named vertices.
data(PCNet)                                  # 18,820-gene human interactome
G <- graph_from_adjacency_matrix(PCNet, mode = "undirected")

Bcells <- which(cell_labels == "B cell")

nets <- run.SCINET(G, A, samples = Bcells,
                   total_subsamples = 100, cells_per_subsample = 10,
                   thread_no = 8)

# Aggregate the subsamples into a mean network and score topological specificity
mu <- Reduce(`+`, nets) / length(nets)
g  <- graph_from_adjacency_matrix(as.matrix(mu), mode = "undirected", weighted = TRUE)
spec <- topo.spec(g, sample_no = 100)
head(sort(spec, decreasing = TRUE), 20)
```

For finer control, the stages are separately callable:

```r
activity <- compute_gene_activities(A = A, samples = Bcells, thread_no = 8)
paired   <- pair.datasets(G, activity)
nets     <- construct_cell_networks(net = as(igraph::as_adjacency_matrix(paired$net),
                                             "TsparseMatrix"),
                                    gene_activities = paired$activity.scores,
                                    thread_no = 8)
```

Worked vignettes are under `demo/`, covering the basic pipeline, the
ACTIONet-decoupled path, and cell-state rather than cell-type networks.

> The `demo/` folder carries roughly 200 MB of pre-computed `.RDS` inputs so
> that the vignettes run as written. Clone with `--depth 1` if you only want
> the package.

## Function reference

### Core, implemented in C++

| Function | Purpose |
| --- | --- |
| `RIN_transform` | Rank-based inverse normal transform of an expression matrix |
| `compute_gene_activities` | Activity scores for a set of cells, with baseline correction |
| `compute_gene_activities_decoupled` | The same, computed from ACTIONet archetypes and their loadings |
| `compute_gene_activities_full` | Activity scores over every cell |
| `subsample_gene_activities` | Repeated subsampling of cells into activity-score columns |
| `subsample_gene_activities_decoupled` | The archetype-decoupled equivalent |
| `construct_cell_networks` | Rescore the edges of a reference network per subsample |
| `construct_cell_networks_noPrior` | All-to-all networks, with no reference network |
| `construct_cell_networks_summary` | Subsample and aggregate in one call |

### R helpers

| Function | Purpose |
| --- | --- |
| `run.SCINET` | End-to-end driver: subsample, pair, and construct |
| `pair.datasets` | Restrict a network and an activity matrix to their shared genes |
| `topo.spec` | Topological specificity score per gene, against a degree-matched null |
| `compute.genesets.compactness` | How tightly a gene set clusters in a weighted network |
| `prioritize.geneset.LOO` | Leave-one-out prioritisation within a gene set |
| `getLCC` | Largest connected component of an adjacency matrix |
| `read.edgelist` | Read a tab-separated edge list into an igraph object |

## Bundled data

| Dataset | Contents |
| --- | --- |
| `PCNet` | Parsimonious Composite Network v1.3, 18,820 genes, from [NDEx](http://www.ndexbio.org/#/network/4de852d9-9908-11e9-bcaf-0ac135e8bacf) |
| `MIPPIE` | Mouse Integrated Protein-Protein Interaction rEference, 10,886 genes |

## Citation

> Mohammadi, S., Davila-Velderrain, J., & Kellis, M. (2019).
> **Reconstruction of Cell-type-Specific Interactomes at Single-Cell
> Resolution.** *Cell Systems*, 9(6), 559-568.e4.
> https://doi.org/10.1016/j.cels.2019.10.007

`CITATION.cff` is included, so GitHub's "Cite this repository" button will
generate BibTeX and APA entries.

## License

GPL (>= 2); see [LICENSE](LICENSE).

## Contributing

Bug reports and pull requests are welcome via
[GitHub issues](https://github.com/shmohammadi86/SCINET/issues).
