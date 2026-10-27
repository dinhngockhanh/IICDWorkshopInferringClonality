# IICD Parameter Estimation Workshop: Inferring Cancer Clonality from DNA

## Installation

1. Install Docker [here](https://docs.docker.com/desktop/) then open it
2. Download this library with

```bash
git clone https://github.com/dinhngockhanh/IICDWorkshopInferringClonality.git
```

3. In Terminal, go to the downloaded folder and run

```bash
docker compose build
```

## Activities during the Workshop

We will first review definitions of tumor characteristics and important metrics when analyzing DNA-sequencing data from [1].
We will then examine each mutation clustering method in detail, using a synthetic tumor sample from the ICGC-TCGA DREAM Challenge [2, 3] as example.
The different options for which algorithm to implement can be listed by running

```bash
docker compose run --rm workshop
```

### PhyloWGS

We can run PhyloWGS [4] with

```bash
docker compose run --rm workshop bash src/run_phylowgs.sh
```

which inputs the mutation table from file `PhyloWGS_ssm.txt` and the copy number profile from file `PhyloWGS_cnv.txt` in folder [data/](https://github.com/dinhngockhanh/IICDWorkshopInferringClonality/tree/main/data).

To save time during the Workshop, file `run_phylowgs.sh` executes PhyloWGS with a smoke test (`MODE = "smoke"`; one MCMC chain with parameters `BURNIN = 5, MCMC = 5, MH = 50`).
Set `MODE = "full"` at the top of [src/run_phylowgs.sh](https://github.com/dinhngockhanh/IICDWorkshopInferringClonality/blob/main/src/run_phylowgs.sh) to run four chains (`BURNIN = 1000, MCMC = 2500, MH = 5000`, seeds 1001–1004). That run may take days. After it finishes, build the one-page view with `python3 src/make_phylowgs_html.py`.
A precomputed page from that four-chain run is [S8-noXY.html](S8-noXY.html) in this directory. Double-click it to open it.

### PyClone-VI

We can run PyClone-VI [5] with

```bash
docker compose run --rm workshop bash src/run_pyclonevi.sh
```

which inputs the mutation table with copy number and sample purity annotation from file `PyCloneVI_mutation_table.tsv` in folder [data/](https://github.com/dinhngockhanh/IICDWorkshopInferringClonality/tree/main/data).

### SciClone

We can run SciClone [6] with

```bash
docker compose run --rm workshop Rscript src/run_sciclone.R
```

which inputs the mutation table from file `sciClone_mutation_table.tsv`, the copy number profile from file `sciClone_cn.tsv` and list of regions with loss of heterogeneity (LOH) to exclude in the analysis from file `sciClone_exclude_loh.bed` in folder [data/](https://github.com/dinhngockhanh/IICDWorkshopInferringClonality/tree/main/data).

### MOBSTER

We can run MOBSTER [7] with

```bash
docker compose run --rm workshop Rscript src/run_mobster.R
```

which inputs the mutation table from file `MOBSTER_mutation_table.csv` in folder [data/](https://github.com/dinhngockhanh/IICDWorkshopInferringClonality/tree/main/data).

### DECODE

We can run DECODE [8] with

```bash
docker compose run --rm workshop Rscript src/run_decode.R
```

which inputs the mutation table from file `DECODE_mutation_table.tsv` in folder [data/](https://github.com/dinhngockhanh/IICDWorkshopInferringClonality/tree/main/data).

## References

1. Tarabichi M, Salcedo A, Deshwar AG, Ni Leathlobhair M, Wintersinger J, Wedge DC, Van Loo P, Morris QD, Boutros PC. A practical guide to cancer subclonal reconstruction from DNA sequencing. Nature methods (2021). https://doi.org/10.1038/s41592-020-01013-2
2. Salcedo A, Tarabichi M, Espiritu SM, Deshwar AG, David M, Wilson NM, Dentro S, Wintersinger JA, Liu LY, Ko M, Sivanandan S. A community effort to create standards for evaluating tumor subclonal reconstruction. Nature biotechnology (2020). https://doi.org/10.1038/s41587-019-0364-z
3. Salcedo A, Tarabichi M, Buchanan A, Espiritu SM, Zhang H, Zhu K, Ou Yang TH, Leshchiner I, Anastassiou D, Guan Y, Jang GH. Crowd-sourced benchmarking of single-sample tumor subclonal reconstruction. Nature biotechnology (2025). https://doi.org/10.1038/s41587-024-02250-y
4. Deshwar AG, Vembu S, Yung CK, Jang GH, Stein L, Morris Q. PhyloWGS: reconstructing subclonal composition and evolution from whole-genome sequencing of tumors. Genome biology (2015). https://doi.org/10.1186/s13059-015-0602-8
5. Gillis S, Roth A. PyClone-VI: scalable inference of clonal population structures using whole genome data. BMC bioinformatics (2020). https://doi.org/10.1186/s12859-020-03919-2
6. Miller CA, White BS, Dees ND, Griffith M, Welch JS, Griffith OL, Vij R, Tomasson MH, Graubert TA, Walter MJ, Ellis MJ. SciClone: inferring clonal architecture and tracking the spatial and temporal patterns of tumor evolution. PLoS computational biology (2014). https://doi.org/10.1371/journal.pcbi.1003665
7. Caravagna G, Heide T, Williams MJ, Zapata L, Nichol D, Chkhaidze K, Cross W, Cresswell GD, Werner B, Acar A, Chesler L. Subclonal reconstruction of tumors by using machine learning and population genetics. Nature genetics (2020). https://doi.org/10.1038/s41588-020-0675-5
8. Chen Y, Jaksik R, Terranova P, El Baghdadi S, Koval A, Kurpas MK, Tavaré S, Kimmel M, Dinh KN. Accurate detection of tumor clonality and ongoing expansion mode from genomic data. bioRxiv (2026). https://doi.org/10.64898/2026.06.15.732415
