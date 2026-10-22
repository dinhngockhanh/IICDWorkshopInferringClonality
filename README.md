#   IICD Parameter Estimation Workshop: Inferring Cancer Clonality from DNA

##  Installation

1.  Install and open [Docker](https://docs.docker.com/desktop/)

2.  Download the library with

```bash
git clone https://github.com/dinhngockhanh/IICDWorkshopInferringClonality.git
```

3.  In Terminal, go to the downloaded folder and run

```bash
docker compose build
```

##  Activities during the presentation

We will review definitions of tumor characteristics in [1].
We will then go through the mechanisms of different mutation clustering methods, including PhyloWGS [2], PyClone-VI [3], SciClone [4], MOBSTER [5] and DECODE [6].
We will compare their performance in ICGC-TCGA DREAM Challenge, a synthetic benchmarking study [7, 8].

##  References
1.  Tarabichi M, Salcedo A, Deshwar AG, Ni Leathlobhair M, Wintersinger J, Wedge DC, Van Loo P, Morris QD, Boutros PC. A practical guide to cancer subclonal reconstruction from DNA sequencing. Nature methods (2021). https://doi.org/10.1038/s41592-020-01013-2
2.  Deshwar AG, Vembu S, Yung CK, Jang GH, Stein L, Morris Q. PhyloWGS: reconstructing subclonal composition and evolution from whole-genome sequencing of tumors. Genome biology (2015). https://doi.org/10.1186/s13059-015-0602-8
3.  Gillis S, Roth A. PyClone-VI: scalable inference of clonal population structures using whole genome data. BMC bioinformatics (2020). https://doi.org/10.1186/s12859-020-03919-2
4.  Miller CA, White BS, Dees ND, Griffith M, Welch JS, Griffith OL, Vij R, Tomasson MH, Graubert TA, Walter MJ, Ellis MJ. SciClone: inferring clonal architecture and tracking the spatial and temporal patterns of tumor evolution. PLoS computational biology (2014). https://doi.org/10.1371/journal.pcbi.1003665
5.  Caravagna G, Heide T, Williams MJ, Zapata L, Nichol D, Chkhaidze K, Cross W, Cresswell GD, Werner B, Acar A, Chesler L. Subclonal reconstruction of tumors by using machine learning and population genetics. Nature genetics (2020). https://doi.org/10.1038/s41588-020-0675-5
6.  Chen Y, Jaksik R, Terranova P, El Baghdadi S, Koval A, Kurpas MK, Tavaré S, Kimmel M, Dinh KN. Accurate detection of tumor clonality and ongoing expansion mode from genomic data. bioRxiv (2026). https://doi.org/10.64898/2026.06.15.732415
7.  Salcedo A, Tarabichi M, Espiritu SM, Deshwar AG, David M, Wilson NM, Dentro S, Wintersinger JA, Liu LY, Ko M, Sivanandan S. A community effort to create standards for evaluating tumor subclonal reconstruction. Nature biotechnology (2020). https://doi.org/10.1038/s41587-019-0364-z
8.  Salcedo A, Tarabichi M, Buchanan A, Espiritu SM, Zhang H, Zhu K, Ou Yang TH, Leshchiner I, Anastassiou D, Guan Y, Jang GH. Crowd-sourced benchmarking of single-sample tumor subclonal reconstruction. Nature biotechnology (2025). https://doi.org/10.1038/s41587-024-02250-y



# S8 亚克隆重建

这个仓库在 DREAM 单样本 `S8-noXY` 上运行 DECODE、MOBSTER、sciClone、PyClone-VI 和 PhyloWGS，并把每个工具的结果收成同一张突变分配表。

算法本体使用已发表的软件。仓库里的脚本只负责调用一次、再提取结果。输入已经放在 `data/`。每个工具吃哪些突变，写在 [docs/data.md](docs/data.md) 和 [docs/algorithms.md](docs/algorithms.md)。

## 运行

需要 Docker。宿主机不必安装 R、conda 或这些 R 包。五个算法分开运行，没有总入口。

```bash
docker compose build
```

第一次构建会安装 R 包、PyClone-VI，并编译 Python 2.7 和 PhyloWGS 的 `mh.o`，时间很长。不带命令的 `docker compose run --rm workshop` 只打印下面这五条命令。

```bash
docker compose run --rm workshop Rscript src/run_decode.R
docker compose run --rm workshop Rscript src/run_mobster.R
docker compose run --rm workshop Rscript src/run_sciclone.R
docker compose run --rm workshop bash src/run_pyclonevi.sh
docker compose run --rm workshop bash src/run_phylowgs.sh
```

每条命令在拟合结果写盘之后，会自己写出 `results/<tool>/clones.tsv`。列是 `purity`、`n_clones`、`clone`、`ccf`。DECODE 和 PhyloWGS 给出自己估计的 purity。PyClone-VI 把 purity 当作输入，所以写 `NA`。MOBSTER 和 sciClone 的 `purity` 是各 clone CCF 里的最大值。

同一轮还会写 PDF。MOBSTER 和 sciClone 用软件自带的拟合图。DECODE 写模型选择图，以及它选中的那一个模型的 SFS。PyClone-VI 不返回拟合密度，图是堆积的观测 VAF，Cluster 0 在最下方。PhyloWGS 没有 SFS 拟合图，脚本画的是和 `clones.tsv` 同一棵树的细胞比例；官方 `witness` 网页不在这条命令里。

Apptainer 使用同一个镜像，把上面的命令换成容器里的那一段：

```bash
apptainer build workshop.sif docker-daemon://clonality-workshop:s8
ulimit -s unlimited
apptainer exec --bind "$PWD":/work --pwd /work workshop.sif Rscript src/run_mobster.R
```

DECODE 使用包的默认粒子数，在约 8000 个二倍体突变上可能要数小时。PyClone-VI 通常是分钟级，第一次还会编译 numba。MOBSTER 和 sciClone 建议按 16 GB 内存准备。PhyloWGS 默认只跑短链冒烟。

要把 PhyloWGS 改成正式链，把 `src/run_phylowgs.sh` 开头的 `MODE` 改成 `"full"`，再运行这个脚本。正式链是 burn-in 1000、2500 次 MCMC、1000 次 MH、一条链。这一档以小时到天计。冒烟链和正式链不要放在一起比较。

## 目录

```text
data/DECODE_mutation_table.tsv
data/DECODE_mutation_table_with_multiplicity.csv   全部突变，只用于 DECODE 的 soft assignment
data/MOBSTER_mutation_table.csv
data/sciClone_mutation_table.tsv
data/sciClone_cn.tsv
data/sciClone_exclude_loh.bed
data/PyCloneVI_mutation_table.tsv
data/PhyloWGS_ssm.txt
data/PhyloWGS_cnv.txt
data/PhyloWGS_ssm_map.tsv
data/mutations.tsv      全部突变清单，提取 clone 表时不使用
data/sample.tsv         purity 和 ploidy；MOBSTER、sciClone 用它把 VAF 换成 CCF
src/run_*.R             DECODE、MOBSTER、sciClone；参数和 clone 提取都在脚本里
src/run_*.sh            PyClone-VI、PhyloWGS；参数和 clone 提取都在脚本里
results/                运行后生成，不提交
```

每个工具的摘要在 `results/<tool>/clones.tsv`。

## 环境为什么不是两个锁文件

`uv.lock` 锁定 Python 3.12 和从 GitHub 安装的 PyClone-VI 0.2.0。`renv.lock` 在 `rocker/r-ver:4.5.1` 里生成，锁定 DECODE、MOBSTER、sciClone 及其依赖。学生不在自己的系统上执行这两个锁。

这两个文件安装不了 PhyloWGS。它需要 Python 2.7.18、旧版 numpy/scipy/ete2，以及用 GSL 编译的 `mh.o`。这些版本写在 `Dockerfile` 的 ARG 里。镜像同时提供 R、编译器和系统库。学生每次只选一个算法运行。

刷新 R 锁文件时，用备课脚本，不要在 macOS 上生成：

```bash
docker build -f docker/Dockerfile.renv-bootstrap -t clonality-renv-bootstrap .
docker run --rm -v "$PWD":/out clonality-renv-bootstrap
```

## 提取

每个运行脚本在保存拟合结果后，用本文件里的函数写 `results/<tool>/clones.tsv`，并写出上面的 PDF。这张表不计算分数。DECODE 另外用 `DECODE_mutation_assignment` 给每个突变写 soft assignment。规则在 [docs/algorithms.md](docs/algorithms.md)。参数写在对应脚本的开头。
