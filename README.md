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

算法本体使用已发表的软件。仓库里的脚本只负责准备输入、调用一次、再提取结果。每个工具吃哪些突变，写在 [docs/data.md](docs/data.md) 和 [docs/algorithms.md](docs/algorithms.md)。

## 运行

需要 Docker。宿主机不必安装 R、conda 或这些 R 包。

```bash
docker compose build
docker compose run --rm workshop
```

第一次构建会安装 R 包、PyClone-VI，并编译 Python 2.7 和 PhyloWGS 的 `mh.o`，时间很长。构建完成后，重复运行就是上面的第二条命令。

Apptainer 使用同一个镜像：

```bash
apptainer build workshop.sif docker-daemon://clonality-workshop:s8
ulimit -s unlimited
apptainer exec --bind "$PWD":/work --pwd /work workshop.sif bash src/run_all.sh
```

`src/run_all.sh` 会先重建 `data/prepared/`，然后按这个顺序运行：PyClone-VI、MOBSTER、sciClone、PhyloWGS 冒烟、DECODE。每个工具结束后立刻提取它自己的表。前面的工具失败不会阻止后面的工具；失败名单写在 `results/failures.txt`。

DECODE 使用包的默认粒子数，在约 8000 个二倍体突变上可能要数小时。PyClone-VI 通常是分钟级，第一次还会编译 numba。MOBSTER 和 sciClone 建议按 16 GB 内存准备。PhyloWGS 默认只跑短链冒烟。

要把 PhyloWGS 改成正式链，编辑 `config/params.toml`，把 `phylowgs.mode` 设为 `"full"`，再运行：

```bash
docker compose run --rm workshop bash src/run_phylowgs.sh
docker compose run --rm workshop python3 src/extract_results.py --tool phylowgs
```

正式链是 burn-in 1000、2500 次 MCMC、1000 次 MH、一条链。这一档以小时到天计。冒烟结果的 `note` 列会标 `mcmc=smoke`，不要拿去和正式链比较。

## 目录

```text
data/raw/          三份原始文件
data/prepared/     00_prepare.py 生成的规范表和五个输入
src/00_prepare.py  过滤和格式转换
src/run_*.R        DECODE、MOBSTER、sciClone
src/run_*.sh       PyClone-VI、PhyloWGS
src/extract_results.py
src/extract_r.R    读取 RDS，供上面的提取脚本调用
src/score.py       评分占位
config/params.toml 种子和各工具参数
results/           运行后生成，不提交
```

统一结果表在 `results/assignments/`。列是 `mutation_id,chrom,pos,tool,included,cluster,ccf,note`。某个工具没有使用的突变，`included` 为 `false`。

## 环境为什么不是两个锁文件

`uv.lock` 锁定 Python 3.12 和从 GitHub 安装的 PyClone-VI 0.2.0。`renv.lock` 在 `rocker/r-ver:4.5.1` 里生成，锁定 DECODE、MOBSTER、sciClone 及其依赖。学生不在自己的系统上执行这两个锁。

这两个文件安装不了 PhyloWGS。它需要 Python 2.7.18、旧版 numpy/scipy/ete2，以及用 GSL 编译的 `mh.o`。这些版本写在 `Dockerfile` 的 ARG 里。镜像同时提供 R、编译器和系统库，所以一条 `docker compose` 命令就够了。

刷新 R 锁文件时，用备课脚本，不要在 macOS 上生成：

```bash
docker build -f docker/Dockerfile.renv-bootstrap -t clonality-renv-bootstrap .
docker run --rm -v "$PWD":/out clonality-renv-bootstrap
```

## 评分

`src/score.py` 只列出已经生成的分配表，然后说明评分规则还没有给出。它不计算 SMC-Het，也不会因为缺少分数而让整次运行失败。
