# 五个算法各看什么

上游程序保持原样。学生要读的是 `src/run_*.R`、`src/run_*.sh` 里的那一次调用。输入已经在 `data/` 里，过滤规则见 `docs/data.md`。每个脚本开头写着自己的参数。五个脚本互不调用。每个脚本在拟合结果写盘之后，用本文件里的函数提取 `results/<tool>/clones.tsv`，并写出该工具的 PDF。

正式参数与集群上的单样本耗时测试一致。PhyloWGS 默认不是那次耗时测试的长度。

## DECODE

输入是 `data/DECODE_mutation_table.tsv` 里二倍体突变的 `Ref_count` 和 `Alt_count`。调用是：

```r
DECODE(
    sample_id = sample_id,
    mutation_table = mutations[, c("Ref_count", "Alt_count")],
    neutral_tail = FALSE,
    compute_parallel = FALSE
)
```

粒子数用包的默认值 `rep(1000, 10)`，脚本里不覆盖。约 8000 个二倍体突变上，这一步可能要数小时。

拟合结果给的是每个簇的 VAF 后验。提取不给每条突变贴标签，也不调用 GRITIC。`purity = 2 * mean(Cluster_1_freq)`。clone 数是 `best_N_clusters`。第 k 个 clone 的 `ccf = mean(Cluster_k_freq) / mean(Cluster_1_freq)`。

图用包里的函数，画在 RDS 写完之后。`DECODE_plot_model_selection` 写出 `S8-noXY_DECODE_model_selection.pdf`：读数阈值、各簇数的选择准则，以及每个 tail 和簇数配置的 SFS。`DECODE_plot_SFS` 写出 `S8-noXY_DECODE_SFS.pdf`，只画 DECODE 自己选中的那个模型。`filename` 不带 `.pdf`，函数自己补上。宽高用包的默认值。

`clones.tsv` 写完之后，`DECODE_mutation_assignment` 给全部突变做 soft assignment。输入是 `data/DECODE_mutation_table_with_multiplicity.csv`，拟合不读这张表。观察 VAF 先换成二倍体尺度：`normalized_VAF = VAF * (c * purity + 2 * (1 - purity)) / (2 * m)`，其中 `c` 是 major 加 minor，`m` 是 multiplicity，purity 与 `clones.tsv` 相同。函数用默认参数，也就是选中的模型和 `inference_A`。结果是 `results/decode/S8-noXY_mutation_assignment.csv`，多出来的列是 `prob_cluster_*`；若选中的模型带 tail，还有 `prob_tail`。这是概率，不是硬标签。

## MOBSTER

输入是 `data/MOBSTER_mutation_table.csv`：二倍体、multiplicity 为 1、VAF 在 (0, 1) 的突变。调用 `mobster_fit(..., K = 1:3, tail = c(TRUE, FALSE), parallel = FALSE)`。`mobster_fit` 会用它自己的 `seed` 参数重置随机数，所以配置里的种子是传进这个参数的，而不是只靠外层的 `set.seed`。

安装来源是 `caravagnalab/mobster` 的 commit `85c898f087b46e15f79144bf57bdc019678e2481`。旧文档里的 `caravagn/mobster` 已经转到这个仓库。

拟合不读 purity。Tail 不是 clone。其余每个簇在二倍体、multiplicity 为 1 时，`ccf = 2 * mean(VAF) / purity`，这里换算用的 purity 来自 `data/sample.tsv` 的 cellularity。写进摘要的 `purity` 是这些 CCF 的最大值。`src/run_mobster.R` 在 RDS 和 runtime 写完后提取这张表。

## sciClone

输入是 `data/sciClone_mutation_table.tsv`，与 MOBSTER 是同一批突变。VAF 写成 0–100 的百分数，染色体不带 `chr`。同时传入 `data/sciClone_cn.tsv` 和 `data/sciClone_exclude_loh.bed`。调用：

```r
sciClone(
    vafs = vafs,
    copyNumberCalls = copy_number,
    sampleNames = sample_id,
    regionsToExclude = loh,
    minimumDepth = 10,
    doClustering = TRUE,
    doClusteringAlongMargins = FALSE
)
```

`cluster.means` 是样本数乘簇数的矩阵，数值是 0–1 的 VAF。clone 数是这个矩阵的簇数。`ccf = 2 * cluster_mean / purity`，换算用 `data/sample.tsv` 的 cellularity。摘要里的 `purity` 是这些 CCF 的最大值。深度过滤发生在 sciClone 内部，不另记到摘要里。`src/run_sciclone.R` 在 RDS 和 runtime 写完后提取这张表。

## PyClone-VI

输入是制表符分隔的 `data/PyCloneVI_mutation_table.tsv`，其中 `tumour_content` 就是 purity。`pyclone-vi fit -c 10 -d beta-binomial -r 10`，然后 `write-results-file`。因为 purity 是输入，摘要里的 `purity` 写 `NA`。clone 数是结果里不同的 `cluster_id` 个数，每个 clone 的 `ccf` 是该簇的 `cellular_prevalence`。版本是 Git 标签 0.2.0，commit `07306831a9a48275ccfe43fe42bbe18d7370bc72`。这个包不在 PyPI 上。

PyClone-VI 不返回混合密度，所以没有参数化的 SFS 曲线。`S8-noXY_sfs.pdf` 是观测 VAF 的堆积直方图，`VAF = alt_counts / (ref_counts + alt_counts)`。同一根柱子里 Cluster 0 在最下方，簇号更大的在上面。图在结果 TSV 写完之后用镜像里的 `ggplot2` 画。

## PhyloWGS

输入是 `data/PhyloWGS_ssm.txt` 和 `data/PhyloWGS_cnv.txt`。样本 purity 不是输入。`src/run_phylowgs.sh` 开头的 `MODE` 默认是 `"smoke"`：`-B 5 -s 5 -i 50`，一条链，种子 1001。冒烟只检查安装和输入，不能拿来和正式结果比较。

正式链把同一个脚本里的 `MODE` 改成 `"full"`，对应 `-B 1000 -s 2500 -i 1000`，仍然是一条链。这一档在集群上大约要十几个小时。集群耗时测试的 4 条链、7 档 MH 迭代不在默认命令里。

`summ.json` 是 gzip。提取时选 `tree_densities` 最大的那棵树。编号 0 是不含突变的正常根，不算 clone。其余每个群体是一个 clone，`ccf` 是它的 `cellular_prevalence`。`purity` 是正常根直接下方那个癌克隆的 `cellular_prevalence`；若这棵树有多个直接子克隆，则把它们的 prevalence 相加。

这个 commit 没有可直接写出的 SFS 拟合 PDF。官方可视化是 `witness/` 网页，要改文件名、解压并起一个 HTTP 服务，不放进学生命令。脚本改画同一棵树：`results/phylowgs/<mode>/S8-noXY_tree.pdf`。节点是群体，边来自 `structure`，标签是 `cellular_prevalence`，群体 0 标成 normal。图题带 `MODE` 和树编号。冒烟和正式链各写在自己的目录里。

Python 2.7.18、numpy 1.16.6、scipy 1.2.3、ete2 2.3.10，以及 commit `2205be16e9ad68426e837229b593f47e36b3547e`，都写在 Dockerfile 的 ARG 里。`mh.o` 在镜像构建时用 GSL 编译。
