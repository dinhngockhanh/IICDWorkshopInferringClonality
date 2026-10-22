# 五个算法各看什么

上游程序保持原样。学生要读的是 `src/00_prepare.py` 里的过滤，以及 `src/run_*.R`、`src/run_*.sh` 里的那一次调用。参数在 `config/params.toml`。

正式参数与集群上的单样本耗时测试一致。PhyloWGS 默认不是那次耗时测试的长度。

## DECODE

输入是二倍体突变的 `Ref_count` 和 `Alt_count`。调用是：

```r
DECODE(
    sample_id = sample_id,
    mutation_table = mutations[, c("Ref_count", "Alt_count")],
    neutral_tail = FALSE,
    compute_parallel = FALSE
)
```

粒子数用包的默认值 `rep(1000, 10)`，脚本里不覆盖。约 8000 个二倍体突变上，这一步可能要数小时。

拟合结果给的是每个簇的 VAF 后验，不给每条突变的标签。`src/extract_r.R` 在拟合完成之后做这件事：用 `2 * mean(Cluster_1_freq)` 作为 DECODE 自己的 purity，用表里已经有的 multiplicity 把 VAF 调到二倍体尺度，再按 SFS 分箱比例做硬分配。这一步不调用 GRITIC。

## MOBSTER

输入是二倍体、multiplicity 为 1、VAF 在 (0, 1) 的突变。调用 `mobster_fit(..., K = 1:3, tail = c(TRUE, FALSE), parallel = FALSE)`。`mobster_fit` 会用它自己的 `seed` 参数重置随机数，所以配置里的种子是传进这个参数的，而不是只靠外层的 `set.seed`。

安装来源是 `caravagnalab/mobster` 的 commit `85c898f087b46e15f79144bf57bdc019678e2481`。旧文档里的 `caravagn/mobster` 已经转到这个仓库。

统一表里的 `ccf` 用该簇的平均 VAF：二倍体、multiplicity 为 1 时，`ccf = 2 * mean(VAF) / purity`。purity 用 cellularity。Tail 行的 `cluster` 保留 MOBSTER 的标签，`ccf` 留空，`note` 为 `tail`。

## sciClone

输入文件与 MOBSTER 是同一批突变。VAF 写成 0–100 的百分数，染色体不带 `chr`。同时传入 `cn.tsv` 和 LOH bed。调用：

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

`cluster.means` 是样本数乘簇数的矩阵，数值是 0–1 的 VAF。统一表里的 `ccf = 2 * cluster_mean / purity`。深度不足或被拷贝数过滤掉的突变 `included = false`，`note` 为 `filtered by sciClone`。

## PyClone-VI

`pyclone-vi fit -c 10 -d beta-binomial -r 10`，然后 `write-results-file`。`cluster` 用结果里的 `cluster_id`，`ccf` 用 `cellular_prevalence`。版本是 Git 标签 0.2.0，commit `07306831a9a48275ccfe43fe42bbe18d7370bc72`。这个包不在 PyPI 上。

## PhyloWGS

默认 `phylowgs.mode = "smoke"`：`-B 5 -s 5 -i 50`，一条链，种子 1001。冒烟只检查安装和输入，不能拿来和正式结果比较。

正式链把 `config/params.toml` 改成 `mode = "full"`，对应 `-B 1000 -s 2500 -i 1000`，仍然是一条链。这一档在集群上大约要十几个小时。集群耗时测试的 4 条链、7 档 MH 迭代不在默认命令里。

`mutass.json` 是一个 zip，里面每个树一个 JSON。`summ.json` 是 gzip。提取时选 `tree_densities` 最大的那棵树。`cluster` 是该树里的群体编号，`ccf` 是该群体的 `cellular_prevalence`。`note` 写成 `mcmc=smoke;tree=<id>` 或 `mcmc=full;tree=<id>`。

Python 2.7.18、numpy 1.16.6、scipy 1.2.3、ete2 2.3.10，以及 commit `2205be16e9ad68426e837229b593f47e36b3547e`，都写在 Dockerfile 的 ARG 里。`mh.o` 在镜像构建时用 GSL 编译。
