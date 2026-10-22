# S8-noXY 输入是怎么来的

学生不生成这些文件。`data/` 里已经是各算法能直接读的表：

- `DECODE_mutation_table.tsv`
- `DECODE_mutation_table_with_multiplicity.csv`：全部突变，只给 DECODE 的 soft assignment 用。拟合不读它
- `MOBSTER_mutation_table.csv`
- `sciClone_mutation_table.tsv`、`sciClone_cn.tsv`、`sciClone_exclude_loh.bed`
- `PyCloneVI_mutation_table.tsv`
- `PhyloWGS_ssm.txt`、`PhyloWGS_cnv.txt`、`PhyloWGS_ssm_map.tsv`
- `mutations.tsv`：全部突变清单，不是某个算法的输入，提取 clone 表时也不读它
- `sample.tsv`：purity 和 ploidy。MOBSTER 和 sciClone 只用它把簇的 VAF 换成 CCF

备课材料里的三份原始文件是突变表、Battenberg 和 cellularity。突变表有 10751 行，列是 `Chromosome,Position,major_cn,minor_cn,Ref_count,Alt_count,multiplicity`，没有深度为 0 的行。Battenberg 有 27 个区段，只用 solution A，也就是 `nMaj1_A`、`nMin1_A`、`frac1_A`。`cellularity = 0.683766856994925` 当作 purity，`ploidy = 1.8958699316085`，写在 `data/sample.tsv`。

区段表和计数核对留在备课目录，不进入学生的 `data/`。

`mutations.tsv` 的 `mutation_id` 是 `染色体:位置`。染色体去掉可选的 `chr` 前缀，位置用 `int(float(...))`。拷贝数计数：

- (1, 1)：7898
- (1, 0)：1930
- (2, 1)：918
- (0, 0)：5

multiplicity 为 1 的有 10501 个，为 2 的有 250 个。二倍体行的 multiplicity 全部是 1。

落在 Battenberg 区段里的突变，其 `major_cn`/`minor_cn` 与 solution A 一致。solution B–F 不使用。

## 区段外面的突变

137 个突变落在所有 Battenberg 区段之外。其中 111 个是二倍体。PhyloWGS 对没有列入 CNV 的 SSM 本来就按二倍体处理，所以这 111 个不另建区段。

其余非二倍体突变里，25 个被收成 6 个补充区段（备课区段表的 `source = gap_fill`，`frac = 1`）。只有同一染色体、同一拷贝数、并且中间没有别的突变或 Battenberg 区段时，才会并成一段。这样补充区段不会把旁边的二倍体突变盖进去。

5 个 (0, 0) 突变不进入 PyClone-VI，也不进入 PhyloWGS 的 CNV。它们是：

- `4:132604652`
- `4:132691244`
- `22:16208440`
- `22:16620565`
- `22:17049725`

对应的 solution A 区段本身也是 (0, 0)，同样不建成 CNV。这些突变仍然留在 SSM 表里，PhyloWGS 会把它们当成二倍体。备课核对文件记录了这些 id。

## 三个亚克隆区段

`frac1_A` 不是 1 的区段有三个：

- 4 号 `132581378–132890699`，拷贝数 (0, 0)，比例约 0.601。不进 CNV。
- 13 号 `19589846–114999838`，拷贝数 (1, 0)，比例约 0.668。CNV 的 `cell_prev` 用这个比例。
- 22 号 `16600330–17208423`，拷贝数 (0, 0)，比例约 0.554。不进 CNV。

旧的 TCGA 脚本把 `cell_prev` 写成 1.0，是因为 ASCAT 区段没有细胞比例。这里 Battenberg 有比例，所以用 `frac1_A`。

## 各算法实际用多少行

- DECODE 拟合：二倍体 7898 行。不看 multiplicity，也不看 purity。
- DECODE assignment：全部 10751 行，列是 `mutation_id,Chromosome,Position,major_cn,minor_cn,Ref_count,Alt_count,multiplicity`。`normalized_VAF` 不在这张表里，运行时用拟合得到的 purity 计算。
- MOBSTER 和 sciClone 的输入文件：二倍体、multiplicity 为 1、且 VAF 在 (0, 1) 的 7898 行。sciClone 运行时还会用 `minimumDepth = 10` 丢掉其中深度低于 10 的 6 行。
- PyClone-VI：multiplicity 为 1 且总拷贝数不是 0，10496 行。`normal_cn = 2`，`tumour_content` 固定为 cellularity。
- PhyloWGS：10751 个 SSM，14 个 CNV。CNV 只保留总拷贝数不是 2、不是 0、且至少含一个 SSM 的区段。`a` 是参考读数，`d` 是深度。区段没有自己的覆盖度，`a`/`d` 用段内 SSM 的参考读数之和与深度之和。`mu_r = 0.999`，`mu_v = 0.5`。
