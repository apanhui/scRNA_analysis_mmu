# ref_index_dir ：参考的索引
# *_fq_dir：每个样本的原始数据路径
# cellranger v5.0,0
cellranger count --id dbdb --transcriptome ref_index_dir --disable-ui --fastqs dbdb_fq_dir --sample dbdb --chemistry SC3Pv3 --include-introns --jobmode sge --maxjobs 5 --mempercore 6
cellranger count --id dbdb-KD --transcriptome ref_index_dir --disable-ui --fastqs dbdb-KD_fq_dir --sample dbdb-KD --chemistry SC3Pv3 --include-introns --jobmode sge --maxjobs 5 --mempercore 6
cellranger count --id dbm --transcriptome ref_index_dir --disable-ui --fastqs dbm_fq_dir --sample dbm --chemistry SC3Pv3 --include-introns --jobmode sge --maxjobs 5 --mempercore 6
