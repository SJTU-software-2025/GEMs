# 测试用脚本，可忽略
# 验证Accession（如AEYC01000002.1）是否完整存在于从https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000194155.1/下载得到的GCA_000194155.1_ASM19415v1_genomic.fna中
# 文件路径等细节请根据实际情况修改
hisat2-build ./GCA_000194155.1_ASM19415v1_genomic.fna  ./GCA_ASM19415v1_genomic_index
hisat2 -x ./GCA_ASM19415v1_genomic_index -f ./AEYC01000002.1_sequence.fasta -S ./GCA_000194155.1_ASM19415v1_genomic.sam
less ./GCA_000194155.1_ASM19415v1_genomic.sam

# samtools faidx GCA_000194155.1_ASM19415v1_genomic.fna GL877878.1:460814-485300