"""
从https://www.ncbi.nlm.nih.gov/Traces/wgs/AEYC01?display=contigs&page=1下载相应的Accession后
根据Gene列的标注提取出目标基因序列
"""

# pip install biopython

from Bio import SeqIO

# 读取FASTA文件
fasta_file = "../data/AEYC01000002.1_sequence.fasta"
for record in SeqIO.parse(fasta_file, "fasta"):
    print(f"ID: {record.id}")
    print(f"Description: {record.description}")

# 基因序列：record.seq[5978..6592]
