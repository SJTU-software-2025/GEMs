# Saccharopolyspora spinosa 基因组尺度代谢模型（Genome scale metabolic model, GSMM）

## Overview

工业微生物是指能够通过工业规模培养获得特定产物或达到特定应用的微生物，其应用涉及工业、农业和医药等诸多领域。
近年来，基因组编辑、代谢工程及机器学习等人工智能技术在微生物菌株改造和选育中应用愈广，
使得重建代谢网络，整合来自其他天然动植物的异源性反应途径，甚至人为从头设计成为可能，显著提高了工业微生物的生产性能。

**基因组尺度代谢模型（GSMM）——细胞的“数字实验室”**

如果把细胞比作一个微型工厂，**基因组尺度代谢模型（GSMM）** 就是科学家为这个工厂建立的“数字孪生系统”。它能通过计算机模拟，预测细胞在不同环境下如何生长、消耗营养、生产能量或有用物质，就像在虚拟世界里做实验一样。

**如何构建这个模型？**

科学家首先从细胞的基因信息入手。基因组中的基因相当于“工厂的蓝图”，指导合成蛋白质（工人）和酶（工具）。通过分析这些基因，可以列出细胞能进行的所有化学反应（生产线）。例如：分解葡萄糖获取能量、合成氨基酸构建细胞等。最终，这些反应会连接成一张巨大的“代谢网络地图”，涵盖成千上万的物质转换步骤。

**GSMM有什么用？**  
1. **生物制造**：比如设计微生物高效生产胰岛素、生物燃料或可降解塑料。通过模拟不同基因改造方案，能快速筛选出最佳“工厂改造计划”。  
2. **医学研究**：分析癌细胞为何疯狂生长？模型可以对比健康细胞和癌细胞的代谢差异，帮助找到精准治疗的靶点。  
3. **环境治理**：预测微生物分解污染物（如石油泄漏）的最佳条件，指导环境修复。  

**优势与挑战**  
- ✅ **省时省钱**：用计算机模拟替代大量实验室试错，原本需数月的实验可能缩短到几小时。  
- ❗ **依赖数据质量**：若基因注释错误（比如误判某个酶的功能），模型预测会“失准”。  
- 🌍 **复杂环境模拟**：实验室条件可控，但真实环境（如人体肠道）的复杂因素仍需进一步研究。

未来，随着人工智能和基因测序技术的发展，GSMM将更精准地指导合成生物学、个性化医疗等领域。这个“细胞模拟器”正在打开生命科学的新可能。

## Project Introduction

本项目构建了 Saccharopolyspora spinosa 的基因组尺度代谢模型，在其代谢空间中进行基因靶点搜索，以多杀菌素关键合成途径为基础分析基因表达情况，研究多杀菌素的高效生物合成，为细胞工厂设计提供理性指导。

---
## Dependencies
- conda 23.10.0 or higher
- Python 3.10
- COBRApy
- tqdm
- zeep
- biopython

---

## Installation
```bash
# 创建新的环境
conda create -n gems python=3.10
conda activate gems

# 安装依赖
pip install cobra tqdm zeep biopython
```

---

## Usage
```bash
# 处理代谢物
python src/metabolite.py

# 添加反应
python src/reaction.py

# 添加边界反应
python src/boundary.py

# 已经存有建立的模型，可直接运行最后的solve
python src/solve.py
```

---

### 边界设定参考表
| 反应类型  | 默认 lower_bound | 默认 upper_bound | 允许负通量 | 允许正通量 | 典型应用场景               |
|-----------|------------------|------------------|------------|------------|----------------------------|
| Exchange  | -1000           | 1000            | ✅          | ✅          | 培养基物质交换             |
| Demand    | 0                | 1000            | ❌          | ✅          | 生物量合成需求             |
| Sink      | -10             | 10              | ✅          | ✅          | 模型拓扑闭合/瞬态平衡      |

---

## Source
https://github.com/SJTU-software-2025/GEMs

---

## Reference

### Documentation
- [Documentation for COBRApy](https://cobrapy.readthedocs.io/en/latest/)
- [KEGG API Manual](https://www.kegg.jp/kegg/rest/keggapi.html)
- [BRENDA SOAP access](https://www.brenda-enzymes.org/soap.php)
- [Biopython API documentation](https://www.osgeo.cn/biopython/)

### Paper
- Wang, X., Zhang, C., Wang, M. et al. Genome-scale metabolic network reconstruction of Saccharopolyspora spinosa for Spinosad Production improvement. Microb Cell Fact 13, 41 (2014). https://doi.org/10.1186/1475-2859-13-41
- Orth, J., Thiele, I. & Palsson, B. What is flux balance analysis?. Nat Biotechnol 28, 245–248 (2010). https://doi.org/10.1038/nbt.1614
- Nogales, J. (2014). A Practical Protocol for Genome-Scale Metabolic Reconstructions. In: McGenity, T., Timmis, K., Nogales , B. (eds) Hydrocarbon and Lipid Microbiology Protocols. Springer Protocols Handbooks. Springer, Berlin, Heidelberg. https://doi.org/10.1007/8623_2014_12
- Ebrahim, A., Lerman, J.A., Palsson, B.O. et al. COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Syst Biol 7, 74 (2013). https://doi.org/10.1186/1752-0509-7-74