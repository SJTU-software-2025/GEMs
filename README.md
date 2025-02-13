# Saccharopolyspora spinosa 基因组尺度代谢模型 (GEM)

## 项目简介
本项目构建了 Saccharopolyspora spinosa 的基因组尺度代谢模型，用于研究多杀菌素的生物合成。

## 环境要求
- Python 3.10
- COBRApy
- pandas
- tqdm

## 安装
```bash
# 创建新的环境
conda create -n gems python=3.10
conda activate gems

# 安装依赖
pip install cobra pandas tqdm
```

## 使用
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
