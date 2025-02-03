import os
import pickle
from cobra import Metabolite
import pandas as pd
from KEGG import get_chemical_formula
import tqdm

# 临时文件路径
TEMP_FILE = "data/metabolites_temp.pkl"
FINAL_FILE = "data/metabolites.pkl"

# 读取已有进度
metabolites = {}
if os.path.exists(TEMP_FILE):
    with open(TEMP_FILE, "rb") as f:
        metabolites = pickle.load(f)
    print(f"从临时文件恢复了 {len(metabolites)} 个代谢物")

# 读取并处理数据
metabolites_df = pd.read_csv("data/metabolites.csv")
# 过滤掉已处理的行
remaining_df = metabolites_df[~metabolites_df["ID"].isin(metabolites.keys())]

try:
    for _, row in tqdm.tqdm(remaining_df.iterrows(), 
                           total=len(remaining_df),
                           desc="Processing metabolites"):
        metabolite_id = row["ID"]
        metabolite_name = row["Name"]
        compartment = row["Compartment"]
        
        metabolite = Metabolite(
            id=metabolite_id,
            name=metabolite_name,
            compartment=compartment,
            formula=get_chemical_formula(metabolite_id.split("_")[0])
        )
        metabolites[metabolite_id] = metabolite
        
        # 每处理10个保存一次
        if len(metabolites) % 10 == 0:
            with open(TEMP_FILE, "wb") as f:
                pickle.dump(metabolites, f)

    # 全部处理完成，保存最终文件
    with open(FINAL_FILE, "wb") as f:
        pickle.dump(metabolites, f)
    # 删除临时文件
    if os.path.exists(TEMP_FILE):
        os.remove(TEMP_FILE)
        
except Exception as e:
    # 发生错误时保存当前进度
    with open(TEMP_FILE, "wb") as f:
        pickle.dump(metabolites, f)
    print(f"发生错误: {str(e)}")
    print(f"进度已保存到 {TEMP_FILE}")