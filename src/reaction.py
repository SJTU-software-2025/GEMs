from cobra import Metabolite, Reaction, Model
from cobra.io import load_json_model, save_json_model
import pickle
import pandas as pd
import re
import tqdm
import os


def clean_model_data(model):
    """清理模型中的NaN值"""
    for reaction in model.reactions:
        if pd.isna(reaction.name):
            reaction.name = ""
        if pd.isna(reaction.gene_reaction_rule):
            reaction.gene_reaction_rule = ""
        if pd.isna(reaction.subsystem):
            reaction.subsystem = ""
        
    for metabolite in model.metabolites:
        if pd.isna(metabolite.name):
            metabolite.name = ""
        if pd.isna(metabolite.formula):
            metabolite.formula = ""

# 临时文件路径
TEMP_FILE = "model/model_temp.pkl"

# 加载或创建新模型
if os.path.exists(TEMP_FILE):
    with open(TEMP_FILE, "rb") as f:
        model = pickle.load(f)
    print(f"从临时文件恢复了模型，包含 {len(model.reactions)} 个反应")
else:
    # 加载 metabolites.pkl
    with open("data/metabolites.pkl", "rb") as f:
        metabolites = pickle.load(f)
    print(f"加载了 {len(metabolites)} 个代谢物")

    # 创建模型
    model = Model("Saccharopolyspora_spinosa")
    
    # 添加代谢物
    print("添加代谢物到模型...")
    for _, metabolite in tqdm.tqdm(metabolites.items()):
        model.add_metabolites(metabolite)

# 加载 reactions.csv
reactions_df = pd.read_csv("data/reactions.csv")
# 过滤已添加的反应
existing_reactions = set(r.id for r in model.reactions)
remaining_df = reactions_df[~reactions_df["Reaction ID"].isin(existing_reactions)]

try:
    # 添加反应
    for _, row in tqdm.tqdm(remaining_df.iterrows(), 
                           total=len(remaining_df),
                           desc="添加反应"):
        try:
            reaction_id = row["Reaction ID"]
            # id 以 "_b" 结尾的反应表示是交换反应 ##### 待修改 #####
            reaction_name = row["Reaction (Long)"]
            reaction_subsystem = row["Subsystem"]
            # 判断反应是否可逆
            if row["Reversible"]:
                lower_bound = -1000.0
                upper_bound = 1000.0
            else:
                lower_bound = 0.0
                upper_bound = 1000.0
            # 创建反应对象
            reaction = Reaction(
                id=reaction_id,
                name=reaction_name,
                subsystem=reaction_subsystem,
                lower_bound=lower_bound,
                upper_bound=upper_bound
            )

            # 添加反应物和产物
            reaction_short = row["Reaction (short)"]
            # 2 C00027 <=> C00007 + 2 C00001 形如这样的反应 <=> 前面赋值负，后面赋值正
            reactants, products = reaction_short.split("<=>")
            # 处理反应物
            for reactant in reactants.split("+"):
                reactant = reactant.strip()
                # 解析系数和代谢物ID
                if ' ' in reactant:
                    coeff, metabolite_id = reactant.split(' ', 1)
                    coeff = float(coeff) * -1  # 反应物系数为负
                else:
                    coeff, metabolite_id = -1, reactant  # 默认系数为-1
                
                # 如果metabolite_id中不含有"_"，则添加"_c"，表示细胞质
                if not any(metabolite_id.endswith(suffix) for suffix in ["_c", "_e", "_b", "_p"]):
                    metabolite_id += "_c"
                metabolite = metabolites[metabolite_id]
                reaction.add_metabolites({metabolite: coeff})

            # 处理产物
            for product in products.split("+"):
                product = product.strip()
                # 解析系数和代谢物ID
                if ' ' in product:
                    coeff, metabolite_id = product.split(' ', 1)
                    coeff = float(coeff)  # 产物系数为正
                else:
                    coeff, metabolite_id = 1, product  # 默认系数为1
                    
                # 如果metabolite_id最后不含有"_"，则添加"_c"，表示细胞质
                if not any(metabolite_id.endswith(suffix) for suffix in ["_c", "_e", "_b", "_p"]):
                    metabolite_id += "_c"

                
                metabolite = metabolites[metabolite_id]
                reaction.add_metabolites({metabolite: coeff})

            # 添加反应规则
            if pd.isnull(row["Gene"]):
                pass
            else:
                reaction_rule = row["Gene"]
                # 将形如 NZ_AEYC01000083.1:41224..42981 or ilvH or NZ_AEYC01000001.1:283645..285390 的基因规则转化为 ‘( NZ_AEYC01000083.1:41224..42981 or ilvH or NZ_AEYC01000001.1:283645..285390 )’
                # 添加时先将complement(XX)重新编码为C-XX
                reaction_rule = re.sub(r'complement\((\S+)\)', r'C-\1', reaction_rule)
                reaction_rule = '(' + reaction_rule + ')'
                reaction.gene_reaction_rule = reaction_rule

            model.add_reactions([reaction])
            
            # 每添加10个反应保存一次
            if len(model.reactions) % 10 == 0:
                with open(TEMP_FILE, "wb") as f:
                    pickle.dump(model, f)
                    
        except Exception as e:
            print(f"处理反应 {reaction_id} 时出错: {str(e)}")
            continue

    # 保存最终模型
    clean_model_data(model)
    save_json_model(model, f"model/{model.id}.json")
    
    # 删除临时文件
    if os.path.exists(TEMP_FILE):
        os.remove(TEMP_FILE)
        
    print(f"模型构建完成，包含:")
    print(f"- {len(model.metabolites)} 个代谢物")
    print(f"- {len(model.reactions)} 个反应")
    
except Exception as e:
    # 发生错误时保存当前进度
    with open(TEMP_FILE, "wb") as f:
        pickle.dump(model, f)
    print(f"发生错误: {str(e)}")
    print(f"进度已保存到 {TEMP_FILE}")