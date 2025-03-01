# 定义边界反应
from cobra.io import load_json_model, save_json_model
from tqdm import tqdm
from cobra import Configuration
Configuration().solver = "glpk"


def main():
    # 读取模型
    model = load_json_model("../model/Saccharopolyspora_spinosa.json")
    # 遍历所有的代谢物，如果代谢物以 _e 结尾，就添加一个边界反应
    for metabolite in tqdm(model.metabolites, desc="Adding boundary reactions"):
        if metabolite.id.endswith("_e"):
            print(f"Adding boundary reaction for {metabolite.id}")
            model.add_boundary(metabolite, type="exchange")
    # 保存模型
    save_json_model(model, "../model/Saccharopolyspora_spinosa_boundaries.json")
# metabolite：C_Spinosad_e
# exchange reaction：EX_C_Spinosad_e

if __name__ == "__main__":
    main()
