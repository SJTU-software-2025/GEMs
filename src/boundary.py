# 定义边界反应
from cobra import Metabolite, Reaction, Model
from cobra.io import load_json_model, save_json_model
import tqdm

def main():
    # 读取模型
    model = load_json_model("model/Saccharopolyspora_spinosa.json")
    # 遍历所有的代谢物，如果代谢物以 _e 结尾，就添加一个边界反应
    for metabolite in tqdm.tqdm(model.metabolites, desc="Adding boundary reactions"):
        if metabolite.id.endswith("_e"):
            print(f"Adding boundary reaction for {metabolite.id}")
            model.add_boundary(metabolite, type="exchange")
    # 保存模型
    save_json_model(model, "model/Saccharopolyspora_spinosa_boundaries.json")

if __name__ == "__main__":
    main()