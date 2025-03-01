from cobra import Configuration
from cobra.io import load_json_model


def set_production_objective(model, target_metabolite):
    """设置目标代谢物的最大化生产"""
    # 重置所有目标
    model.objective = {}

    # 设置新的目标
    target_metabolite_id = target_metabolite + "_e"
    target_reaction_id = "EX_" + target_metabolite + "_e"
    if model.metabolites.get_by_id(target_metabolite_id):
        export_reaction = model.reactions.get_by_id(target_reaction_id)
        model.objective = export_reaction
        print(f"目标设置为: {target_metabolite} 的最大化生产")
    else:
        raise ValueError(f"未找到目标代谢物: {target_metabolite}")

# R_Spinosad	Spinosad = Spinosad_e	C_Spinosad <=> C_Spinosad_e
# EX_C_Spinosad_e                       C_Spinosad_e <=>

def define_medium(model, medium_dict):
    """定义培养基条件"""
    # 首先关闭所有交换反应
    for reaction in model.exchanges:
        reaction.bounds = (0, 0)
    
    # 设置指定的培养基组分
    for metabolite, (lb, ub) in medium_dict.items():
        try:
            ex_reaction = model.reactions.get_by_id(f"EX_{metabolite}_e")
            ex_reaction.bounds = (lb, ub)
            print(f"设置培养基组分: {metabolite} [{lb}, {ub}]")
        except KeyError:
            print(f"警告: 未找到代谢物 {metabolite} 的交换反应")


def main():
    # 设置求解器
    Configuration().solver = "glpk"

    # 读取模型
    print("加载模型...")
    model = load_json_model("../model/Saccharopolyspora_spinosa_boundaries.json")
    print(f"模型加载完成，包含:")
    print(f"- {len(model.metabolites)} 个代谢物")
    print(f"- {len(model.reactions)} 个反应")
    print(f"- {len(model.genes)} 个基因")

    # 设置目标产物
    set_production_objective(model, 'C_Spinosad')

    # 运行FBA
    try:
        solution = model.optimize()
        print("\n优化结果:")
        print(f"目标值: {solution.objective_value:.6f}")
        print(f"求解状态: {solution.status}")

        # 了解Spinosad的输入输出行为
        print("\nSpinosad生产情况:")
        print(model.metabolites.C_Spinosad_c.summary())
        print(model.metabolites.C_Spinosad_e.summary())

        # 了解主要的能量(C00002<->atp)生产和消耗反应
        print(model.metabolites.C00002_c.summary())

        # 输出主要通量
        print("\n主要代谢物通量:")
        for ex in model.exchanges:
            if abs(solution.fluxes[ex.id]) > 1e-6:
                print(f"{ex.name}: {solution.fluxes[ex.id]:.6f}")
                
    except Exception as e:
        print(f"优化过程中出现错误: {str(e)}")


if __name__ == "__main__":
    main()
