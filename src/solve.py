from cobra import Configuration
from cobra.io import load_json_model
import pandas as pd

def set_production_objective(model, target_metabolite):
    """设置目标代谢物的最大化生产"""
    # 重置所有目标
    model.objective = {}
    
    # 设置新的目标
    try:
        target = model.metabolites.get_by_id(target_metabolite)
        # 创建输出反应
        export_reaction = model.reactions.get_by_id(f"R_Spinosad_b")
        model.objective = export_reaction
        print(f"目标设置为: {target.name} 的最大化生产")
    except KeyError:
        raise ValueError(f"未找到目标代谢物: {target_metabolite}")

def define_medium(model, medium_dict):
    """定义培养基条件"""
    # 首先关闭所有交换反应
    for reaction in model.exchanges:
        reaction.bounds = (0, 0)
    
    # 设置指定的培养基组分
    for metabolite, (lb, ub) in medium_dict.items():
        try:
            ex_reaction = model.reactions.get_by_id(f"EX_{metabolite}")
            ex_reaction.bounds = (lb, ub)
            print(f"设置培养基组分: {metabolite} [{lb}, {ub}]")
        except KeyError:
            print(f"警告: 未找到代谢物 {metabolite} 的交换反应")

def main():
    # 设置求解器
    Configuration().solver = "glpk"

    # 读取模型
    print("加载模型...")
    model = load_json_model("model/Saccharopolyspora_spinosa.json")
    print(f"模型加载完成，包含:")
    print(f"- {len(model.metabolites)} 个代谢物")
    print(f"- {len(model.reactions)} 个反应")
    print(f"- {len(model.genes)} 个基因")

    # 设置目标产物
    set_production_objective(model, 'C_Spinosad_b')

    # 定义培养基条件
    medium = {
        'glucose': (-10, 0),  # 葡萄糖摄入
        'o2': (-20, 0),      # 氧气摄入
        'nh4': (-5, 0),      # 氨摄入
        'pi': (-5, 0)        # 无机磷摄入
    }
    define_medium(model, medium)

    # 运行FBA
    try:
        solution = model.optimize()
        print("\n优化结果:")
        print(f"目标值: {solution.objective_value:.6f}")
        print(f"求解状态: {solution.status}")
        
        # 输出主要通量
        print("\n主要代谢物通量:")
        for ex in model.exchanges:
            if abs(solution.fluxes[ex.id]) > 1e-6:
                print(f"{ex.name}: {solution.fluxes[ex.id]:.6f}")
                
    except Exception as e:
        print(f"优化过程中出现错误: {str(e)}")

if __name__ == "__main__":
    main()