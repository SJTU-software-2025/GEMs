from cobra.io import load_json_model
from cobra import Reaction, Metabolite
from escher import Builder
import pandas as pd
import re
from pathlib import Path
import requests
import time
from collections import Counter
from tqdm import tqdm
import json

# ai写的，更新代谢物信息用，现在还有bug，回头改
# 加载模型
model = load_json_model("integrated_model.json")

# 1. 找出主干反应
def analyze_core_reactions(model):
    # 获取所有反应的通量
    solution = model.optimize()
    active_reactions = [rxn.id for rxn in model.reactions 
                       if abs(solution.fluxes[rxn.id]) > 1e-10]
    return active_reactions

# 2. GPR分析
def analyze_gpr(model):
    gpr_info = []
    for reaction in model.reactions:
        if reaction.genes:
            gpr_info.append({
                'Reaction': reaction.id,
                'Genes': ', '.join([g.id for g in reaction.genes]),
                'GPR_rule': reaction.gene_reaction_rule
            })
    return pd.DataFrame(gpr_info)

# 3. 更新代谢物信息
def update_metabolite_info(model, metabolite_id, formula, charge, molecular_weight):
    if metabolite_id in model.metabolites:
        met = model.metabolites.get_by_id(metabolite_id)
        met.formula = formula
        met.charge = charge
        met.molecular_weight = molecular_weight

# 4. 设置目标反应（多杀菌素生产）
def set_production_objective(model, product_id):
    # 假设 product_id 是多杀菌素的输出反应
    model.objective = product_id
    model.objective_direction = 'max'

# 5. 定义交换反应
def define_medium(model, medium_components):
    # medium_components 是一个字典，包含培养基成分和其上下限
    for metabolite, bounds in medium_components.items():
        if f'EX_{metabolite}' in model.reactions:
            exchange_rxn = model.reactions.get_by_id(f'EX_{metabolite}')
            exchange_rxn.bounds = bounds

def load_metabolite_data():
    """加载代谢物数据并返回DataFrame"""
    metabolites_file = Path('data/metabolites.csv')
    if metabolites_file.exists():
        return pd.read_csv(metabolites_file)
    return None

def fetch_kegg_compound(compound_id):
    """从KEGG API获取化合物信息"""
    base_url = "http://rest.kegg.jp/get/cpd:"
    try:
        # 移除'C'前缀并补充0到5位数
        compound_number = compound_id.lstrip('C').zfill(5)
        
        # 添加超时设置
        response = requests.get(f"{base_url}{compound_number}", timeout=5)
        
        if response.status_code == 200:
            data = response.text
            # 解析KEGG响应
            formula = None
            exact_mass = None
            mol_weight = None
            charge = 0  # 默认电荷为0
            
            for line in data.split('\n'):
                if line.startswith('FORMULA'):
                    formula = line.split()[1]
                elif line.startswith('EXACT_MASS'):
                    exact_mass = float(line.split()[1])
                elif line.startswith('MOL_WEIGHT'):
                    mol_weight = float(line.split()[1])
            
            return {
                'formula': formula,
                'weight': mol_weight or exact_mass,
                'charge': charge
            }
        elif response.status_code == 404:
            print(f"\n代谢物 {compound_id} 在KEGG数据库中未找到")
            return None
    except requests.Timeout:
        print(f"\n获取 {compound_id} 时请求超时")
    except Exception as e:
        print(f"\n获取KEGG信息时出错 ({compound_id}): {str(e)}")
    return None

def validate_chemical_formula(formula):
    """验证化学式的正确性"""
    if not formula:
        return False
    
    # 基本格式验证
    pattern = r'^([A-Z][a-z]?\d*)*$'
    if not re.match(pattern, formula):
        return False
    
    # 解析化学式
    elements = {}
    element_pattern = r'([A-Z][a-z]?)(\d*)'
    for element, count in re.findall(element_pattern, formula):
        count = int(count) if count else 1
        if element not in ['C', 'H', 'O', 'N', 'P', 'S', 'Fe', 'Co', 'Cu', 'Zn', 'Mg', 'Ca', 'Na', 'K', 'Cl']:
            return False
        elements[element] = elements.get(element, 0) + count
    
    # 基本规则验证
    if not elements:
        return False
    
    return True

def update_model_metabolites(model, metabolites_df):
    """更新模型中的代谢物信息"""
    # 创建缓存目录
    cache_dir = Path('cache')
    cache_dir.mkdir(exist_ok=True)
    cache_file = cache_dir / 'metabolites_cache.json'
    
    def save_cache():
        """保存当前的代谢物信息到缓存文件"""
        with open(cache_file, 'w', encoding='utf-8') as f:
            json.dump(common_metabolites, f, indent=2)
    
    def is_metabolite_complete(met):
        """检查代谢物信息是否完整"""
        return all([
            hasattr(met, 'formula') and met.formula,
            hasattr(met, 'charge') and met.charge is not None,
            hasattr(met, 'molecular_weight') and met.molecular_weight
        ])
    
    # 加载缓存的代谢物信息
    if cache_file.exists():
        with open(cache_file, 'r', encoding='utf-8') as f:
            common_metabolites = json.load(f)
    else:
        common_metabolites = {
            'C00001': {'formula': 'H2O', 'charge': 0, 'weight': 18.015},
            'C00002': {'formula': 'C10H16N5O13P3', 'charge': -4, 'weight': 507.18},
            'C00003': {'formula': 'C21H28N7O14P2', 'charge': -1, 'weight': 663.43},
            'C00004': {'formula': 'C21H29N7O14P2', 'charge': -2, 'weight': 664.44},
            'C00005': {'formula': 'C21H30N7O17P3', 'charge': -4, 'weight': 744.41},
            'C00006': {'formula': 'C21H29N7O17P3', 'charge': -3, 'weight': 743.41},
            'C00007': {'formula': 'O2', 'charge': 0, 'weight': 31.998},
            'C00008': {'formula': 'C10H15N5O10P2', 'charge': -3, 'weight': 427.20},
            'C00009': {'formula': 'H3PO4', 'charge': -2, 'weight': 97.977},
            'C00010': {'formula': 'C21H36N7O16P3S', 'charge': -4, 'weight': 767.53}
        }
    
    updated_count = 0
    kegg_updated_count = 0
    error_count = 0
    skipped_count = 0
    save_interval = 10  # 每更新10个代谢物保存一次
    
    # 首先统计需要从KEGG获取的代谢物
    need_kegg_query = []
    for met in model.metabolites:
        kegg_id = met.id.split('_')[0]
        if not is_metabolite_complete(met) and kegg_id not in common_metabolites:
            need_kegg_query.append((met, kegg_id))
    
    total_kegg_queries = len(need_kegg_query)
    print(f"\n需要从KEGG获取的代谢物数量: {total_kegg_queries}")
    
    # 使用两个进度条
    main_pbar = tqdm(total=len(model.metabolites), desc="总体进度")
    kegg_pbar = tqdm(total=total_kegg_queries, desc="KEGG查询进度", leave=False)
    
    try:
        # 首先处理已有信息的代谢物
        for met in model.metabolites:
            kegg_id = met.id.split('_')[0]
            
            # 检查是否已有完整信息
            if is_metabolite_complete(met):
                skipped_count += 1
                main_pbar.update(1)
                continue
            
            # 尝试从缓存中获取信息
            if kegg_id in common_metabolites:
                info = common_metabolites[kegg_id]
                if validate_chemical_formula(info['formula']):
                    met.formula = info['formula']
                    met.charge = info['charge']
                    met.molecular_weight = info['weight']
                    updated_count += 1
                else:
                    error_count += 1
                    print(f"\n警告: 代谢物 {kegg_id} 的化学式无效: {info['formula']}")
                main_pbar.update(1)
        
        # 然后处理需要从KEGG获取的代谢物
        for met, kegg_id in need_kegg_query:
            # 更新两个进度条的描述
            main_pbar.set_postfix({
                '本地更新': updated_count,
                'KEGG更新': kegg_updated_count,
                '已跳过': skipped_count,
                '错误': error_count
            })
            kegg_pbar.set_postfix({'当前查询': kegg_id})
            
            # 从KEGG获取信息
            kegg_info = fetch_kegg_compound(kegg_id)
            if kegg_info and validate_chemical_formula(kegg_info['formula']):
                met.formula = kegg_info['formula']
                met.charge = kegg_info['charge']
                met.molecular_weight = kegg_info['weight']
                kegg_updated_count += 1
                common_metabolites[kegg_id] = kegg_info
            
            # 更新进度条
            main_pbar.update(1)
            kegg_pbar.update(1)
            
            # 定期保存缓存
            if (updated_count + kegg_updated_count) % save_interval == 0:
                save_cache()
            
            time.sleep(0.1)  # 避免过于频繁请求KEGG API
                
    except KeyboardInterrupt:
        print("\n检测到程序中断，正在保存当前进度...")
        save_cache()
        raise
    except Exception as e:
        print(f"\n发生错误: {str(e)}")
        save_cache()
        raise
    finally:
        # 确保在结束时保存缓存
        save_cache()
        main_pbar.close()
        kegg_pbar.close()
    
    print(f"\n更新统计:")
    print(f"- 从本地更新: {updated_count} 个代谢物")
    print(f"- 从KEGG更新: {kegg_updated_count} 个代谢物")
    print(f"- 已跳过: {skipped_count} 个代谢物")
    print(f"- 出错: {error_count} 个代谢物")
    
    return updated_count + kegg_updated_count

# 主要分析流程
def main():
    # 1. 加载模型
    print("正在加载模型...")
    model = load_json_model("integrated_model.json")
    
    # 2. 加载代谢物数据
    print("正在加载代谢物数据...")
    metabolites_df = load_metabolite_data()
    
    # 3. 更新代谢物信息
    print("正在更新代谢物信息...")
    updated_count = update_model_metabolites(model, metabolites_df)
    print(f"已更新 {updated_count} 个代谢物的信息")
    
    # 4. 分析主干反应
    core_reactions = analyze_core_reactions(model)
    print(f"发现 {len(core_reactions)} 个主干反应")
    
    # 5. GPR分析
    gpr_df = analyze_gpr(model)
    gpr_df.to_csv('gpr_analysis.csv', index=False)
    print("GPR分析已保存到 gpr_analysis.csv")
    
    # 6. 设置多杀菌素生产为目标
    set_production_objective(model, 'C_Spinosad_e')
    
    # 7. 设置培养基条件
    medium = {
        'glucose': (-10, 0),
        'o2': (-20, 0),
        'nh4': (-5, 0),
        'pi': (-5, 0)
    }
    define_medium(model, medium)
    
    # 8. 优化模型
    try:
        solution = model.optimize()
        print(f"多杀菌素产量: {solution.objective_value}")
    except Exception as e:
        print(f"优化过程中出现错误: {str(e)}")
    
    # 9. 保存更新后的模型
    from cobra.io import save_json_model
    save_json_model(model, "updated_model.json")
    print("更新后的模型已保存为 updated_model.json")
    
    # 10. 可视化
    builder = Builder(
        map_name="e_coli_core.Core metabolism",
        model=model
    )
    builder.save_html('metabolic_network.html')
    print("可视化已保存为 metabolic_network.html")

if __name__ == "__main__":
    main()