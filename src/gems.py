import cobra
from cobra import Model, Reaction, Metabolite
import pandas as pd

# 代谢物和反应数据文件路径
metabolite_file = "../data/metabolites.csv"  # 代谢物文件路径
reaction_file = "../data/reactions.csv"  # 反应文件路径

# 读取代谢物数据
metabolite_data = pd.read_csv(metabolite_file)

# 构建代谢物ID到名称的映射，并去掉 "_c" 后缀
metabolite_mapping = {}
for _, row in metabolite_data.iterrows():
    metabolite_id = row["ID"].replace("_c", "")  # 去掉后缀 "_c"，_c表示细胞质，但反应中默认没有_c了
    metabolite_name = row["Name"]
    compartment = row["Compartment"]
    
    # 存储代谢物映射关系
    metabolite_mapping[metabolite_id] = {
        "name": metabolite_name,
        "compartment": compartment
    }

# 报错：Gurobi 许可证已过期，导致无法创建模型。cobrapy 默认会尝试使用 Gurobi 求解器
# conda install -c conda-forge glpk
# 切换到免费的求解器，如 GLPK

from cobra import Configuration
Configuration().solver = "glpk"

# 初始化模型
model = Model("Integrated_Model")

# 将代谢物添加到模型中
for met_id, met_info in metabolite_mapping.items():
    metabolite = Metabolite(
        id=met_id,
        name=met_info["name"],
        compartment=met_info["compartment"]
    )
    model.add_metabolites([metabolite])

# 读取反应数据
reaction_data = pd.read_csv(reaction_file)

# 遍历反应数据，构建并添加反应
for _, row in reaction_data.iterrows():
    reaction_name = row['Reaction ID'].replace(' ', '-') #-?

    # 创建Reaction对象
    reaction = Reaction(reaction_name)
    reaction.name = reaction_name
    
    # 设置反应是否可逆
    reaction.lower_bound = -1000.0 if row['Reversible'] == 'TRUE' else 0.0
    reaction.upper_bound = 1000.0
    
    # 解析化学反应式，添加底物和产物
    reactants_products = row['Reaction (short)'].split('<=>')
    reactants = reactants_products[0].strip().split('+')  # 底物
    products = reactants_products[1].strip().split('+')  # 产物
    
    # 添加底物
    for reactant in reactants:
        reactant = reactant.strip()  # 去掉空格
        # 检查是否包含化学计量系数（如 2 C00013）
        if ' ' in reactant:
            coeff, r_id = reactant.split(' ', 1) #1?
            coeff = float(coeff)  # 转换为数值
        else:
            coeff, r_id = 1.0, reactant  # 如果没有化学计量系数，默认为1
        
        # 去掉代谢物ID中的 "_c" 后缀
        r_id = r_id.replace("_c", "") #?
        
        # 获取代谢物对象
        metabolite = model.metabolites.get_by_id(r_id)
        reaction.add_metabolites({metabolite: -coeff})  # 负值表示底物
    
    # 添加产物
    for product in products:
        product = product.strip()  # 去掉空格
        # 检查是否包含化学计量系数（如 2 C00009）
        if ' ' in product:
            coeff, p_id = product.split(' ', 1)
            coeff = float(coeff)  # 转换为数值
        else:
            coeff, p_id = 1.0, product  # 如果没有化学计量系数，默认为1
        
        # 去掉代谢物ID中的 "_c" 后缀
        p_id = p_id.replace("_c", "")
        
        # 获取代谢物对象
        metabolite = model.metabolites.get_by_id(p_id)
        reaction.add_metabolites({metabolite: coeff})  # 正值表示产物
    
    # 添加到模型
    if len(reaction.id) >= 256: #why?
        print(f"反应ID过长，已截断：{reaction.id}")
        reaction.id = reaction.id[:100]

    gene_reaction_rule = row['Gene'] if not pd.isna(row['Gene']) else ""

    if gene_reaction_rule != "":
        # 如果基因反应规则包含 'or'，则遍历每个基因并替换 complement()
        if 'or' in gene_reaction_rule:
            # 遍历每个基因（以 'or' 为分隔符）
            genes = gene_reaction_rule.split(' or ')
            # 对每个基因进行处理，替换 complement()
            for i, gene in enumerate(genes):  # 使用 enumerate 来获取索引和元素
                if gene.startswith('complement'):
                    genes[i] = 'C-' + gene[11:-1]  # 直接更新列表中的元素
            # 将处理后的基因重新连接成完整的基因反应规则
            gene_reaction_rule = ' or '.join(genes)
        else:
            if gene_reaction_rule.startswith('complement'):
                gene_reaction_rule = 'C-' + gene_reaction_rule[11:-1]

    # 将最终的基因反应规则赋值给 reaction.gene_reaction_rule
    reaction.gene_reaction_rule = gene_reaction_rule

    model.add_reactions([reaction])
    
    # # 可选：附加信息（基因、酶、子系统）
    # reaction.gene_reaction_rule = row['Gene'] if not pd.isna(row['Gene']) else ""
    # reaction.annotation = {
    #     "enzyme": row['Enzyme'],
    #     "subsystem": row['Subsystem']
    # }

# 保存模型
cobra.io.save_json_model(model, "../model/integrated_model.json")
print("模型构建完成并保存为 integrated_model.json")
