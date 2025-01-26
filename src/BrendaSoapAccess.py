"""
EC number & Substrate -> Turnover number API
Reference: https://www.brenda-enzymes.org/soap.php
运行：python3 BrendaSoapAccess.py
"""
#!/usr/bin/python
from zeep import Client
import hashlib
import pandas as pd

wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
# Web服务描述语言文件的URL，指明了SOAP服务的结构。

password = hashlib.sha256("Zrp031030".encode("utf-8")).hexdigest()
# 使用SHA-256加密对密码进行加密

client = Client(wsdl)
# 使用Zeep库初始化一个SOAP客户端

"""
API test
"""

# parameters = ( "Zrp1030@sjtu.edu.cn",password,"ecNumber*1.1.1.1","organism*Homo sapiens","kmValue*",
#               "kmValueMaximum*","substrate*","commentary*","ligandStructureId*","literature*" )
# # 定义SOAP请求中需要传递的参数
#
# resultArray = client.service.getKmValue(*parameters)
# # 调用getKmValue方法并获取结果
#
# print (resultArray)

"""
主函数
"""

# 为什么要搜两次再把结果合并：
# 假设反应A+B->C+D 会出现这样的情况：使用EC编号+SubstrateA搜不到结果，使用EC编号+SubstrateB能搜到结果

# 现有一个dataframe，有两列，分别叫TurnoverNumber_Reactant1AsSubstrate和TurnoverNumber_Reactant2AsSubstrate，合并这两列为新的一列：TurnoverNumber
# 暂定的合并规则为：如果第一列非NaN，则取第一列的值，如何第一列是NaN，第二列有值，则取第二列的值，如何都是NaN，则结果也为NaN

reaction_file = "../data/reactions.csv"
reaction_data = pd.read_csv(reaction_file)
reaction_data_test = reaction_data.copy()

turnover_numbers_Reactant1 = []

# 遍历reaction_data["Enzyme"]中的每一个EC编号
for _, row in reaction_data.iterrows():
    ec_number = row['Enzyme']
    # 判断 ec_number 是否为 NaN
    if pd.isna(ec_number):
        turnover_numbers_Reactant1.append(None)
        continue

    # 一个EC编号的最大长度：X.XX.XXX.XXX 小于等于12个字符
    # https://www.brenda-enzymes.org/all_enzymes.php
    # 部分反应有多个酶催化，如何处理？
    # 此处先跳过
    if len(ec_number) > 12:
        turnover_numbers_Reactant1.append(None)
        continue

    # 未处理可逆反应
    # 当反应物有多个时，如何处理？
    # 此处暂时取第一个反应物作为底物
    reactants_products = row['Reaction (Long)'].split(' = ')
    reactants = reactants_products[0].split(' + ')

    substrate = reactants[0]
    # 检查是否包含化学计量系数（如 2 Pyruvate）
    if ' ' in substrate and row['Reaction ID'] != 'R01488' and row['Reaction ID'] != 'R01724':
        parts = substrate.split(' ', 1)
        # 检查第一个空格前的部分是否是数字
        # 如trans,trans-Farnesyl diphosphate和2 Geranylgeranyl diphosphate
        if parts[0].isdigit():  # 如果是数字
            substrate = parts[1]  # 反应物名称
        else:
            substrate = substrate  # 反应物名称直接为整个字符串

    """
    热补丁存在特殊情况的行
    """
    if row['Reaction ID'] == 'R01488':
        substrate = 'NADP+'
    if row['Reaction ID'] == 'R01724':
        substrate = 'Nicotinate'

    parameters = (
        "Zrp1030@sjtu.edu.cn",
        password,  # 密码
        f"ecNumber*{ec_number}",  # 将EC编号加入查询参数
        "turnoverNumber*",
        "turnoverNumberMaximum*",
        f"substrate*{substrate}",  # 将底物加入查询参数
        "commentary*",
        "organism*",
        "ligandStructureId*",
        "literature*"
    )

    resultArrayReactant1 = client.service.getTurnoverNumber(*parameters)

    # 从结果中提取turnoverNumber并添加到turnover_numbers列表中
    # 形如：
    # 4.1.1.47 0.25 glyoxylate mutant V51D
    # 4.1.1.47 2.9 glyoxylate mutant V51E
    # 4.1.1.47 18.5 glyoxylate mutant V51S
    # 4.1.1.47 18.9 glyoxylate wild-type enzyme
    # 4.1.1.47 19.7 glyoxylate mutant E520
    # 这样的搜索结果如何处理？
    # 此处暂时取第一个结果
    if resultArrayReactant1:
        turnover_number = resultArrayReactant1[0].turnoverNumber if hasattr(resultArrayReactant1[0],
                                                                            'turnoverNumber') else None  # 获取第一个结果的turnoverNumber
    else:
        turnover_number = None  # 如果没有返回结果，则赋值为None

    turnover_numbers_Reactant1.append(turnover_number)

reaction_data_test['TurnoverNumber_Reactant1AsSubstrate'] = turnover_numbers_Reactant1
print(reaction_data_test['TurnoverNumber_Reactant1AsSubstrate'].count())

turnover_numbers_Reactant2 = []  # 修改为新的列表名

# 遍历 reaction_data["Enzyme"] 中的每一个 EC 编号
for _, row in reaction_data.iterrows():
    ec_number = row['Enzyme']

    # 判断 ec_number 是否为 NaN
    if pd.isna(ec_number):
        turnover_numbers_Reactant2.append(None)
        continue

    # 一个 EC 编号的最大长度：X.XX.XXX.XXX 小于等于 12 个字符
    # https://www.brenda-enzymes.org/all_enzymes.php
    # 部分反应有多个酶催化，如何处理？
    # 此处先跳过
    if len(ec_number) > 12:
        turnover_numbers_Reactant2.append(None)
        continue

    # 未处理可逆反应
    # 当反应物有多个时，如何处理？
    # 此处暂时取第二个反应物作为底物（如果反应物大于等于2个）
    reactants_products = row['Reaction (Long)'].split(' = ')
    reactants = reactants_products[0].split(' + ')

    # 选择第二个反应物作为底物（如果反应物大于等于2个）
    if len(reactants) >= 2:
        substrate = reactants[1]
    else:
        substrate = reactants[0]  # 如果只有一个反应物，选择第一个反应物作为底物

    # 检查是否包含化学计量系数（如 2 Pyruvate）
    if ' ' in substrate and row['Reaction ID'] != 'R01488' and row['Reaction ID'] != 'R01724':
        parts = substrate.split(' ', 1)
        # 检查第一个空格前的部分是否是数字
        if parts[0].isdigit():  # 如果是数字
            substrate = parts[1]  # 反应物名称
        else:
            substrate = substrate  # 反应物名称直接为整个字符串

    # 热补丁存在特殊情况的行
    if row['Reaction ID'] == 'R01488':
        substrate = 'NADP+'
    if row['Reaction ID'] == 'R01724':
        substrate = 'Nicotinate'

    # 设置参数
    parameters = (
        "Zrp1030@sjtu.edu.cn",
        password,  # 密码
        f"ecNumber*{ec_number}",  # 将 EC 编号加入查询参数
        "turnoverNumber*",
        "turnoverNumberMaximum*",
        f"substrate*{substrate}",  # 将底物加入查询参数
        "commentary*",
        "organism*",
        "ligandStructureId*",
        "literature*"
    )

    resultArrayReactant2 = client.service.getTurnoverNumber(*parameters)  # 修改为 resultArrayReactant2

    # 从结果中提取 turnoverNumber 并添加到 turnover_numbersReactant2 列表中
    if resultArrayReactant2:
        turnover_number = resultArrayReactant2[0].turnoverNumber if hasattr(resultArrayReactant2[0],
                                                                            'turnoverNumber') else None
    else:
        turnover_number = None  # 如果没有返回结果，则赋值为 None

    turnover_numbers_Reactant2.append(turnover_number)  # 修改为 turnover_numbersReactant2

reaction_data_test['TurnoverNumber_Reactant2AsSubstrate'] = turnover_numbers_Reactant2
reaction_data_test['TurnoverNumber_Reactant2AsSubstrate'].count()

# 将列表转换为 DataFrame
df = pd.DataFrame(turnover_numbers_Reactant1, columns=['TurnoverNumber'])
# 将 DataFrame 存储为 Excel 文件
df.to_excel('../data/使用第一个反应物+EC编号搜索到的酶kcat值.xlsx', index=False, engine='openpyxl')  # 不保存索引
df2 = pd.DataFrame(turnover_numbers_Reactant2, columns=['TurnoverNumber'])
df2.to_excel('../data/使用第二个反应物+EC编号搜索到的酶kcat值.xlsx', index=False, engine='openpyxl')  # 不保存索引

# 使用 combine_first 来合并两列
reaction_data_test['TurnoverNumber'] = reaction_data_test['TurnoverNumber_Reactant1AsSubstrate'].combine_first(reaction_data_test['TurnoverNumber_Reactant2AsSubstrate'])
print(reaction_data_test['TurnoverNumber'].count())

reaction_data_test = reaction_data_test.drop(columns=['TurnoverNumber_Reactant1AsSubstrate', 'TurnoverNumber_Reactant2AsSubstrate'], axis=1)
reaction_data_test.to_csv('../data/reactions+kcat.csv', index=False)
print("Done!")