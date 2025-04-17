"""
BRENDA: TURNOVER NUMBER [1/s] 	SUBSTRATE
Sabio-RK:

https://cobrapy.readthedocs.io/en/latest/media.html:
fluxes have the unit mmol / [gDW h] (concentration per gram dry weight of cells and hour)

1kDa表示相对分子质量为1000的分子
酶的分子量为 50kDa = 50000g/mol(典型值)
假设总酶量占细胞干重的 40% (即 0.4g酶/gDW)。
总酶池约束为 0.008 mmol酶/gDW
"""


from cobra.io import load_model
from cobra import Configuration

# 设置求解器
Configuration().solver = "glpk"

model = load_model("textbook")

for reaction in model.reactions:
    # 每个反应对应一个酶
    # ub设置为总酶池约束
    enzyme_var = model.problem.Variable(f"Enzyme_{reaction.id}", lb=0, ub=8)
    model.add_cons_vars(enzyme_var)

    # 中位数约0.5（单位：1/s）
    kcat = 1800 # (单位：1/h)
    constraint_max = model.problem.Constraint(
        reaction.flux_expression - kcat * enzyme_var,
        ub=0,  # v_i - kcat*E_i ≤ 0 → v_i ≤ kcat*E_i
        name=f"enzyme_constraint_{reaction.id}_max"
    )
    constraint_min = model.problem.Constraint(
        reaction.flux_expression + kcat * enzyme_var,
        lb=0,  # v_i + kcat*E_i ≥ 0 → v_i ≥ -kcat*E_i
        name=f"enzyme_constraint_{reaction.id}_min"
    )
    model.add_cons_vars([constraint_max, constraint_min])

total_protein_content = 8  # 单位：0.008mmol/gDW（假设值）× 1000 = 8mmol/gDW (弥补求解器精度不足)
enzyme_vars = [ var for var in model.variables if var.name.startswith("Enzyme_")]

total_enzyme_expr = sum(enzyme_vars)
total_enzyme_constraint = model.problem.Constraint(
    total_enzyme_expr,
    lb=0,
    ub=total_protein_content,
    name="total_enzyme_constraint"
)
model.add_cons_vars(total_enzyme_constraint)

solution = model.optimize()
print(solution.fluxes)
print(model.summary())
