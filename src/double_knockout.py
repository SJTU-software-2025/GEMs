from cobra import Configuration
from cobra.io import load_json_model
from solve import set_production_objective
import pandas as pd
from cobra.flux_analysis import (double_gene_deletion, double_reaction_deletion)

# COBRApy: double_gene_deletion(cobra_model, cobra_model.genes[-5:]).round(4)
# .round(4)用于将数值四舍五入到小数点后四位的操作。它应用于 DataFrame 中的数值列。
# \binom{5}{2} + 4 = 10 + 4 = 14


def knock_out_double_genes_and_reactions(model):
    """
    COBRApy封装的方法，遍历进行反应和基因的双删除
    Unable to specify the reaction ID or gene ID
    """
    double_reaction_deletion_results = double_reaction_deletion(model)
    double_reaction_deletion_results.to_csv('../data/double_reaction_deletion_results.csv')
    print(f"结果已保存到: '../data/double_reaction_deletion_results.csv'")
    double_gene_deletion_results = double_gene_deletion(model)
    double_gene_deletion_results.to_csv('../data/double_gene_deletion_results.csv')
    print(f"结果已保存到: '../data/double_gene_deletion_results.csv'")


def find_essential_double_reaction(model, production_threshold=0.01):
    """
    Find essential double reaction
    param production_threshold: The minimum production change threshold to consider a reaction essential.
    return: A list of essential two combinations of reactions.
    """
    essential_reactions = []
    objective_value = model.optimize().objective_value
    double_reaction_deletion_results = pd.read_csv('../data/double_reaction_deletion_results.csv')
    for _, row in double_reaction_deletion_results.iterrows():
        if row['status'] != 'optimal':
            continue
        if abs((row['growth'] - objective_value) / objective_value) >= production_threshold:
            essential_reactions.append(row['ids'])
    return essential_reactions


def find_essential_double_gene(model, production_threshold=0.01):
    """
    Find essential double gene
    param production_threshold: The minimum production change threshold to consider a gene essential.
    return: A list of essential two combinations of genes
    """
    essential_genes = []
    objective_value = model.optimize().objective_value
    double_gene_deletion_results = pd.read_csv('../data/double_gene_deletion_results.csv')
    for _, row in double_gene_deletion_results.iterrows():
        if row['status'] != 'optimal':
            continue
        if abs((row['growth'] - objective_value) / objective_value) >= production_threshold:
            essential_genes.append(row['ids'])
    return essential_genes


def main():
    Configuration().solver = "glpk"
    model = load_json_model("../model/Saccharopolyspora_spinosa_boundaries.json")
    set_production_objective(model, 'C_Spinosad')

    knock_out_double_genes_and_reactions(model)

    essential_reactions = find_essential_double_reaction(model)
    print(f"Essential reactions: {essential_reactions}")
    essential_genes = find_essential_double_gene(model)
    print(f"Essential genes: {essential_genes}")


if __name__ == "__main__":
    main()
