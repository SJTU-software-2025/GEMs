from cobra import Configuration
from cobra.io import load_json_model
from solve import set_production_objective
import pandas as pd
from cobra.flux_analysis import (single_gene_deletion, single_reaction_deletion)

# 创建一个空的 DataFrame 用于存储反应ID和优化结果
knockout_single_reaction_result = []
knockout_single_gene_results = []


def knock_out_single_reaction(model, reaction_id):
    """
    Knock out single reaction
    Able to specify the reaction ID
    """
    with model:
        model.reactions.get_by_id(reaction_id).knock_out()
        solution = model.optimize()
        knockout_single_reaction_result.append({'Reaction ID': reaction_id, 'Objective Value': solution.objective_value})


def knock_out_single_gene(model, gene_id):
    """
    Knock out single gene
    Able to specify the gene ID
    """
    with model:
        model.genes.get_by_id(gene_id).knock_out()
        solution = model.optimize()
        knockout_single_gene_results.append({'Gene ID': gene_id, 'Objective Value': solution.objective_value})


def knock_out_single_genes_and_reactions(model):
    """
    COBRApy封装的方法，遍历计算单个反应和基因的删除效果
    Unable to specify the reaction ID or gene ID
    """
    single_reaction_deletion_results = single_reaction_deletion(model)
    single_reaction_deletion_results.to_csv('../data/single_reaction_deletion_results.csv')
    print(f"结果已保存到: '../data/single_reaction_deletion_results.csv'")
    single_gene_deletion_results = single_gene_deletion(model)
    single_gene_deletion_results.to_csv('../data/single_gene_deletion_results.csv')
    print(f"结果已保存到: '../data/single_gene_deletion_results.csv'")


def find_essential_single_reaction(model, production_threshold=0.01):
    """
    Find essential single reaction
    param production_threshold: The minimum production change threshold to consider a reaction essential.
    return: A list of essential reactions.
    """
    essential_reactions = []
    objective_value = model.optimize().objective_value
    single_reaction_deletion_results = pd.read_csv('../data/single_reaction_deletion_results.csv')
    for _, row in single_reaction_deletion_results.iterrows():
        if row['status'] != 'optimal':
            continue
        if abs((row['growth'] - objective_value) / objective_value) >= production_threshold:
            essential_reactions.append(row['ids'])
    return essential_reactions


def find_essential_single_gene(model, production_threshold=0.01):
    """
    Find essential single gene
    param production_threshold: The minimum production change threshold to consider a gene essential.
    return: A list of essential genes
    """
    essential_genes = []
    objective_value = model.optimize().objective_value
    single_gene_deletion_results = pd.read_csv('../data/single_gene_deletion_results.csv')
    for _, row in single_gene_deletion_results.iterrows():
        if row['status'] != 'optimal':
            continue
        if abs((row['growth'] - objective_value) / objective_value) >= production_threshold:
            essential_genes.append(row['ids'])
    return essential_genes


def main():
    Configuration().solver = "glpk"
    model = load_json_model("../model/Saccharopolyspora_spinosa_boundaries.json")
    set_production_objective(model, 'C_Spinosad')

    for reaction in model.reactions:
        knock_out_single_reaction(model, reaction.id)

    # 将结果转换为 DataFrame
    knockout_single_reaction_df_result = pd.DataFrame(knockout_single_reaction_result)
    # 将 DataFrame 存储为 CSV 文件
    knockout_single_reaction_df_result.to_csv('../data/knockout_single_reaction_result.csv', index=False)
    print(f"结果已保存到: '../data/knockout_single_reaction_result.csv'")

    for gene in model.genes:
        knock_out_single_gene(model, gene.id)

    knockout_single_gene_df_results = pd.DataFrame(knockout_single_gene_results)
    knockout_single_gene_df_results.to_csv('../data/knockout_single_gene_result.csv', index=False)
    print(f"结果已保存到: '../data/knockout_single_gene_result.csv'")

    knock_out_single_genes_and_reactions(model)

    essential_reactions = find_essential_single_reaction(model)
    print(f"Essential reactions: {essential_reactions}")
    essential_genes = find_essential_single_gene(model)
    print(f"Essential genes: {essential_genes}")


if __name__ == "__main__":
    main()
