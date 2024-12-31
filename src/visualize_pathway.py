from cobra.io import load_json_model
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from pathlib import Path
import plotly.graph_objects as go
import pandas as pd


# 这个文件是ai写的，用来可视化代谢通路，检查得不是很仔细，临时用一用，后面还是会用escher来画的
# 常见的辅因子和小分子代谢物
COMMON_METABOLITES = {
    'C00001',  # H2O
    'C00080',  # H+
    'C00007',  # O2
    'C00011',  # CO2
    'C00013',  # Diphosphate
    'C00009',  # Orthophosphate
    'C00008',  # ADP
    'C00002',  # ATP
    'C00003',  # NAD+
    'C00004',  # NADH
    'C00005',  # NADPH
    'C00006',  # NADP+
    'C00010',  # CoA
}

def create_metabolic_network(model, remove_common_metabolites=True):
    """
    创建代谢网络图，支持过滤常见代谢物
    
    参数:
    - model: cobra模型
    - remove_common_metabolites: 是否移除常见代谢物
    
    返回:
    - G: networkx图对象
    - node_types: 节点类型字典
    """
    # 创建有向图
    G = nx.DiGraph()
    node_types = {}
    
    # 添加所有反应和代谢物
    for reaction in model.reactions:
        # 跳过交换反应和生物量反应
        if reaction.id.startswith('EX_') or reaction.id.startswith('DM_') or reaction.id.startswith('BIOMASS_'):
            continue
            
        # 添加反应节点
        G.add_node(reaction.id, type='reaction')
        node_types[reaction.id] = 'reaction'
        
        # 添加代谢物节点和边
        for substrate in reaction.reactants:
            # 如果是常见代谢物且选择移除，则跳过
            if remove_common_metabolites and substrate.id.split('_')[0] in COMMON_METABOLITES:
                continue
                
            G.add_node(substrate.id, type='metabolite')
            node_types[substrate.id] = 'metabolite'
            G.add_edge(substrate.id, reaction.id)
        
        for product in reaction.products:
            # 如果是常见代谢物且选择移除，则跳过
            if remove_common_metabolites and product.id.split('_')[0] in COMMON_METABOLITES:
                continue
                
            G.add_node(product.id, type='metabolite')
            node_types[product.id] = 'metabolite'
            G.add_edge(reaction.id, product.id)
    
    return G, node_types

def create_interactive_visualization(model, output_file):
    """
    创建交互式的代谢网络可视化
    
    参数:
    - model: cobra模型
    - output_file: 输出文件路径
    """
    # 创建网络图
    G, node_types = create_metabolic_network(model)
    
    # 使用Kamada-Kawai布局
    print("计算节点布局...")
    pos = nx.kamada_kawai_layout(G)
    
    # 创建节点轨迹
    node_x = []
    node_y = []
    node_text = []
    node_color = []
    
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(node.split('_')[0])  # 简化标签
        if node_types[node] == 'reaction':
            node_color.append('lightblue')
        else:
            node_color.append('lightgreen')
    
    # 创建边轨迹
    edge_x = []
    edge_y = []
    
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    # 创建节点轨迹
    nodes_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=node_text,
        textposition="top center",
        textfont=dict(size=8),
        marker=dict(
            color=node_color,
            size=10,
            line_width=2))
    
    # 创建边轨迹
    edges_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode='lines',
        line=dict(width=0.5, color='#888'),
        hoverinfo='none'
    )
    
    # 创建图形
    fig = go.Figure(data=[edges_trace, nodes_trace],
                   layout=go.Layout(
                       title='代谢网络图 (可交互)',
                       titlefont_size=16,
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20,l=5,r=5,t=40),
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                   )
    
    # 保存为HTML文件
    print(f"正在保存交互式网络图...")
    fig.write_html(output_file)
    
    print(f"网络图已保存至: {output_file}")
    print(f"节点总数: {len(G.nodes())}")
    print(f"边总数: {len(G.edges())}")

def main():
    try:
        # 加载模型
        print("正在加载模型...")
        model = load_json_model("integrated_model.json")
        
        # 创建输出目录
        output_dir = Path('output')
        output_dir.mkdir(exist_ok=True)
        
        # 生成交互式网络图
        print("\n生成交互式代谢网络图...")
        create_interactive_visualization(model, output_dir / 'metabolic_network.html')
        
    except Exception as e:
        print(f"\n程序执行出错: {str(e)}")
        raise

if __name__ == "__main__":
    main() 