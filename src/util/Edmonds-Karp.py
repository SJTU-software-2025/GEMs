from collections import deque

def bfs(residual_graph, source, sink, parent):
    num_nodes = len(residual_graph)
    visited = [False] * num_nodes
    queue = deque()
    
    queue.append(source)
    visited[source] = True
    parent[source] = -1  # 源节点无前驱
    
    while queue:
        u = queue.popleft()
        
        for v in range(num_nodes):
            if not visited[v] and residual_graph[u][v] > 0:
                queue.append(v)
                visited[v] = True
                parent[v] = u
                if v == sink:
                    return True  # 找到汇点，存在增广路径
    return False  # 未找到增广路径

def edmonds_karp(graph, source, sink):
    num_nodes = len(graph)
    # 初始化残留网络（直接复制输入图）
    residual_graph = [row[:] for row in graph]
    # parent数组用于记录增广路径中每个节点的前驱节点
    # 初始化为-1表示所有节点都还未访问，源节点的前驱将保持为-1
    parent = [-1] * num_nodes
    max_flow = 0
    
    # 当存在从源到汇的增广路径时
    while bfs(residual_graph, source, sink, parent):
        # 找到增广路径上的最小残留容量
        path_flow = float('inf')
        v = sink
        while v != source:
            u = parent[v]
            path_flow = min(path_flow, residual_graph[u][v])
            v = u
        
        # 更新残留网络和反向边
        v = sink
        while v != source:
            u = parent[v]
            residual_graph[u][v] -= path_flow  # 减少正向边容量
            residual_graph[v][u] += path_flow  # 增加反向边容量
            v = u
        
        max_flow += path_flow
        parent = [-1] * num_nodes  # 重置父节点数组
    
    return max_flow, residual_graph  # 返回流量和残留网络

def get_edge_flows(original_graph, residual_graph):
    flows = {}
    for u in range(len(original_graph)):
        for v in range(len(original_graph[0])):
            if original_graph[u][v] > 0:  # 只处理原始存在的边
                flow = original_graph[u][v] - residual_graph[u][v]
                # if flow > 0:
                flows[(u, v)] = flow
    return flows

# 示例用法
if __name__ == "__main__":
    # 邻接矩阵表示图（有向图）
    graph = [
        [0, 3, 2, 0],  # 节点0（源）
        [0, 0, 5, 2],  # 节点1
        [0, 0, 0, 3],  # 节点2
        [0, 0, 0, 0]   # 节点3（汇）
    ]
    source = 0
    sink = 3
    
    max_flow, residual = edmonds_karp(graph, source, sink)
    flows = get_edge_flows(graph, residual)
    
    print(f"最大流量: {max_flow}")
    print("边流量分配:")
    for (u, v), f in flows.items():
        print(f"边 {u}→{v}: {f}/{graph[u][v]}")