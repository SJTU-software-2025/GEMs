from collections import deque

def bfs(graph, start):
    """
    广度优先搜索算法实现
    :param graph: 图的邻接表表示（字典或列表）
    :param start: 起始节点
    :return: 按BFS顺序访问的节点列表
    """
    visited = set()          # 记录已访问节点
    queue = deque()          # 使用双端队列实现队列
    traversal_order = []     # 存储遍历顺序
    
    queue.append(start)      # 将起始节点加入队列
    visited.add(start)       # 标记为已访问
    
    while queue:
        node = queue.popleft()   # 从队列左侧取出节点（先进先出）
        traversal_order.append(node)
        
        # 遍历当前节点的所有邻居
        for neighbor in graph[node]:
            if neighbor not in visited:
                queue.append(neighbor)
                visited.add(neighbor)
    
    return traversal_order


# 示例用法
if __name__ == "__main__":
    # 图的邻接表表示（字典形式）
    graph = {
        'A': ['B', 'C'],
        'B': ['A', 'D', 'E'],
        'C': ['A', 'F'],
        'D': ['B'],
        'E': ['B', 'F'],
        'F': ['C', 'E']
    }

    print("BFS遍历顺序:", bfs(graph, 'A'))
    # 输出 : ['A', 'B', 'C', 'D', 'E', 'F']