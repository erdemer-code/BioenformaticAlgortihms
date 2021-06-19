import random


def string_composition(text, k):
    kmers = [text[i:i + k] for i in range(len(text) - k + 1)]
    kmers.sort()
    return kmers


def prefix_suffix_check(kmer1, kmer2):
    window = 0

    for i in range(0, len(kmer1) - window + 1):
        for j in range(0, len(kmer2) - window + 1):
            if kmer1[i:i + window] == kmer2[j:j + window]:
                window += 1
    return window


def overlap_graph(kmers):
    kmers.pop(0)
    graph = {}
    for i in kmers:
        graph[i] = []

    for key1 in graph.keys():
        max_score = 0
        value = ""
        for key2 in kmers:
            score = 0
            if key1 != key2:
                score = prefix_suffix_check(key1, key2)
                if max_score < score:
                    max_score = score
                    value = key2

        graph[key1] = value
    g_list = []
    for g in graph:
        g_list.append(g)

    g_list.sort()

    for gl in g_list:
        print(f"{gl} --> {graph[gl]}")


def de_brujin_graph(kmers):
    res = {}
    edges = []
    nodes = set()
    k = len(kmers[0])
    for kmer in kmers:
        for i in range(len(kmer) - k + 1):
            edges.append((kmer[i:+k - 1], kmer[i + 1:i + k]))
            nodes.add(kmer[i:i + k - 1])
            nodes.add(kmer[i + 1:i + k])

    nodes = sorted(nodes)
    for n in nodes:
        res[n] = []
    for e in edges:
        res[e[0]].append(e[1])

    for r in res:
        print(f"{r} --> {', '.join(res[r])}")

#
# def check_degree(degrees):
#     for d in degrees.values():
#         if d != 0:
#             return False
#     return True
#
#
# def find_eulerian_cycle(graph):
#     degrees = {}
#     path = []
#     # random_start = 6
#     random_start = random.choice(list(graph.keys()))
#     if len(graph[random_start]) > 1:
#         path.append((random_start, graph[random_start][random.randint(0, len(graph[random_start]) - 1)]))
#     else:
#         path.append((random_start, graph[random_start][0]))
#     for g in graph.keys():
#         degrees[g] = 0
#     while True:
#         for _ in range(50):
#             next = path[-1][1]
#             if len(graph[next]) > 1:
#                 path.append((next, graph[next][random.randint(0, len(graph[next]) - 1)]))
#             else:
#                 path.append((next, graph[next][0]))
#             # path.append((next, random.choice(graph[next])))
#             path = list(dict.fromkeys(path))
#         for p in path:
#             degrees[p[0]] += 1
#             degrees[p[1]] -= 1
#         if check_degree(degrees):
#             break
#     # print(path)
#     res = ""
#     for i in range(0, len(path)):
#         if i < len(path) - 1:
#             res += str(path[i][0]) + "-->"
#         else:
#             res += str(path[i][0]) + "-->" + str(path[i][1])
#     return res, path
#
#
# def euler_cycle_iterator(graph):
#     res = None
#     path = None
#     res, path = find_eulerian_cycle(graph)
#     idx = [f[0] for f in path]
#     idx = list(dict.fromkeys(idx))
#     idx.sort()
#     full_path = set()
#     while idx != list(range(len(graph))):
#         res, path = find_eulerian_cycle(graph)
#         idx = [f[0] for f in path]
#         idx = list(dict.fromkeys(idx))
#         idx.sort()
#
#     return res
#

def eulerian_cycle(graph):
    edges = {}

    for i in range(len(graph)):
        edges[i] = len(graph[i])

    path = []
    res = []

    rand_start = random.randint(0, len(graph) - 1)
    path.append(rand_start)
    node = rand_start

    while len(path) > 0:

        if edges[node]:
            path.append(node)

            next = graph[node][-1]

            edges[node] -= 1
            graph[node].pop()
            node = next
        else:
            res.append(node)
            node = path[-1]
            path.pop()

    result_str = ""
    res.reverse()

    for i in range(len(res)):
        if i < len(res) - 1:
            result_str += str(res[i]) + " --> "
        else:
            result_str += str(res[i])

    return result_str


print("String Composition Problem: ")
print(string_composition("TATGGGGTGC", 3))
print("Overlap Graph Problem: ")
overlap_graph(["ATGCG", "GCATG", "CATGC", "AGGCA", "GGCAT"])
print("De Bruijn Graph from k-mers Problem: ")
de_brujin_graph(["GAGG", "CAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG"])
print("Eulerian Cycle Problem: ")
print(eulerian_cycle({
    0: [3],
    1: [0],
    2: [1, 6],
    3: [2],
    4: [2],
    5: [4],
    6: [5, 8],
    7: [9],
    8: [7],
    9: [6]}))
