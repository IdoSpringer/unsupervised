import pickle
import time
import networkx as nx


def check_similarity(geno1, geno2):
    count_similar = 0
    # geno1 = geno1.split('~')
    # geno2 = geno2.split('~')
    for i in range(0, 10, 2):
        if geno1[i] == geno2[i]:
            if geno1[i + 1] == geno2[i + 1]:
                count_similar += 2
            else:
                count_similar += 1
        elif geno1[i + 1] == geno2[i + 1]:
            count_similar += 1
        elif geno1[i] == geno2[i + 1]:
            count_similar += 1
        elif geno1[i + 1] == geno2[i]:
            count_similar += 1
        if count_similar < i:
            return 0
    return count_similar


def create_graph(f_data, weights, resolution):
    t = time.time()
    list_9_nodes = []
    races_dict = {}
    graph = nx.Graph()
    with open(f_data) as f_in:
        for line in f_in:
            id, geno, race1, race2, _, _ = line.strip().split(',')
            if geno not in graph.nodes():
                graph.add_node(geno, ids={}, alleles=10)
                alleles_list = geno.split('~')
                for i in range(len(alleles_list)):
                    alleles_9_list = alleles_list[0:i] + alleles_list[i+1:10]
                    geno_9 = '~'.join(alleles_9_list)
                    if geno_9 in graph.nodes():
                        adjs = graph.adj[geno_9]
                        for adj in adjs:
                            if geno != adj:
                                graph.add_edge(geno, adj, weight=weights[1])
                        graph.add_edge(geno, geno_9, weight=9)
                    else:
                        graph.add_node(geno_9, alleles=9)
                        list_9_nodes.append(geno_9)
                        graph.add_edge(geno, geno_9, weight=9)
            graph.nodes[geno]['ids'][id] = (race1, race2)
            # graph.nodes[geno]['alleles'] = geno.split('~')
            races_dict[race1] = races_dict.get(race1, 0) + 1
            races_dict[race2] = races_dict.get(race2, 0) + 1
    print(len(races_dict))
    dict_count = {}
    print(len(graph.nodes()))
    for node in list_9_nodes:
        # if graph.node[node]['alleles'] == 9:
            graph.remove_node(node)
    # dict_genos = graph.nodes()
    """
    t = time.time()
    for idx, key1 in enumerate(graph.nodes()):
        #val1 =dict_genos[key1]['alleles']
        val1 = dict_genos[key1]
        if (idx % 10000) == 0:
            print("index: ", idx)
            print("total time:", time.time() - t)


    # for key1 in graph.nodes():
        last_key = key1
        start = False
        for key2 in graph.nodes():
            if start:
                #similarity = check_simlarity(val1,dict_genos[key2]['alleles'])
                similarity = check_simlarity(val1, dict_genos[key2])
                if similarity == 8:
                    graph.add_edge(key1, key2, weight=weights[2])
                    dict_count[similarity] = dict_count.get(similarity, 0) + 1
                elif similarity == 9:
                    graph.add_edge(key1, key2, weight=weights[1])
                    dict_count[similarity] = dict_count.get(similarity, 0) + 1
            else:
               if key2 == last_key:
                   start = True
    """
    print(dict_count)
    pickle.dump(graph, open("./grid_pkl/graph_"+str(resolution)+"_"+str(weights[0]), "wb"))
    pickle.dump(races_dict, open("./grid_pkl/races_dict_"+str(resolution)+"_"+str(weights[0]), "wb"))
    # pickle.dump(dict_ids, open("./grid_pkl/dict_ids_"+str(resolution)+"_"+str(weights[0]), "wb"))
    print("graph creation time: ", time.time()-t)
    return graph, races_dict  # , dict_ids


def update_graph(weights):
    graph = pickle.load(open("./graph", "rb"))
    # for u,v in graph.edges:
    for u, v, weight in graph.edges.data('weight'):
        if graph[u][v]['weight'] == 8:
            graph[u][v]['weight']=weights[2]
        elif graph[u][v]['weight'] == 9:
            graph[u][v]['weight'] = weights[1]
        elif graph[u][v]['weight'] == 10:
            graph[u][v]['weight'] = weights[0]
    return graph


def load_graph(graph_file, race_dict_file):
    # get graph and race_dict
    # graph_file = 'unsupervised/geno/pkl/graph'
    # race_dict_file = 'unsupervised/geno/pkl/races_dict'
    with open(graph_file, 'rb') as handle:
        graph = pickle.load(handle)
    with open(race_dict_file, 'rb') as handle:
        race_dict = pickle.load(handle)
    return graph, race_dict


def subgraph(graph, number_of_nodes):
    graph = nx.convert_node_labels_to_integers(graph)
    subgraph = graph.subgraph(range(10000))
    return subgraph
