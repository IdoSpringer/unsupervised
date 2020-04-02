import pickle
import time
from math import e
import matplotlib
import networkx as nx
import community
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cProfile
# import matplotlib.patches as mpatches



def check_simlarity(geno1, geno2):
    count_similar = 0

    #geno1 = geno1.split('~')
    #geno2 = geno2.split('~')

    for i in range(0,10, 2):
        if geno1[i] == geno2[i]:
            if geno1[i+1] == geno2[i+1]:
                count_similar +=2
            else:
                count_similar +=1
        elif geno1[i+1] == geno2[i+1]:
            count_similar += 1
        elif  geno1[i] == geno2[i+1]:
            count_similar += 1
        elif  geno1[i+1] == geno2[i]:
            count_similar += 1

        if count_similar < i:
            return 0

    return count_similar

def create_graph(f_data, weights, resolution):
    t=time.time()
    list_9_nodes = []
    races_dict = {}
    graph = nx.Graph()
    with open(f_data) as f_in:
        for line in f_in:
            id, geno, race1, race2, _, _ = line.strip().split(',')

            if not geno in graph.nodes():
                graph.add_node(geno, ids={}, alleles = 10)
                alleles_list = geno.split('~')
                for i in range(len(alleles_list)):
                    alleles_9_list = alleles_list[0:i] + alleles_list[i+1:10]
                    geno_9 = ('~').join(alleles_9_list)
                    if geno_9 in graph.nodes():
                        adjs = graph.adj[geno_9]
                        for adj in adjs:
                            if geno != adj:
                                graph.add_edge(geno, adj, weight=weights[1])
                        graph.add_edge(geno, geno_9, weight=9)
                    else:
                        graph.add_node(geno_9,  alleles=9)
                        list_9_nodes.append(geno_9)
                        graph.add_edge(geno, geno_9, weight=9)

            graph.nodes[geno]['ids'][id] = (race1, race2)

            #graph.nodes[geno]['alleles'] = geno.split('~')

            races_dict[race1] = races_dict.get(race1, 0) + 1
            races_dict[race2] = races_dict.get(race2, 0) + 1

    print(len(races_dict))
    dict_count = {}
    print(len(graph.nodes()))

    for node in list_9_nodes:
        #if graph.node[node]['alleles'] == 9:
            graph.remove_node(node)

    #dict_genos = graph.nodes()
    """t = time.time()
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
                   start = True"""

    print(dict_count)
    pickle.dump(graph, open("./grid_pkl/graph_"+str(resolution)+"_"+str(weights[0]), "wb"))
    pickle.dump(races_dict, open("./grid_pkl/races_dict_"+str(resolution)+"_"+str(weights[0]), "wb"))
    #pickle.dump(dict_ids, open("./grid_pkl/dict_ids_"+str(resolution)+"_"+str(weights[0]), "wb"))
    print("graph creation time: ", time.time()-t)
    return graph, races_dict#, dict_ids


def community_detection(graph, reso, weights):
    t=time.time()
    partition = community.best_partition(graph, resolution=reso)
    pickle.dump(partition, open("./grid_pkl/cd_"+str(reso)+"_"+str(weights[0]), "wb"))
    print("communiry detection time:", time.time()-t)
    return partition


def community_nodes(cd, resolution, weights):
    t=time.time()
    com_nodes = dict()
    for node, com in cd.items():
        com_nodes[com] = com_nodes.get(com, set()) | set([node])
    pickle.dump(com_nodes,open("./grid_pkl/com_nodes_"+str(resolution)+"_"+str(weights[0]),"wb"))
    print("community nodes time:", time.time()-t)
    return com_nodes

def com_races(com_nodes, dict_genos):
    t = time.time()
    com_race_count = {}
    for com, genos in com_nodes.items():
        com_race_count[com] = dict()
        for geno in genos:
            for id, races in dict_genos[geno]['ids'].items():
                concat = races[0] + '_' + races[1]
                #concat =  dict_ids[id]['races'][0]+"_"+ dict_ids[id]['races'][1]
                com_race_count[com][concat] = com_race_count[com].get(concat, 0) + 1
    print("com races time:", time.time() - t)
    return com_race_count


def most_pop(com_race_count,f_results,weights,resolution):
    t=time.time()
    count_maxs=0
    count_all=0
    top_race = dict()
    percent_top_race = dict()
    com_size = dict()
    com_tags=dict()
    for com, races in com_race_count.items():
        dict_apperances = {}
        for race_count in races:
            race1, race2 = race_count.split('_')
            dict_apperances[race1] = dict_apperances.get(race1, 0) + races[race_count]*0.5
            dict_apperances[race2] = dict_apperances.get(race2, 0) + races[race_count] * 0.5
        most_common = sorted(dict_apperances.items(), key=lambda x: x[1], reverse=True)
        com_tags[com]=most_common
        top_race[com] = most_common[0][0]
        """if len(most_common) > 1:
            if most_common[0][1] == most_common[1][1]:
                top_race[com] = ('_').join(sorted[ most_common[0][0],  most_common[1][0]])"""
        com_size[com]=sum(races.values())
        percent_top_race[com]=most_common[0][1]/com_size[com]
        if com_size[com]>2:
            count_maxs+=most_common[0][1]
            count_all+=com_size[com]
    if count_all>0:
        total = count_maxs/count_all
    else:
        total=0
    print("total:",total," | max:", count_maxs, " | sum all:", count_all," | weights:", weights[0]," | resolution:", resolution)
    f_results.write(str(resolution) + '\t' + str(weights[0]) + '\t' + str(total*100) + '\n')
    print("most pop", time.time()-t)
    return top_race, com_size, percent_top_race, com_tags

def plot_pop(percent_top_race,com_size, resolution, weights):
    t= time.time()

    size, accuracy = [], []
    for com in percent_top_race:
        # plt.scatter(com_size[com], percent_top_race[com], color='blue', marker='*')
        size.append(com_size[com])
        accuracy.append(percent_top_race[com])

    plt.figure()
    plt.scatter(size, accuracy, color='blue', marker='*')
    plt.xlabel("Community Size")
    plt.tight_layout()
    # plt.savefig("./grid_plots/most pop race percent"+str(resolution)+"_"+str(weights[0])+".png")
    plt.show()
    print("plot pop time:", time.time()-t)
    return

def update_graph(weights):
    graph = pickle.load(open("./graph", "rb"))
    # for u,v in graph.edges:
    for u, v, weight in graph.edges.data('weight'):
        if graph[u][v]['weight']==8:
            graph[u][v]['weight']=weights[2]
        elif graph[u][v]['weight']==9:
            graph[u][v]['weight'] = weights[1]
        elif graph[u][v]['weight']==10:
            graph[u][v]['weight'] = weights[0]
    return graph

def plot_communitues(top_race, com_size, resolution, weights):
    t= time.time()
    bins = np.logspace(0, 16, num=17, base=2.0)
    com_bin = {}
    paint_list = []
    bin_list, com_per, all, communityid = [], [], [],[]

    for com, size in com_size.items():
        if com in top_race:
            for i in range(0, 16):
                b1 = bins[i]
                b2 = bins[i + 1]
                if size >= b1 and size < b2:
                    if size>0:
                        com_bin[b2] = com_bin.get(b2, 0)+1
                        paint_list.append(top_race[com])
                        communityid.append(str(com))
                        bin_list.append(b2)

    map = matplotlib.colors.ListedColormap(['black', 'red', 'tan', 'deepskyblue', 'blue', 'springgreen', 'pink', 'yellow', 'silver',
                   'darkorange', 'cyan', 'deeppink', 'salmon', 'teal', 'mediumpurple','gold','orchid','darkgreen','green','orange','navy','brown','rebeccapurple',
                   'darkviolet','lavender','grey','cornflowerblue','aquamarine','lightsteelblue'])

    partial = []
    for bin in bin_list:
        partial.append(1 / com_bin[bin])

    all.append(bin_list)
    all.append(partial)  # paint of com
    all.append(communityid)  # unique id
    all.append(paint_list)
    set_paint=set(paint_list)
    print("len coms painted", set_paint)

    rows = zip(all[0], all[1], all[2], all[3])
    headers = ['box', 'Value', 'id','paint_com']

    df = pd.DataFrame(rows, columns=headers)
    df = df.sort_values(by=['box', 'paint_com'])
    df.pivot_table(index='box', columns='paint_com', values='Value',
                   aggfunc='sum', fill_value=0).plot.bar(stacked=True, legend=True, colormap=map)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),fancybox=True, shadow=True,prop={'size': 6})
    plt.xlabel('Community Size', fontsize='large')
    plt.tight_layout()
    plt.savefig("./grid_plots/color com" + str(resolution ) + "_" + str(weights[0] ) + ".png")
    print("plot communities time: ", time.time()-t)
    # plt.show()
    e=0

def confusion_matrix(top_race,nodes_data,cd, races,weights,resolution):
    t= time.time()
    races_list= list(races.keys())
    matrix = np.zeros((len(races),len(races)))
    #for id in nodes_data:
    for geno in nodes_data:
        # if com_size[cd[geno]] >2:
        pred = top_race[cd[geno]]
        for id, races_id in nodes_data[geno]['ids'].items():
            real_race1, real_race2 = races_id[0],races_id[1]
            matrix[races_list.index(real_race1),races_list.index(pred)]+=0.5
            matrix[races_list.index(real_race2),races_list.index(pred)] += 0.5

    row_sums = matrix.sum(axis=1, dtype='float')
    new_matrix = np.zeros((len(races),len(races)))
    for i, (row, row_sum) in enumerate(zip(matrix, row_sums)):
        if row_sum==0:
            new_matrix[i, :]=0
        else:
            new_matrix[i, :] = row / row_sum

    new_matrix = np.around(new_matrix, 3)
    plt.matshow(new_matrix);
    plt.colorbar()
    plt.tight_layout()
    # plt.savefig("./grid_plots/cofution_matrix_" + str(resolution) + "_" + str(weights[0]) + ".png")
    plt.show()
    print("confusion matrix time:", time.time()-t)
    return

def ent_calculation(freq_d, base=None):
    base = 2 if base is None else base
    value, counts = freq_d.keys(), freq_d.values()
    counts = np.array(list(counts))
    # value, counts = np.unique(freq_l, return_counts=True)
    norm_counts = counts / sum(counts)
    return -(norm_counts * np.log(norm_counts) / np.log(base)).sum()


def check_entropy(com_size, com_labels):
    t= time.time()
    com_entropy_dict =dict()
    size, ent = [], []
    for com in com_labels:
        if com_size[com]>5:
            freqs = {}
            for item in com_labels[com]:
                for tag in item:
                    freqs[tag] = freqs.get(tag, 0) + item[1]
            x = ent_calculation(freqs)
            # x = ent_calculation(com_labels[com])
            com_entropy_dict[com] = x
    for com in com_entropy_dict:
        size.append(com_size[com])
        ent.append(com_entropy_dict[com])
    plt.figure()
    plt.scatter(size, ent, color='blue', marker='*', label='Community Entropy')
    plt.title('Entropy')
    plt.xlabel("Community Size")
    plt.ylabel("Entropy")
    plt.legend()
    plt.tight_layout()
    plt.savefig("./grid_plots/entropy" + str(resolution) + "_" + str(weights[0]) + ".png")
    print("entropy calculation time:", time.time()-t)
    # plt.show()
    return

if __name__ == "__main__":
    pr = cProfile.Profile()
    pr.enable()
    t=time.time()
    f_data = 'new_data_geno.txt'


    # weights_list = [[4,1,0.3],[2,1,0.5]]
    weights_list = [[10,1, 0.2]]
    resolution_list = [1]
    # resolution_list = np.logspace(0.001, 1, num=5, base=2)
    f_results = open('results.csv','w')
    for weights in weights_list:
        for resolution in resolution_list:

            # graph, race_dict = create_graph(f_data,weights, resolution)
            graph = pickle.load(open("./graph", "rb"))
            race_dict = pickle.load(open("./races_dict","rb"))
            print(len(race_dict))
            # graph = update_graph(weights)

            num_of_comp = nx.number_connected_components(graph)
            print("number of comp:", num_of_comp)
            # print(race_dict)

            # cd = community_detection(graph, resolution, weights)
            cd = pickle.load(open("./cd", "rb"))

            # com_nodes = community_nodes(cd,resolution, weights)
            com_nodes = pickle.load(open("./com_nodes", "rb"))
            print("number of communities:", len(com_nodes))
            com_race_count = com_races(com_nodes, graph.nodes())
            print(com_race_count)

            top_race, com_size, percent_top_race, com_labels = most_pop(com_race_count,f_results,weights,resolution)
            values=set(top_race.values())
            print("num of diffrent com tags", len(values))

            plot_pop(percent_top_race, com_size, resolution, weights)
            plot_communitues(top_race, com_size, resolution, weights)
            confusion_matrix(top_race, graph.nodes, cd, race_dict,weights, resolution)
            check_entropy(com_size, com_labels)
    print("total time:", time.time()-t)
    pr.disable()
    pr.print_stats(sort="time")

