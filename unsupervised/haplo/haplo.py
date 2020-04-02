import pickle
import time
from math import e

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import community
#import joblib

def check_simlarity(haps1, haps2):
    list_count_similar = [0, 0]

    for i in range(2):
        for j in range(2):
            hap1 = haps1[i].split('~')
            hap2 = haps2[j].split('~')
            count_similar = 0
            for locus in range(len(hap1)):
                if hap1[locus] == hap2[locus]:
                    count_similar += 1
            if count_similar >= 4:
                list_count_similar[(i+j)%2] += count_similar

    return max(list_count_similar)

def create_graph(f_data,weights, resolution):
    t = time.time()
    #dict_ids = {}
    races_dict = {}
    graph = nx.Graph()
    with open(f_data) as f_in:
        for line in f_in:
            h =  line.strip().split(',')
            id, hap_race1, hap_race2, _, _ = line.strip().split(',')
            hap1, race1 = hap_race1.split(';')
            hap2, race2 = hap_race2.split(';')
            #dict_ids[id] = [(hap1, hap2), (race1, race2)]

            haps =  (',').join(sorted([hap1, hap2]))
            if not haps in graph.nodes():
                graph.add_node(haps, ids={})
            graph.node[haps]['ids'][id] = (race1, race2)

            races_dict[race1] = races_dict.get(race1, 0) + 1
            races_dict[race2] = races_dict.get(race2, 0) + 1

    print(len(races_dict))
    dict_count = {}
    l = 0

    for key1 in graph.nodes():
        last_key = key1
        start = False
        for key2 in graph.nodes():
            if start:
                l+=1
                similarity = check_simlarity(key1.split(','), key2.split(','))
                if similarity == 4:
                    graph.add_edge(key1,key2 , weight = weights[3])
                if similarity == 5:
                    graph.add_edge(key1,key2 , weight = weights[2])
                if similarity == 9:
                    graph.add_edge(key1,key2 , weight = weights[1])
                if similarity == 10:
                    graph.add_edge(key1,key2 , weight = weights[0])

                dict_count[similarity] = dict_count.get(similarity, 0) + 1
            else:
               if key2 == last_key:
                   start = True
    print(dict_count)
    print(l)
    pickle.dump(graph, open("./grid_pkl/graph_" + str(resolution) + "_" + str(weights[0]), "wb"))
    pickle.dump(races_dict, open("./grid_pkl/races_dict_" + str(resolution) + "_" + str(weights[0]), "wb"))
    #joblib.dump(dict_ids, open("./grid_pkl/dict_ids_" + str(resolution*100) + "_" + str(weights[0]*100), "wb"))
    print("graph creation time: ", time.time() - t)
    return graph, races_dict#, dict_ids



def community_detection(graph,reso, weights):
    t=time.time()
    partition = community.best_partition(graph,resolution=reso)
    # pickle.dump(partition, open("./cd_haplo.pkl", "wb"))
    pickle.dump(partition, open("./grid_pkl/cd_" + str(reso) + "_" + str(weights[0]), "wb"))
    print("communiry detection time:", time.time() - t)
    return partition


def community_nodes(cd, resolution, weights):
    t = time.time()
    com_nodes = dict()
    for node, com in cd.items():
        com_nodes[com] = com_nodes.get(com, set()) | set([node])
    # pickle.dump(com_nodes,open("./com_nodes_haplo","wb"))
    pickle.dump(com_nodes, open("./grid_pkl/com_nodes_" + str(resolution) + "_" + str(weights[0]), "wb"))
    print("community nodes time:", time.time() - t)
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
    for com, races in com_race_count.items():
        dict_apperances = {}
        for race_count in races:
            race1, race2 = race_count.split('_')
            dict_apperances[race1] = dict_apperances.get(race1, 0) + races[race_count]*0.5
            dict_apperances[race2] = dict_apperances.get(race2, 0) + races[race_count] * 0.5
        most_common = sorted(dict_apperances.items(), key=lambda x: x[1], reverse=True)
        top_race[com] = most_common[0][0]
        """if len(most_common) > 1:
            if most_common[0][1] == most_common[1][1]:
                top_race[com] = ('_').join(sorted[ most_common[0][0],  most_common[1][0]])"""
        com_size[com]=sum(races.values())
        percent_top_race[com]=most_common[0][1]/com_size[com]
        if com_size[com]>5:
            count_maxs+=most_common[0][1]
            count_all+=com_size[com]
    if count_all>0:
        total = count_maxs/count_all
    else:
        total=0
    print("total:",total," | max:", count_maxs, " | sum all:", count_all," | weights:", weights," | resolution:", resolution)
    f_results.write(str(resolution) + '\t' + str(weights[0]) + '\t' + str(total*100) + '\n')
    print("most pop", time.time()-t)
    return top_race, com_size, percent_top_race

def plot_pop(percent_top_race,com_size, resolution, weights):
    plt.figure()
    for com in percent_top_race:
        plt.scatter(com_size[com], percent_top_race[com], color='blue', marker='*')
    plt.xlabel("Community Size")
    plt.savefig("./grid_plots/most pop race percent"+str(resolution)+"_"+str(weights[0])+".png")
    return

def update_graph(weights):
    graph=pickle.load(open("./graph_haplo_new","rb"))
    for u, v, weight in graph.edges.data('weight'):
        if graph[u][v]['weight']==0.2:
            graph[u][v]['weight']=weights[3]
        elif graph[u][v]['weight']==0.4:
            graph[u][v]['weight'] = weights[2]
        elif graph[u][v]['weight']==0.7:
            graph[u][v]['weight'] = weights[1]
        elif graph[u][v]['weight']==1:
            graph[u][v]['weight'] = weights[0]
            print("1")
    return graph

def plot_communitues(top_race, com_size, resolution, weights):
    bins = np.logspace(0, 12, num=13, base=2.0)
    com_bin = {}
    paint_list = []
    bin_list, com_per, all, communityid = [], [], [],[]

    for com, size in com_size.items():
        if com in top_race:
            for i in range(0, 12):
                b1 = bins[i]
                b2 = bins[i + 1]
                if size >= b1 and size < b2:
                    if size>5:
                        com_bin[b2] = com_bin.get(b2, 0)+1
                        paint_list.append(top_race[com])
                        communityid.append(str(com))
                        bin_list.append(b2)

    map = matplotlib.colors.ListedColormap(['black', 'red', 'tan', 'deepskyblue', 'blue', 'springgreen', 'pink', 'yellow', 'silver',
                   'darkorange', 'cyan', 'deeppink', 'salmon', 'teal', 'mediumpurple','gold','orchid','rebeccapurple','green','orange','navy','brown','darkgreen',
                   'lavender','crimson','lightgreen','grey','cornflowerblue','aquamarine','lightsteelblue'])

    partial = []
    for bin in bin_list:
        partial.append(1 / com_bin[bin])

    all.append(bin_list)
    all.append(partial)  # paint of com
    all.append(communityid)  # unique id
    all.append(paint_list)

    rows = zip(all[0], all[1], all[2], all[3])
    headers = ['box', 'Value', 'id','paint_com']

    df = pd.DataFrame(rows, columns=headers)
    df = df.sort_values(by=['box', 'paint_com'])
    df.pivot_table(index='box', columns='paint_com', values='Value',
                   aggfunc='sum', fill_value=0).plot.bar(stacked=True, legend=True, colormap=map)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),fancybox=True, shadow=True)
    plt.xlabel('Community Size', fontsize='large')
    plt.savefig("./grid_plots/color com" + str(resolution ) + "_" + str(weights[0] ) + ".png")



def confusion_matrix(top_race,nodes_data,cd, races,weights,resolution):
    races_list= list(races.keys())
    matrix = np.zeros((len(races),len(races)))
    for id in nodes_data:
        real_race1, real_race2 = nodes_data[id]['races'][0],nodes_data[id]['races'][1]
        pred = top_race[cd[id]]
        matrix[races_list.index(pred),races_list.index(real_race1)]+=0.5
        matrix[races_list.index(pred),races_list.index(real_race2)] += 0.5

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
    plt.savefig("./grid_plots/cofution_matrix_" + str(resolution) + "_" + str(weights[0]) + ".png")

def ent_calculation(freq_d, base=None):
    base = e if base is None else base
    value, counts = freq_d.keys(), freq_d.values()
    counts = np.array(list(counts))
    # value, counts = np.unique(freq_l, return_counts=True)
    norm_counts = counts / sum(counts)
    return -(norm_counts * np.log(norm_counts) / np.log(base)).sum()


def check_entropy(com_size, com_labels):
    com_entropy_dict =dict()
    size, ent = [], []
    for com in com_labels:
        if com_size[com]>5:
            freqs = {}
            for item in com_labels[com]:
                for tag in item:
                    freqs[tag] = freqs.get(tag, 0) + item[1]
            x = ent_calculation(freqs)
            com_entropy_dict[com]=x
    for com in com_entropy_dict:
        size.append(com_size[com])
        ent.append(com_entropy_dict[com])
    plt.figure()
    plt.scatter(size, ent, color='blue', marker='*', label='Community Entropy')
    plt.title('Entropy')
    plt.xlabel("Community Size")
    plt.ylabel("Entropy")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    t = time.time()
    f_data = 'don_pat_haplotypes.csv'
    weights = [12, 10, 4, 0.2]
    # resolution_list = [1, 0.75, 0., 0.25]
    resolution_list=[1]
    f_results = open('results.csv', 'w')
    for resolution in resolution_list:
        #f_data = 'input'
        # graph, race_dict = create_graph(f_data,weights,resolution)
        # graph=pickle.load(open("./graph_haplo_new","rb"))
        graph = update_graph(weights)
        num_of_comp = nx.number_connected_components(graph)
        print("number of comp:", num_of_comp)
        # race_dict=pickle.load(open("./pkl","rb"))
        # print(race_dict)
        cd = community_detection(graph,resolution, weights)
        # cd = pickle.load(open("./cd.pkl", "rb"))
        com_nodes = community_nodes(cd,resolution, weights)
        print('num clusters: ', len(com_nodes))
        com_race_count = com_races(com_nodes, graph.nodes())
        print(com_race_count)
        top_race, com_size, percent_top_race = most_pop(com_race_count, f_results, weights, resolution)
        plot_pop(percent_top_race, com_size, resolution, weights)
    print("total time:", time.time() - t)