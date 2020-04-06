import pickle
import time
import matplotlib
import networkx as nx
import community
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cProfile
import json
import pathlib



#create graph: nodes of genotypes, edges- between 2 genotypes with 9 identical alleles
def create_graph(f_data, weights, resolution):
    t=time.time()
    list_9_nodes = []
    races_dict = {}
    graph = nx.Graph()
    with open(f_data) as f_in:
        for line in f_in:
            id, geno, race1, race2 = line.strip().split(',')[0:4]
            geno = geno.replace('+', '~').replace('^', '~')

            #if id with a genotype that is not yet in the graph
            if not geno in graph.nodes():
                graph.add_node(geno, ids={}, alleles = 10)
                alleles_list = geno.split('~')
                #create temporal nodes of all combination of 9-alleles
                for i in range(len(alleles_list)):
                    alleles_9_list = alleles_list[0:i] + alleles_list[i+1:10]
                    geno_9 = ('~').join(alleles_9_list)
                    #if 9-alleles already in the graph - then have a genotype-1 match
                    if geno_9 in graph.nodes():
                        adjs = graph.adj[geno_9]
                        for adj in adjs:
                            if geno != adj:
                                graph.add_edge(geno, adj, weight=weights[0])
                        graph.add_edge(geno, geno_9, weight=9)
                    else:
                        graph.add_node(geno_9,  alleles=9)
                        list_9_nodes.append(geno_9)
                        graph.add_edge(geno, geno_9, weight=9)

            graph.nodes[geno]['ids'][id] = (race1, race2)

            #count how much from each race
            races_dict[race1] = races_dict.get(race1, 0) + 1
            races_dict[race2] = races_dict.get(race2, 0) + 1

    print(len(races_dict))
    dict_count = {}
    print(len(graph.nodes()))

    #remove temporal nodes
    for node in list_9_nodes:
        graph.remove_node(node)

    print(dict_count)
    #pickle.dump(graph, open("./grid_pkl/graph_"+str(resolution)+"_"+str(weights[0]), "wb"))
    #pickle.dump(races_dict, open("./grid_pkl/races_dict_"+str(resolution)+"_"+str(weights[0]), "wb"))
    print("graph creation time: ", time.time()-t)
    return graph, races_dict


#run community detection on graph
def community_detection(graph, reso, weights):
    t=time.time()
    partition = community.best_partition(graph, resolution=reso)
    pickle.dump(partition, open("./grid_pkl/cd_"+str(reso)+"_"+str(weights[0]), "wb"))
    print("communiry detection time:", time.time()-t)
    return partition


#convert community format- save to dict of community_number:list of genos
def community_nodes(cd, resolution, weights):
    t=time.time()
    com_nodes = dict()
    for node, com in cd.items():
        com_nodes[com] = com_nodes.get(com, set()) | set([node])
    pickle.dump(com_nodes,open("./grid_pkl/com_nodes_"+str(resolution)+"_"+str(weights[0]),"wb"))
    print("community nodes time:", time.time()-t)
    return com_nodes

#which races in each community, and number of appearances
def com_races(com_nodes, dict_genos):
    t = time.time()
    com_race_count = {}
    for com, genos in com_nodes.items():
        com_race_count[com] = dict()
        for geno in genos:
            for id, races in dict_genos[geno]['ids'].items():
                concat = races[0] + '_' + races[1]
                com_race_count[com][concat] = com_race_count[com].get(concat, 0) + 1
    print("com races time:", time.time() - t)
    return com_race_count

#for each community that contains more then 2 ids find:
# max: the number of ids from the race that appears the most
# count_all: number of all ids that in the community
# total: max/count_all
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


#option to change the wight of graph that read from pickle
def update_graph(weights):
    graph = pickle.load(open("./graph", "rb"))
    for u, v, weight in graph.edges.data('weight'):
        if graph[u][v]['weight']==8:
            graph[u][v]['weight']=weights[2]
        elif graph[u][v]['weight']==9:
            graph[u][v]['weight'] = weights[1]
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
    plt.xlabel("Predicted race")
    plt.ylabel("Real race")
    plt.savefig("./grid_plots/cofution_matrix_" + str(resolution) + "_" + str(weights[0]) + ".png")
    #plt.show()
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
    #f_data = 'new_data_geno.txt'

    pathlib.Path('grid_pkl').mkdir(parents=False, exist_ok=True)
    pathlib.Path('grid_plots').mkdir(parents=False, exist_ok=True)

    # Read configuration file and load properties
    with open('params_file.json') as f:
        params = json.load(f)

    config = {
        "f_data": params.get("data_file", 'new_data_geno.txt'),
        "weights": params.get("weights", [1]),
        "reso": params.get("resolution", 1),
    }

    weights_list = [config["weights"]]
    resolution_list = [config["reso"]]
    f_results = open('results.csv','w')
    for weights in weights_list:
        for resolution in resolution_list:

            graph, race_dict = create_graph(config["f_data"],weights, resolution)
            #graph = pickle.load(open("./grid_pkl/graph_"+str(resolution)+"_"+str(weights[0]), "rb"))
            #race_dict = pickle.load(open("./grid_pkl/races_dict_"+str(resolution)+"_"+str(weights[0]), "rb"))
            print(len(graph.nodes()))
            print(len(race_dict))
            # graph = update_graph(weights)

            num_of_comp = nx.number_connected_components(graph)
            print("number of comp:", num_of_comp)
            # print(race_dict)

            cd = community_detection(graph, resolution, weights)
            #cd = pickle.load(open("./grid_pkl/cd_"+str(resolution)+"_"+str(weights[0]), "rb"))

            com_nodes = community_nodes(cd,resolution, weights)
            #com_nodes = pickle.load(open("./grid_pkl/com_nodes_"+str(resolution)+"_"+str(weights[0]),"rb"))
            print("number of communities:", len(com_nodes))
            com_race_count = com_races(com_nodes, graph.nodes())
            print(com_race_count)

            top_race, com_size, percent_top_race, com_labels = most_pop(com_race_count,f_results,weights,resolution)
            values=set(top_race.values())
            print("num of diffrent com tags", len(values))

            plot_communitues(top_race, com_size, resolution, weights)
            confusion_matrix(top_race, graph.nodes, cd, race_dict,weights, resolution)
            check_entropy(com_size, com_labels)
    print("total time:", time.time()-t)
    pr.disable()
    #pr.print_stats(sort="time")

