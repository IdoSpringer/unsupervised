import pickle
import networkx as nx
import community

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


def create_graph(f_data):
    dict_ids = {}
    races_dict = {}
    graph = nx.Graph()
    with open(f_data) as f_in:
        for line in f_in:
            id, hap_race1, hap_race2, _, _ = line.strip().split(',')
            hap1, race1 = hap_race1.split(';')
            hap2, race2 = hap_race2.split(';')
            dict_ids[id] = [(hap1, hap2), (race1, race2)]

            graph.add_node(id)#, haplotypes = [hap1, hap2], races = [race1, race2])


            races_dict[race1] = races_dict.get(race1, 0) + 1
            races_dict[race2] = races_dict.get(race2, 0) + 1


    print(len(races_dict))



    for key1, val1 in dict_ids.items():
        for key2, val2 in dict_ids.items():
            if key1 != key2:
                similarity = check_simlarity(val1[0], val2[0])
                if similarity == 4:
                    graph.add_edge(key1,key2 , weight = 3)
                if similarity == 5:
                    graph.add_edge(key1,key2 , weight = 5)
                if similarity == 9:
                    graph.add_edge(key1,key2 , weight = 6)
                if similarity == 10:
                    graph.add_edge(key1,key2 , weight = 7)

    return graph, races_dict, dict_ids



def community_detection(graph):
    partition = community.best_partition(graph)
    pickle.dump(partition, open("./cd_haplo.pkl", "wb"))
    return partition


def community_nodes(cd):
    com_nodes = dict()
    for node, com in cd.items():
        com_nodes[com] = com_nodes.get(com, set()) | set([node])
    pickle.dump(com_nodes,open("./com_nodes_haplo","wb"))
    return com_nodes

def com_races(com_nodes, dict_ids):
    com_race_count = {}
    for com, ids in com_nodes.items():
        com_race_count[com] = dict()
        for id in ids:
            concat =  dict_ids[str(id)][1]+"_"+ dict_ids[str(id)][2]
            com_race_count[com][concat] = com_race_count[com].get(concat, 0) + 1

    return com_race_count


if __name__ == "__main__":
    f_data = 'don_pat.high.hap.freqs'
    #f_data = 'input'
    graph, race_dict, dict_ids = create_graph(f_data)
    print(race_dict)
    cd = community_detection(graph)
    com_nodes = community_nodes(cd)
    com_race_count = com_races(com_nodes, dict_ids)