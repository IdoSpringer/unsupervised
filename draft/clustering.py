import community
from sklearn.cluster import SpectralClustering
from draft.utils import *

# methods:
# spectral clustering
# louvain
# dbscan


def louvain_community_detection(graph, reso, weights):
    t = time.time()
    # louvain algorithm
    partition = community.best_partition(graph, resolution=reso)
    # save result to file
    pickle.dump(partition, open("./grid_pkl/cd_" + str(reso) + "_" + str(weights[0]), "wb"))
    print("community detection time:", time.time() - t)
    return partition


def load_louvain_communities():
    pass


def spectral_clustering(graph, number_of_clusters):
    adjacency_matrix = nx.adjacency_matrix(graph)
    sc = SpectralClustering(number_of_clusters, affinity='precomputed', n_init=100, assign_labels='discretize')
    sc.fit_predict(adjacency_matrix)
    pass


def dbscan():
    pass


if __name__ == '__main__':
    graph_file = 'geno/pkl/graph'
    race_dict_file = 'geno/pkl/races_dict'
    graph, race_dict = load_graph(graph_file, race_dict_file)
    spectral_clustering(graph, number_of_clusters=20)
    pass


# evaluation methods:
# using labels (color clusters by max)
# entropy based