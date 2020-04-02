import pickle
import time
from math import e
import networkx as nx
import community
import matplotlib.pyplot as plt
from matplotlib import pylab
import numpy as np
import pandas as pd


def load_graph():
    # get graph and race_dict
    graph_file = 'unsupervised/geno/pkl/graph'
    race_dict_file = 'unsupervised/geno/pkl/races_dict'
    with open(graph_file, 'rb') as handle:
        graph = pickle.load(handle)
    with open(race_dict_file, 'rb') as handle:
        race_dict = pickle.load(handle)
    return graph, race_dict


def plot_graph(graph, race_dict):
    print(graph.number_of_nodes())
    # print(graph.nodes[:5])
    print(graph.number_of_edges())
    # nx.draw(graph)
    nx.draw_spring(graph)

    # plt.savefig("graph.png")  # save as png
    plt.show()
    pass


def save_graph(graph,file_name):
    # initialze Figure
    plt.figure(num=None, figsize=(20, 20), dpi=80)
    plt.axis('off')
    fig = plt.figure(1)
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos)
    nx.draw_networkx_edges(graph, pos)
    nx.draw_networkx_labels(graph, pos)

    cut = 1.00
    xmax = cut * max(xx for xx, yy in pos.values())
    ymax = cut * max(yy for xx, yy in pos.values())
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)

    plt.savefig(file_name, bbox_inches="tight")
    pylab.close()
    del fig


if __name__ == '__main__':
    graph, race_dict = load_graph()
    graph = nx.convert_node_labels_to_integers(graph)
    subgraph = graph.subgraph(range(10000))
    plot_graph(subgraph, race_dict)
    # save_graph(subgraph, 'graph.png')
