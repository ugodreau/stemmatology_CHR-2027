import os
import functools
import itertools
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.traversal.depth_first_search import dfs_tree
from tqdm import tqdm


def leaves(graph):
    """
    Returns the terminal leaves of a tree

    Parameters
    ----------
    graph : nx.DiGraph()
        tree (directed acyclic graph)

    Returns
    -------
    list : list of node objects
        list of node labels of leaves

    """
    return [node for node in graph.nodes() if graph.out_degree(node) == 0]

def internal_nodes(graph):
    """
    Returns the non-leaves node of a tree

    Parameters
    ----------
    graph : nx.DiGraph()
        tree (directed acyclic graph)

    Returns
    -------
    list : list of node objects
        list of internal nodes of a graph

    """
    return [n for n in graph.nodes() if not n in leaves(graph)]

def root(graph):
    """
    Returns the root of a tree
    
    Parameters
    ----------
    graph : nx.DiGrap()
        tree (directed acyclic graph)

    Returns
    -------
    n : node object
        label of the root of graph

    """
    for n in graph.nodes():
        if graph.in_degree(n) == 0:
            return n


def subtree_size(graph, node):
    """
    Returns the size of the subtree of graph with given root node

    Parameters
    ----------
    graph : nx.DiGraph()
        tree
    node : root of the subtree

    Returns
    -------
    int
        Size of the subtree

    """
    return len(dfs_tree(graph,node).nodes())

def living_nb(g):
    return sum([g.nodes[n]['state'] for n in g.nodes()])

def generate_tree(lda, mu, Tact, Tinact):
    """
    Generate a tree (arbre réel) according to birth death model.

    Parameters
    ----------
    lda : float
        birth rate of new node per node per iteration
    mu : float
        death rate of nodes per node per per iteration
    Nact : int
        number of iterations of the active reproduction phase
    Ninact : int
        number of iterations of the pure death phase (lda is set to 0)

    Returns
    -------
    G : nx.DiGraph()
        networkx graph object of the generated tree with following node attributes:
            'state' : boolean, True if node living at the end of simulation
            'birth_time' : int
            'death_time' : int

    """
    currentID = 0
    G = nx.DiGraph()
    G.add_node(currentID)
    living_nodes = set([0])

    birth_time = {0:0}
    death_time = {}

    pop = 1
    prob_birth = lda / (lda + mu)
    prob_death =  mu / (lda + mu)
    prob_event = lda + mu

    t = 0

    while t < Tact:
        if pop == 0:
            t = Tact
            break
        next_event = np.random.exponential(scale = 1. / (prob_event * pop))
        if t + next_event > Tact:
            t = Tact
            break

        t += next_event
        r = np.random.rand()
        current_node = np.random.choice(list(living_nodes))
        if r < prob_birth:
            currentID += 1
            G.add_node(currentID)
            G.add_edge(current_node, currentID)
            living_nodes.add(currentID)
            pop += 1
            birth_time[currentID] = t
        else:
            living_nodes.remove(current_node)
            pop -= 1
            death_time[current_node] = t
    
    while t < Tact + Tinact:
        if pop == 0:
            t = Tact + Tinact
            break
        next_event = np.random.exponential(scale = 1. / (mu * pop))
        if t + next_event > Tact + Tinact:
            t = Tact + Tinact
            break
        t += next_event
        current_node = np.random.choice(list(living_nodes))
        living_nodes.remove(current_node)
        pop -= 1
        death_time[current_node] = t
    
    living = {n:(n in living_nodes) for n in G.nodes()}
    nx.set_node_attributes(G, living, 'state')
    nx.set_node_attributes(G, birth_time, 'birth_time')
    nx.set_node_attributes(G, death_time, 'death_time')

    return G

def generate_stemma(GG, assign_length=True):
    G = nx.DiGraph(GG)
    living = {n: G.nodes[n]['state'] for n in list(G.nodes())}
    
    # recursivelly remove dead leaves until all terminal nodes are living witnesses
    
    terminal_dead_nodes = [n for n in leaves(G) if not living[n]]
    while terminal_dead_nodes != []:
        for n in terminal_dead_nodes:
            G.remove_node(n)
        terminal_dead_nodes = [n for n in leaves(G) if not living[n]]

    # remove non-branching consecutive dead nodes

    if assign_length:    
        for e in G.edges():
            G.edges[e]['length'] = 0

    unwanted_virtual_nodes = [n for n in list(G.nodes()) if living[n] == False and G.out_degree(n) == 1 and G.in_degree(n) == 1]

    while unwanted_virtual_nodes != []:
        for n in unwanted_virtual_nodes:
            f = list(G.predecessors(n))[0]
            c = list(G.successors(n))[0]
            G.add_edge(f, c)
            G.edges[(f,c)]['length'] = 1 + G.edges[(f,n)]['length'] + G.edges[(n,c)]['length'] 
            G.remove_node(n)
        unwanted_virtual_nodes = [n for n in list(G.nodes()) if living[n] == False and G.in_degree(n) == 1 and G.out_degree(n) ==1]

    if not living[root(G)] and G.out_degree(root(G)) == 1:
        G.remove_node(root(G))
    
    for e in G.edges():
        G.edges[e]['label'] = str(G.edges[e]['length'])
        
    return G

def draw_tree(G, filename):
    """
    Draw a tree generated from generate_tree or generate_stemma as svg file 
    with living nodes colored in red and grey dead nodes

    Parameters
    ----------
    G : nx.DiGraph()
        tree object with at least 'state' attribute given
    filename : string
        

    Returns
    -------
    None.
    """
    living = nx.get_node_attributes(G, 'state')
    color_map = {node: 'red' if state else 'grey' for node, state in zip(living.keys(), living.values())}
    nx.set_node_attributes(G, color_map, name='color')
    nx.nx_pydot.write_dot(G, 'graph.dot')
    os.system('dot -Tsvg graph.dot > {}.svg'.format(filename))

def generate_tree_fixed_size(lda, mu, Ta,Ti, n):
    g = generate_tree(lda, mu, Ta,Ti)
    while living_nb(g) != n:
        g = generate_tree(lda, mu, Ta, Ti)
    return g