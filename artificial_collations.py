from collections import Counter
import numpy as np
import networkx as nx
import birth_death_utils as bd
from tqdm.notebook import tqdm
import random as rnd

def new_sequence(source, inventory, p_new_locus, p_mod, p_inov):
    """Generate a new sequence corresponding to the text produced from a source
    after a single direct copy

    Parameters
    ----------
    source : dict
        sequence of the source manuscript, structured as
        {l1 : {app1.1: <reading>, app1.2: <reading>,...},
        l2 : {app2.1: <reading>,...},
        ...}
    inventory : dict
        dictionary enumerating all readings present in the current state
        of the tradition, for all existing loci in all manuscripts, structured as
        {l1 : {app1.1 :[rdg1, rdg2,...], app1.2 : [...], ...},
        ...}
    p_new_locus : float
        probability per line and unit time of creating a new variant locus,
        introducing a second reading wrt root reading
    p_mod : float
        probability of activating and modifying an already existing variant locus in the
        current state of the tradition per locus and unit time
    p_inov : float
        probability that the modification of activated locus (with rate p_mod)
        introduces a new reading (inovation), otherwise, new reading while be selcted
        among already existing readings at the locus
        

    Returns
    -------
    dict
        sequence of the new manuscript, structured as source
    """
    new_seq = {l:{app:0 for app in inventory[l].keys()} for l in inventory.keys()}
    for l in source.keys():
        for app in source[l].keys():
            new_seq[l][app] = source[l][app]

    
    # creation of a new locus with rate p_new_locus per lines
    for l in source.keys():
        r = rnd.random()
        if r < p_new_locus:
            new_app_index = max(inventory[l].keys()) + 1 if len(inventory[l]) > 0 else 0
            new_seq[l][new_app_index] = 1
            inventory[l][new_app_index] = [0,1]
    
    # mobilization of an already existing locus (potentially elsewhere in trad.)
    # with rate p_mod per all loci in trad. 
    for l in inventory.keys():
        for app in inventory[l].keys():
            
            r = rnd.random()
            if r < p_mod:
                rr = rnd.random()

                # creation of new variant with prob. p_inov
                if rr < p_inov:
                    new_var_value = max(inventory[l][app]) + 1 if len(inventory[l][app]) > 0 else 0
                    new_seq[l][app] = new_var_value
                    inventory[l][app].append(new_var_value)
                
                # reversion/convergence to another existing variant with prob. 1 - p_inov
                else:
                    target_value = rnd.choice([v for v in inventory[l][app] if not v==new_seq[l][app]])
                    new_seq[l][app] = target_value
    
    return new_seq


def generate_collation_from_tree(g, n_sites, p_new_locus, p_mod, p_inov):
    """
    Generate a full collation from a given tree g (nx.DiGraph object) representing
    the _arbre reel_ of the tradition.
    Returns the sequence of all witnesses (surviving manuscripts) and the inventory
    of all readings for all loci/app.
    """
    root = bd.root(g)
    g.nodes[root]['seq'] = {i:{} for i in range(n_sites)}

    inventory = {i:{} for i in range(n_sites)}

    to_assign = list(g.successors(root))

    while len(to_assign) > 0:
        for ms in to_assign:
            parent = list(g.predecessors(ms))[0]
            parent_seq = g.nodes[parent]['seq']
            child_sequence = new_sequence(parent_seq, inventory, p_new_locus, p_mod, p_inov)
            g.nodes[ms]['seq'] = child_sequence
        
        to_assign_new = []
        for ms in to_assign:
            to_assign_new += list(g.successors(ms))
        to_assign = list(to_assign_new)
    
    wits_sequences = {n:g.nodes[n]['seq'] for n in g.nodes() if g.nodes[n]['state']}

    return inventory, wits_sequences


def reduce_collation(inventory,wits_sequences):
    '''Reduce the sequences of witnesses from generate_collation_from_tree
    to the values present in the witnesses only.
    '''
    wits = wits_sequences.keys()
    reduced_sequences = {w:{l:{} for l in inventory.keys()} for w in wits}
    locus_dim = {l:{} for l in inventory.keys()}

    for l in inventory.keys():
        for app in inventory[l].keys():
            all_rdgs = set([seq[l][app] if app in seq[l] else 0 for seq in wits_sequences.values()])
            if len(all_rdgs) > 1:
                reduced_values = {v:i for v,i in zip(list(all_rdgs), range(len(all_rdgs)))}
                locus_dim[l][app] = len(reduced_values)
                for w in wits:
                    if app in wits_sequences[w][l]:
                        reduced_sequences[w][l][app] = reduced_values[wits_sequences[w][l][app]]
                    else:
                        reduced_sequences[w][l][app] = reduced_values[0]
    return locus_dim, reduced_sequences

def get_variants_distributions(p_app, p_mod, p_inov, n_sites, tree_set, n_runs):
    '''Returns the distributions of
        + number of loci/app per verses/l
        + number of values/rdg per loci/app
    for a set of simulated tradition of n_sites verses, averaged over n_trads
    realizations of the process and over all trees in tree_set
    '''
    counter_app = []
    counter_l = []
    for g in tree_set:
        dim_list = []
        for i in range(n_runs):
            inv, seqs_all_mss = generate_collation_from_tree(g, n_sites, p_app, p_mod, p_inov)
            dims, seqs = reduce_collation(inv,seqs_all_mss)
            app_nb_per_l = []
            for l in dims.values():
                app_nb_per_l.append(len(l.keys()))
                for app in l.values():
                    dim_list.append(app)
        
        c_app = Counter(dim_list)
        c_l = Counter(app_nb_per_l)
        counter_app.append(c_app)
        counter_l.append(c_l)

    maxval_app = max((sum([list(c.keys()) for c in counter_app],[])))
    maxval_l = max((sum([list(c.keys()) for c in counter_l],[])))

    X_app = np.arange(2,maxval_app+1)
    X_l = np.arange(0,maxval_l+1)

    Y_mean_app = np.array([])
    Y_mean_l = np.array([])
    
    Y_var_app = np.array([])
    Y_var_l = np.array([])

    for x_l in X_l:
        y_l = []
        for c_l in counter_l:
            y_l.append(c_l[x_l]/sum(list(c_l.values())))

        Y_mean_l = np.append(Y_mean_l, np.mean(y_l))
        Y_var_l = np.append(Y_var_l,np.std(y_l))
    
    for x_app in X_app:
        y_app = []
        for c_app in counter_app:
            y_app.append(c_app[x_app]/sum(list(c_app.values())))

        Y_mean_app = np.append(Y_mean_app, np.mean(y_app))
        Y_var_app = np.append(Y_var_app,np.std(y_app))
    
    return (X_l,Y_mean_l,Y_var_l), (X_app,Y_mean_app, Y_var_app)

def exp_extrapol(X, Y, iw, f_length):
    '''
    Extrapolate a finite probability distribution with best fitting exponential tail
    '''
    X1,Y1 = list(X), list(Y)
    xi = X[-iw:]
    yi = np.log(Y[-iw:])
    alpha, beta = np.polyfit(xi,yi,1)
    for i in range(X[-1]+1, f_length+1):
        X1.append(i)
        Y1.append(np.exp(alpha*i+beta))
    return np.array(X1), np.array(Y1)/sum(Y1)