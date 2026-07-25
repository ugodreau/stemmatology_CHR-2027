import requests, json, time
import os
import numpy as np
import matplotlib.pyplot as plt
import random as rnd
from tqdm import tqdm
import re
import pandas as pd
import matplotlib.pyplot as plt
from lxml import etree as et
from collections import Counter
import birth_death_utils as bd
import artificial_collations as ac
import output_utils as u
import pickle
import warnings
import networkx as nx
warnings.filterwarnings('ignore')



output_dir = os.getcwd()+'/beast_output/'
path = 'simulations_data/simulated_traditions_24_07'

for k in tqdm(range(100)):
    with open(f'{path}/tree_{k}.pickle', 'rb') as f:
        st = pickle.load(f)
    with open(f'{path}/sequences_{k}.pickle', 'rb') as f:
        wits_seqs = pickle.load(f)
    dates = {n:st.nodes[n]['birth_time'] for n in st.nodes() if st.nodes[n]['state']}

    model = u.BeastModel(seqs=wits_seqs,dates=dates,seq_name=f'sequence_{k}', model='fbd',n_iter=4*10**6)
    model.generate_xml().write(f'{path}/fdb_beast_{k}.xml')

    os.system(f'draft/BEAST2/beast/bin/beast -overwrite -prefix {output_dir} {path}/fdb_beast_{k}.xml')