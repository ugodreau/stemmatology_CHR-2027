import os
import birth_death_utils as bd
import networkx as nx
import numpy as np
from lxml import etree as et
import artificial_collations as ac

##### BEAST2 functions #######

MODELS = {'bdsky':'BDSKY_serial.xml'}
XPATH = {'substModel':"/beast/run[@id='mcmc']/distribution[@id='posterior']/distribution[@id ='likelihood']/distribution[@id = 'treeLikelihood.{name}']/siteModel[@id='SiteModel.s:{name}']/substModel[@id='LewisMK.s:{name}']",
    'mcmc': "/beast/run[@id='mcmc']",
    'dates': "run/state/tree/trait[@traitname='date']"
    }

class BeastModel:
    def __init__(self,
        sw = {}, 
        dates = {}, 
        seq_name = 'beast_sequence',
        model = 'bdsky', 
        params = {}, 
        n_iter=10000):

        self.sw = sw
        self.dates = dates
        self.seq_name = seq_name
        self.model = model
        self.params = params
        self.n_iter = n_iter
        
        self.n_sites = len(list(sw.values())[0])
        self.char_nb = self.compute_char_nb(self.sw)
    
    def compute_char_nb(self, sw):
        all_values = []
        for seq in sw.values():
            all_values.append(max(list(map(int, seq))))
        char_nb = max(all_values) + 1
        return char_nb

    def set_attribute(self, xml, element_name, attr, value):
        element_xpath = XPATH[element_name].format(name=self.seq_name)
        element = xml.xpath(element_xpath)[0]
        element.attrib[attr] = value
    
    def generate_xml(self):

        ### load beast config file template ###

        with open(f"beast_templates/{MODELS[self.model]}",'r') as f:
            content = f.read()
            content = content.replace("dummyname", self.seq_name)
        with open('_temp_prototype.xml','w') as f:
            f.write(content)
        with open('_temp_prototype.xml', 'r') as f:
            parser = et.XMLParser(remove_comments=True)
            file = et.parse(f, parser=parser)
        os.system("rm _temp_prototype.xml")

        ### set sequences ####

        for w,s in self.sw.items():
            new_sequence = et.SubElement(file.xpath('data')[0], 'sequence')
            new_sequence.attrib['id'] = f'seq_{w}'
            new_sequence.attrib['spec'] ="Sequence"
            new_sequence.attrib['taxon'] = w
            new_sequence.attrib['totalcount'] = str(self.char_nb)
            new_sequence.attrib['value'] = s

        user_data_type = et.SubElement(file.xpath('data')[0], 'userDataType')
        user_data_type.attrib['id'] = "StandardData.0"
        user_data_type.attrib['spec'] = "beast.base.evolution.datatype.StandardData"
        user_data_type.attrib['ambiguities'] = ""
        user_data_type.attrib['nrOfStates'] = str(self.char_nb)

        subst_model = file.xpath(XPATH['substModel'].format(name=self.seq_name))[0]
        subst_model.attrib['stateNumber'] = str(self.char_nb)


        ### set dates ###

        if self.dates != None:
            dates_str = ''
            for w,d in list(self.dates.items())[:-1]:
                dates_str +=f'{w}={d},'
            dates_str += '{}={}'.format(*list(self.dates.items())[-1])

            print(dates_str)
            self.set_attribute(file, 'dates','value', dates_str)

        ### set iteration number ###

        self.set_attribute(file, 'mcmc', 'chainLength', str(self.n_iter))

        ### set additional parameters ###
        ### .... 
        
        return file

##### Nexus and Newick format functions ######

def to_newick(G, root=None, node_attrs=None, edge_attrs=None, branch_attr="weight"):
    if root is None:
        roots = [n for n, d in G.in_degree() if d == 0]
        if len(roots) != 1:
            raise ValueError(f"Expected exactly one root, found {len(roots)}")
        root = roots[0]

    def build(node, parent):
        children = list(G.successors(node))

        # Recurse into children first
        if children:
            subtrees = ",".join(build(c, node) for c in children)
            label = f"({subtrees}){node}"
        else:
            label = f"{node}"

        # Branch length (from the edge connecting parent -> node)
        branch_len = ""
        if parent is not None and branch_attr:
            bl = G.edges[parent, node].get(branch_attr)
            if bl is not None:
                branch_len = f":{bl}"

        # NHX metadata (node attrs + edge attrs, merged)
        nhx_pairs = []
        if node_attrs:
            for a in node_attrs:
                if a in G.nodes[node]:
                    nhx_pairs.append(f"{a}={G.nodes[node][a]}")
        if edge_attrs and parent is not None:
            for a in edge_attrs:
                if a in G.edges[parent, node]:
                    nhx_pairs.append(f"{a}={G.edges[parent, node][a]}")

        nhx = f"[{':'.join(nhx_pairs)}]" if nhx_pairs else ""

        return f"{label}{branch_len}{nhx}"

    return build(root, None) + ";"

def seq_to_nexus(sw, dates, name, max_val):
    matrix = ''
    for ms, s in sw.items():
        matrix += ms+'|'+str(dates[ms])+'\t'+s+'\n'

    char_set = ''.join(map(str,range(max_val+1)))

    config_s = f'''#NEXUS

    BEGIN DATA;
        DIMENSIONS NTAX={len(sw.keys())} NCHAR={len(s)};
        FORMAT DATATYPE=STANDARD MISSING=? GAP=-  SYMBOLS="{char_set}";
    MATRIX
    {matrix}
    ;
    END;
    '''

    return config_s

def trad_to_nexus(st, seqs):
    wits_seqs = {}
    dates = {}
    for w,s in seqs.items():
        sequence = []
        for l in s.values():
            for app in l.values():
                sequence.append(app)
        wits_seqs[f'ms_{w}'] = ''.join(map(str,sequence))
        dates[f'ms_{w}'] = int(st.nodes[w]['birth_time'])

    max_val = max([max(l.values()) for l in loc_sizes.values() if l != {}])
    nexus_str =seq_to_nexus(wits_seqs, dates, 'test', max_val)
    return nexus_str