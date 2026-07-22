import os
import birth_death_utils as bd
import networkx as nx
import numpy as np
from lxml import etree as et

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