import numpy as np
import networkx as nx
import os
from dynalysis.basics import readline, InputError, AlgorithmError, output_line
from dynalysis.classes import motif

def print_info(name, value):
    if len(name)!=len(value):
        print(name)
        print(value)
        raise InputError('function print_info')
    else:
        res=[]
        for n in range(len(name)):
            if type(value[n])!=list:
                res.append(name[n]+'='+str(value[n]))
            else:
                res.append(name[n]+'='+'_'.join(map(lambda x:str(x), value[n])))
        return res
    
def sort_nodes(nodelist):
    mapdic, res={},[]
    for n in nodelist:
        try:
            ei, num = n.split('_')
        except ValueError:
            ei, num = n, 0
        mapdic[(ei, num)]=n #originally int(num), removed for general names
    for key in sorted(mapdic):
        res.append(mapdic[key])
    return res

class Neuron():
    def __init__(self, name, intrinsic, std, receptor, target, LTP):
        self.name=name
        self.intrinsic=intrinsic
        self.std=std
        self.receptor=receptor
        self.target=target
        self.LTP = LTP
        self.in_properties=['N', 'C', 'Taum', 'RestPot', 'ResetPot', 'Threshold']
        self.std_properties=['STD_pv', 'STD_tD', 'STD_pv1', 'STD_tD', 'STD_pv2', 'STD_tD', 'STD_pv3', 'STD_tD']
        self.rec_properties=['Receptor', 'Tau', 'RevPot', 'FreqExt', 'MeanExtEff', 'MeanExtCon']
        self.target_properties=['TargetPopulation', 'TargetReceptor', 'MeanEff', 'weight', 'STDType']
        self.LTP_properties=['LTP_tLP', 'LTP_PosISI', 'LTP_NegISI']
    def __str__(self):
        res=[self.name]
        for tar in self.target:
            res.append(tar[0]+', '+str(tar[2]))
        return '\n'.join(res)
    def update_properties(self, item):
        pass
    def update_target(self, updatelist):
        cp_target = self.target[:]
        for tar in cp_target:
            if tar[0] in [s[0] for s in updatelist]:
                index=self.target.index(tar)
                self.target.remove(tar)
                for stuff in updatelist:
                    if stuff[0]==tar[0]:
                        updatelist.remove(stuff)
                        self.target.insert(index,stuff)
        for remain in updatelist:
            self.target.append(remain)
        return self.target
    def output(self, filename):
        #Intrinsic and std properties
        res = ['NeuralPopulation:'+str(self.name)]+print_info(self.in_properties, self.intrinsic)+\
        [''] #+print_info(self.std_properties, self.std)+['']
        #LTP properties
        res += print_info(self.LTP_properties, self.LTP)+['']
        #Receptor properties
        for rec in self.receptor:
            res += ['Receptor:'+str(rec[0])]
            res += print_info(self.rec_properties, rec)[1:]+['EndReceptor', '']
        #Target properties
        for tar in self.target:
            res += ['TargetPopulation:' + str(tar[0])]
            res += print_info(self.target_properties, tar)[1:]+['EndTargetPopulation', '']
        #output
        res += ['EndNeuralPopulation','%--------------------------------------']
        
        output = open(filename, 'a')
        for info in res:
            output.write(info+'\n')
            #output_line(filename, info)
            #print >> output, info
        output.close()
        return res
    
def data_specification(kw):
    intrinsic=[kw['N'],kw['C'],kw['Taum'],kw['RestPot'],\
    kw['ResetPot'],kw['Threshold']]
    std=kw['STD']
    receptor=[kw['AMPA'],kw['NMDA'],kw['GABA'],kw['Ach'],kw['Glucl']]
    LTP=kw['LTP']
    return intrinsic, std, receptor, LTP
    
def EI(pre, post, to_deter):
    pre, post = pre.split('_')[0], post.split('_')[0]
    for index in range(4):
        if (pre, post)==[('exc', 'inh'), ('exc', 'exc'), ('inh', 'exc'), ('inh', 'inh')][index]:
            return to_deter[index]

def exec_conf(meanEff,**kwargs):
    #kwargs
    kw = {'ID':2447, 'filename':'network.conf', 'event_dic':{},\
          'update_tar':{}, 'add_connections':[], 'graph':[], 'Mot':None,\
          'suppress':False}
    kw.update(kwargs)
    ID, filename, event_dic = kw['ID'], kw['filename'], kw['event_dic']
    update_tar, add_connections = kw['update_tar'], kw['add_connections']
    graph, Mot, suppress = kw['graph'], kw['Mot'], kw['suppress']
    #CM
    if Mot!=None:
        M=Mot
        graph = nx.DiGraph()
        graph.add_edges_from(M.pairlist+add_connections)
    elif ID!=False or ID==0:
        M=motif(ID)
        graph = nx.DiGraph()
        try:
            graph.add_edges_from(M.pairlist+add_connections)
        except ValueError:
            print('Connectivity Matrix input incorrect, please try again.')
    else:
        graph=kw['graph']
    #data specification:
    kk={'N':1, 'C':0.5, 'Taum':20.0, 'RestPot':-70.0, 'ResetPot':-55.0,\
        'Threshold':-50, 'STD':[0.22, 125, 0.22, 125, 0.22, 125, 0.22, 125],\
        'AMPA':['AMPA','2','0','0','10.5','1'],\
        'NMDA':['NMDA','100','0','0','0','0'],\
        'GABA':['GABA','5','-90','0','0','0'],\
        'Ach':['Ach','20','0','0','3.11627897405','1'],\
        'Glucl':['Glucl','5','-90','0','0','0'], 'meanEff':meanEff,\
        'LTP':[30, 0, 0]}
    del kw['event_dic'] #For faster computation, meaningless algorithm-wise
    del kw['update_tar']
    kk.update(kw)
    change, neulist = {},list(graph.nodes())
    if os.path.isfile('./conf_info.txt'):
        if not suppress: print('conf_info.txt exists, changing data information.')
        data=readline('conf_info.txt')
        for row in data:
            kwcp = kk.copy()
            for e in row.split(' ')[1:]:
                param, val = e.split('=')
                splitval = val.split(',')
                if splitval[0]==val:
                    kwcp[param] = val
                else:
                    kwcp[param] = splitval
            change[row.split(' ')[0]]=kwcp
            neulist.remove(row.split(' ')[0])
    elif event_dic=={}:
        pass
        #print('No specified change in data information.')
    for node in neulist:
        change[node]=kk
    for node in event_dic.keys():
        if not suppress: print('event_dic for '+node+' exists, changing data information.')
        external = change[node].copy()
        external.update(event_dic[node])
        change[node] = external
        
    #.conf
    to_clear = open(filename, 'w')
    to_clear.close()
    for node in sort_nodes(graph.nodes()):
        target = []
        intrinsic, std, receptor, LTP = data_specification(change[node])
        for post in sort_nodes(graph.successors(node)):
            if post not in change[node].keys(): #weight=1, not 7.5
                target.append([post, EI(node,post,['AMPA','AMPA','GABA','GABA']),\
                               EI(node,post,change[node]['meanEff']), 1, EI(node,post,['0', '1', '2', '3'])])
            else:
                target.append([post, EI(node,post,['AMPA','AMPA','GABA','GABA']),\
                               change[node][post].split('\n')[0], 1, EI(node,post,['0', '1', '2', '3'])])
        neuron=Neuron(node, intrinsic, std, receptor, target, LTP)
        #update_tar
        if node in update_tar.keys():
            neuron.update_target(update_tar[node])
            if not suppress: print('For neuron '+node+', target information is updated.')
        #output
        neuron.output(filename)
    return graph