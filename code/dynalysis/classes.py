import os 
import shutil
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import dynalysis.basics as bcs
from time import sleep
from dynalysis.basics import AlgorithmError, readline, InputError, fano,\
to_string, to_float, to_int
from random import choice
from warnings import warn

#=====subfunctions=====#
def assign_neulist(tpe):
    '''
    tpe: str, ''.join([0s, 1s]), where 0s represent exc neurons and 1 inhs,
        or list, e.g. [0, 0, 1, 1]
    returns: 1) list, [post, pre, pre] in str
        2) list, [(pre, post), (pre, post)] in str
        3) list, [(pre, post), (pre, post)] in binary
    '''
    inh, exc, out=1, 1, []
    for t in tpe:
        if int(t)==0:
            out.append('exc_'+str(exc))
            exc+=1
        elif int(t)==1:
            out.append('inh_'+str(inh))
            inh+=1
        else: raise AlgorithmError('func assign neulist')
    return out, [(t, tpe[0]) for t in tpe[1:]], [(s, out[0]) for s in out[1:]]

def sort_by_value(l):
    ''' returns a sorted list of strings, sorted as if it were ints '''
    return to_string(sorted(to_int(l)))

def sort_nodes(nodelist):
    ''' sort by exc then inh, low number then high '''
    mapdic, res={},[]
    for n in nodelist:
        try:
            ei, num = n.split('_')
        except ValueError:
            ei, num = n, 0
        mapdic[(ei, num)]=n
    for key in sorted(mapdic):
        res.append(mapdic[key])
    return res

def sort_pairlist(pairlist):
    mapdic, res={}, []
    for pair in pairlist:
        e1, e2 = pair[0].split('_')
        e3, e4 = pair[1].split('_')
        k=(e1, int(e2), e3, int(e4))
        mapdic[k]=pair
    for k in sorted(mapdic.keys()):
        res.append(mapdic[k])
    return res

def test_pol(node):
    return node.split('_')[0]

def get_ID(node):
    return node.split('_')[1]

def convert_neuname(neuron, name):
    if name=='': return neuron
    ei, num=neuron.split('_')
    return ei+'_'+name+'-'+num

def preneu(target, pairlist):
    prelist=[]
    for pair in pairlist:
        if pair[1]==target: prelist.append(pair[0])
    return prelist

def postneu(target, pairlist):
    prelist=[]
    for pair in pairlist:
        if pair[0]==target: prelist.append(pair[1])
    return prelist

#=====functions regarding dynalysis topfiles====#
def topfiles(f, state='single'):
    import fnmatch
    matches = []
    dynalysis='C:\\Users\\belle\\anaconda3\\Lib\\site-packages\\dynalysis\\topfile'
    for root, dirnames, filenames in os.walk(dynalysis):
        for filename in fnmatch.filter(filenames, f):
            matches.append(os.path.join(root, filename))
    if state=='single': return matches[0]
    elif state=='all': return matches
    elif state=='parent': return matches[0][:-len(f)]
    else: raise InputError('func topfiles')
    
def dynalysis_path():
    return 'C:\\Users\\RoundyPie\\Anaconda2\\lib\\site-packages\\dynalysis'

def copy_flysim(destination=os.getcwd()):
    shutil.copy2(topfiles('flysim.out'), destination)
    
def copy_iter0(destination=os.getcwd()):
    shutil.copy2(topfiles('ss_flysim_iter=0.sh'), destination)

#=====classes=====#
class Motif():
    def __init__(self, birth=0, n=4):
        self.ID=int(birth)
        self.n = n
        self.mm, self.connectome = self.MotherMatrix()
    def MotherMatrix(self):
        n=self.n
        mothermatrix, count=[[] for i in range(n)], 0
        for i in range(n):
            for j in range(n-1):
                mothermatrix[i].append(count)
                count+=1
            mothermatrix[i].insert(i, '0')
        #MM to connectome
        connectome={}
        for i in range(n):
            for j in range(n):
                if mothermatrix[i][j]!='0':
                    connectome[mothermatrix[i][j]]=(self.neulist[i], self.neulist[j])
        return mothermatrix, connectome
    def pre_neuron(self, ind):
        prelist=[row[ind] for row in self.MotherMatrix()[0]]
        prelist.remove('0')
        return prelist
    def post_neuron(self, ind):
        postlist=self.MotherMatrix()[0][ind]
        postlist.remove('0')
        return postlist
    def post_neuron_from_conn(self, ind, connection):
        neulist, connectome=self.neulist, self.MotherMatrix()[1]
        pre, res=neulist[ind], []
        neulist.pop(ind)
        for postneu in neulist:
            neuron=connectome.keys()[connectome.values().index((pre, postneu))]
            if neuron in connection:
                res.append(postneu)
        return res
    def ID_renew(self, motif):
        renew={1:0, 2:3, 3:1, 4:6, 5:2, 6:9, 7:7, 8:4, 9:8, 10:11, 11:10, 12:5}
        return sorted([renew[m] for m in motif])
	
class motif(Motif):
    def __init__(self, birth=0, n=4, listtype='1100', name=''):
        self.ID = int(birth)
        self.name = name
        self.listtype = listtype
        self.typeID = self.update_typeID()
        self.neulist=self.update_neulist()
        self.neunamelist=self.update_neunamelist()
        self.n = len(self.neulist)
        self.mm, self.connectome = self.MotherMatrix()
        self.connlist = self.update_connlist()
        self.pairlist = self.update_pairlist()
        self.motif_character = self.update_motif_character()
    def update_all(self):
        self.neulist=sort_nodes(self.update_neulist())
        self.neunamelist=self.update_neunamelist()
        self.n = len(self.neulist)
        self.mm, self.connectome = self.MotherMatrix()
        self.connlist = self.update_connlist()
        self.pairlist = self.update_pairlist()
        self.motif_character = self.update_motif_character()
    def update_typeID(self):
        tot=0
        listtype=sorted([int(i) for i in self.listtype])
        for c in range(len(listtype)):
            tot+=listtype[-(c+1)]*2**c
        return tot
    def update_neulist(self):
        return sorted(assign_neulist(self.listtype)[0], key=str.lower)
    def update_neunamelist(self):
        l=[]
        for neu in self.neulist: l.append(convert_neuname(neu, self.name))
        return l
    def update_connlist(self):
        def max2(ID, n_start):
            n, deter=n_start, True
            while deter:
                if ID>=2**(n-1) and ID<2**n: return n-1, ID-2**(n-1)
                if n<0: raise AlgorithmError('func update_connection')
                n-=1
        connlist, n_start, ID_left = [], (self.n)**2-self.n, self.ID
        while ID_left>0:
            n_start, ID_left = max2(ID_left, n_start)
            connlist.append(n_start)
        return connlist
    def update_pairlist(self):
        pairlist=[]
        for conn in self.connlist:
            pairlist.append(self.connectome[conn])
        return sort_pairlist(pairlist)
    def update_motif_character(self):
        node_character = {}
        for neu in self.neulist:
            node_character[neu]=['self-'+neu.split('_')[0]]
        for pair in self.pairlist:
            node_character[pair[0]].append('to-'+pair[1].split('_')[0])
            node_character[pair[1]].append('from-'+pair[0].split('_')[0])
        return node_character
    def assign_ID(self, connlist): #Interactive
        if type(connlist)==int: self.ID=connlist
        else:
            ID, dic, resconn=0, self.connectome, []
            for conn in connlist:
                if type(conn)==int: upp=conn
                else: upp=dic.keys()[dic.values().index(conn)]
                resconn.append(upp)
                ID+=2**upp
            self.ID, self.connlist = ID, resconn
        self.update_all()
    def assign_neulist(self, listtype): #Interactive
        self.listtype=listtype
        self.update_all()
    def assign_pairlist(self, pairlist):
        self.pairlist = sort_pairlist(pairlist)
        neulist, connlist, listtype=[], [], []
        for pair in pairlist:
            for item in pair:
                if item not in neulist: neulist.append(item)
        for neu in neulist:
            if neu[:3]=='exc': listtype.append(0)
            elif neu[:3]=='inh': listtype.append(1)
            else: raise AlgorithmError('class motif: func assign_pairlist')
        self.listtype = ''.join(to_string(listtype))
        self.neulist = sort_nodes(neulist)
        self.n = len(self.neulist)
        self.mm, self.connectome = self.MotherMatrix()
        for pair in pairlist:
            conn=self.connectome.keys()[self.connectome.values().index(pair)]
            connlist.append(conn)
        self.connlist = sorted(connlist) 
        self.ID = sum([2**i for i in connlist])
        self.update_all()
    def draw(self): #Interactive ##too ugly
        if self.listtype=='0011': ##
            G=nx.DiGraph()
            G.add_edges_from(self.pairlist)
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.set_xticks([])
            ax1.set_yticks([])
            position = {'exc_1':(50, -100), 'exc_2':(50, -50), 'inh_1':(75, -75), 'inh_2':(100, -75)}
            ns = np.array([900, 900, 900, 900])
            nc = np.array(['g', 'r', 'r', 'g'])
            nx.draw_networkx(G, pos=position, node_size=ns, node_color=nc, ax=ax1)
            plt.savefig(str(self.ID)+'_motif.png')
            plt.clf()
        else:
            G=nx.DiGraph()
            G.add_edges_from(self.pairlist)
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.set_xticks([])
            ax1.set_yticks([])
            nx.draw(G)
            plt.savefig(str(self.ID)+'_motif.png')
            plt.clf()
    def pre_neuron_from_conn(self, post):
        neulist, connectome=self.neulist[:], self.MotherMatrix()[1]
        res=[]
        neulist.remove(post)
        for preneu in neulist:
            neuron=connectome.keys()[connectome.values().index((preneu, post))]
            if neuron in self.connlist:
                res.append(preneu)
        return sort_nodes(res)
    def output(self, dataname, namename, tpe_list):
        out_data=open(dataname, 'w')
        out_name=open(namename, 'w')
        max_exc=self.listtype.count('0')
        for pair in self.pairlist:
            if test_pol(pair[0])=='exc': p1 = str(int(get_ID(pair[0]))-1)
            if test_pol(pair[0])=='inh': p1 = str(int(get_ID(pair[0]))-1+max_exc)
            if test_pol(pair[1])=='exc': p2 = str(int(get_ID(pair[1]))-1)
            if test_pol(pair[1])=='inh': p2 = str(int(get_ID(pair[1]))-1+max_exc)
            out_data.write(p1+' '+p2+' 1\n')
        for t in range(4):
            out_name.write(tpe_list[t]+'-F-00000'+str(t)+'\n')
        out_data.close()
        out_name.close()

class entrydic():
    def __init__(self, dic={}):
        self.dic = dic
    def __len__(self):
        return len(self.dic)
    def column_len(self):
        return len(self.dic.values()[0])
    def add_row(self, key, row):
        self.dic[key]=row
    def add_column(self, foreign_dic):
        if len(foreign_dic.values())!=len(self):
            raise InputError('Column length incorrect.')
        else:
            for k in self.dic.keys():
                for item in foreign_dic[k]:
                    self.dic[k].append(item)
    def mean(self):
        mean_dic = {}
        for k in self.dic.keys():
            mean_dic[k] = np.mean(self.dic[k])
        return mean_dic
    def fano(self):
        fano_dic = []
        for k in self.dic.keys():
            fano_dic[k] = fano(self.dic[k])
        return fano_dic
    def specific_neu(self, neu_num):
        spec = {}
        for k in self.dic.keys():
            spec[k] = self.dic[k][neu_num]
        return spec
    def multi_neu(self):
        multi = [{} for m in range(self.column_len())]
        for k in self.dic.keys():
            for i in range(self.column_len()):
                multi[i][k] = self.dic[k][i]
        return multi
    def __str__(self):
        return repr(self.rows)
    
class branch():
    def __init__(self, name, parentpath):
        self.name = str(name)
        self.parentlist = self.link_to_list(parentpath) if type(parentpath)==str\
        else parentpath
        self.parentlink = self.list_to_link(parentpath) if type(parentpath)==list\
        else parentpath
        self.pathlist = self.parentlist + [self.name]
        self.pathlink = self.list_to_link(self.pathlist)
    def mkdir(self):
        if not os.path.exists(self.pathlink): os.makedirs(self.pathlink)
        sleep(0.2)
    def link_to_list(self, path):
        pathlist = path.split('\\')
        while '' in pathlist:
            pathlist.remove('')
        return pathlist
    def list_to_link(self, pathlist):
        return '\\'.join(pathlist)
    def cp_to(self, filename, destination): #Interactive
        filename = self.pathlink+'\\'+filename
        if type(destination)==str: shutil.copy2(filename, destination)
        elif type(destination)==int: shutil.copy2(filename, self.list_to_link(destination))
        else: raise InputError('func cp_to')
    def cp_from(self, filename, source): #Interactive
        os.chdir(source)
        shutil.copy2(filename, self.pathlink)
        os.chdir(self.pathlink)
    def move_to(self, filename, destination): #Interactive
        if type(destination)==str: shutil.move(self.pathlink+'\\'+filename, destination+'\\'+filename)
        elif type(destination)==int: shutil.move(self.pathlink+'\\'+filename, self.list_to_link(self.pathfile+[filename]))
        else: raise InputError('func move_to')
    def move_from(self, filename, source): #Interactive
        os.chdir(source)
        shutil.move(source+'\\'+filename, self.pathlink+'\\'+filename)
        os.chdir(self.pathlink)
    def move_self(self, destination): #Interactive
        os.chdir(self.parentlink)
        shutil.move(self.pathlink, destination+'\\'+self.name)
        os.chdir(self.parentlink)
        #Update all attributes
        self.parentlist = self.link_to_list(destination)
        self.parentlink = destination
        self.pathlist = self.parentlist + [self.name]
        self.pathlink = self.list_to_link(self.pathlist)
    def move_content_to(self, destination):
        for item in os.listdir(self.pathlink):
            shutil.copy2(item, destination)
    def get_sublist(self):
        return os.listdir(self.pathlink)
    
class entry():
    def __init__(self, key, value):
        self.keytype=type(key)
        if self.keytype==str:self.key=[key[k] for k in range(len(key))]
        elif self.keytype==list: self.key=[[val[v] for v in range(len(val))] for val in key]
        else: raise InputError('class entry')
        self.value=[[val[v] for v in range(len(val))] for val in value]
    def process_str(self, item, row):
        for p in range(len(item)/2):
            row=row.split(item[2*p])[int(item[2*p+1])]
        return row
    def process_list(self, item, row):
        res=[]
        for val in item:
            cprow=row
            for p in range(len(val)/2):
                cprow=cprow.split(val[2*p])[int(val[2*p+1])]
            res.append(cprow)
        return tuple(res) #because keys cannot be lists
    def process_all(self, row):
        return row.split('\n')[0].split(' ')
    def readdata(self, filename, quick_format=False, start=0, tofloat=False):
        data, dic=readline(filename)[start:], {}
        for row in data:
            if self.keytype==str: k=self.process_str(self.key, row)
            else: k=self.process_list(self.key, row)
            if not quick_format: v=self.process_list(self.value, row)
            else: v=self.process_all(row)[1:]
            if k not in dic.keys(): dic[k]=[]
            if not tofloat: dic[k].append(v)
            else: dic[k].append(to_float(v))
        return dic
    def merge_data(self, instant, filename, filename2, **kwargs):
        #kwargs
        kw={'qf1':False, 'qf2':False}
        kw.update(kwargs)
        quick_format1, quick_format2=kw['qf1'], kw['qf2']
        #main
        dic1=self.readdata(filename, quick_format=quick_format1)
        dic2=instant.readdata(filename2, quick_format=quick_format2)
        dic3={}
        for key in dic1.keys():
            try:
                dic3[key]=[dic1[key][0][0], dic2[key][0][0]]
            except KeyError:
                print('dic2 has no key '+str(key))
        return dic3

class output_dic():
    def __init__(self, filename, to_output):
        self.filename=filename
        self.output=to_output
    def outputdata(self, readwrite='w'):
        outfile=open(self.filename, readwrite)
        for key in sort_by_value(self.output.keys()):
            val=self.output[key]
            if type(val)==list: line=str(key)+' '+' '.join(self.output[key])
            elif type(val)==(str or int or float): line=str(key)+' '+str(val)
            else: raise InputError('func output')
            outfile.write(line)
            outfile.write('\n')
        outfile.close()
        
class outhelper():
    def __init__(self, outfilename):
        self.outfilename=outfilename
    def output(self, subject, rc='row'):
        outfile=open(self.outfilename, 'a')
        if type(subject)!=(list or tuple):
            outfile.write(str(subject)+'\n')
        else:
            if type(subject[0])!=(list or tuple) and rc=='row':
                outfile.write(' '.join(to_string(subject))+'\n')
            else:
                for row in subject: self.output(row, rc)
        outfile.close()
        
class network():
    def __init__(self, nodes, pairlist):
        self.nodes = nodes
        self.pairlist = pairlist
        self.max_exc, self.max_inh = self.update_maxs()
        self.size = len(nodes)
        self.connsize = len(self.pairlist)
        self.redundant_warning(self.nodes)
    def update_maxs(self):
        max_exc, max_inh = 0, 0
        for node in self.nodes:
            polarity, num = node.split('_')
            if polarity == 'exc' and int(num)>max_exc: max_exc=int(num)
            if polarity == 'inh' and int(num)>max_inh: max_inh=int(num)
        return max_exc, max_inh
    def update_ID(self, oldID):
        polarity, num = oldID.split('_')
        if polarity == 'exc': return 'exc_'+str(int(num)+self.max_exc)
        if polarity == 'inh': return 'inh_'+str(int(num)+self.max_inh)
    def add_node(self, num): #Interactive
        for n in range(num):
            self.max_exc+=1
            self.nodes.append('exc_'+str(self.max_exc))
        self.max_exc, self.max_inh = self.update_maxs()
        self.size = len(self.nodes)
        self.redundant_warning(self.nodes)
    def add_motif(self, motif): #Interactive
        for node in motif.neulist:
            self.nodes.append(self.update_ID(node))
        for pair in motif.pairlist:
            self.pairlist.append((self.update_ID(pair[0]), self.update_ID(pair[1])))
        self.max_exc, self.max_inh = self.update_maxs()
        self.size = len(self.nodes)
        self.connsize = len(self.pairlist)
        self.redundant_warning(self.nodes)
        self.redundant_warning(self.pairlist)
    def add_connection(self, pair): #Interactive
        for neuron in pair:
            if neuron not in self.nodes: self.nodes.append(neuron)
        self.pairlist.append(pair)
        self.size = len(self.nodes)
        self.connsize = len(self.pairlist)
        self.redundant_warning(self.nodes)
        self.redundant_warning(self.pairlist)
    def add_random_connection(self, iteration=1, allow_recurrent=False, tolerance=0.1): #Interactive
        if iteration >= self.size*(self.size-1)-len(self.pairlist)-self.size*tolerance:
            raise bcs.AlgorithmError('class network: iteration is too high.') 
        for it in range(int(iteration)):
            node1, node2 = choice(self.nodes), choice(self.nodes)
            if (node1, node2) not in self.pairlist:
                if not allow_recurrent and (node2, node1) not in self.pairlist:
                    self.add_connection((node1, node2))
                if allow_recurrent: self.add_connection((node1, node2))
        self.connsize = len(self.pairlist)
        self.redundant_warning(self.nodes)
        self.redundant_warning(self.pairlist)
    def remove_connection(self, pair): #Interactive
        self.pairlist.remove(pair)
        self.connsize = len(self.pairlist)
    def remove_random_connection(self, iteration=1): #Interactive
        for it in range(int(iteration)):
            pair = choice(self.pairlist)
            self.pairlist.remove(pair)
        self.connsize = len(self.pairlist)
    def redundant_warning(self, nodes):
        if len(list(set(nodes)))!=len(nodes): raise bcs.AlgorithmError('Redundant nodes detected!')
    def draw(self): #Interactive
        G=nx.DiGraph()
        G.add_edges_from(self.pairlist)
        pos={}
        for node in self.nodes:
            pos[node]=(int(node.split('_')[1])//25, int(node.split('_')[1])%25)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_xticks([])
        ax1.set_yticks([])
        nx.draw_networkx(G, ax=ax1, node_size=2, pos=pos, with_labels=False, width=0.2)
        plt.savefig('network.png', dpi=200)
        plt.clf()
    def output(self, dataname, namename, exctype='VGlut', inhtype='Gad1'): #Interactive
        out_data = open(dataname, 'w')
        out_name = open(namename, 'w')
        for pair in self.pairlist:
            if test_pol(pair[0])=='exc': p1 = str(int(get_ID(pair[0]))-1)
            if test_pol(pair[0])=='inh': p1 = str(int(get_ID(pair[0]))-1+self.max_exc)
            if test_pol(pair[1])=='exc': p2 = str(int(get_ID(pair[1]))-1)
            if test_pol(pair[1])=='inh': p2 = str(int(get_ID(pair[1]))-1+self.max_exc)
            out_data.write(p1+' '+p2+' 1\n')
        pending = []
        for node in self.nodes:
            polarity, num = node.split('_')
            if polarity=='exc':
                num = str(int(num)-1)
                out_name.write(exctype+'-F-'+'0'*(6-len(num))+num+'\n')
            if polarity=='inh':
                num = str(int(num)-1+self.max_exc)
                pending.append(inhtype+'-F-'+'0'*(6-len(num))+num+'\n')
        for pend in pending:
            out_name.write(pend)
        out_data.close()
        out_name.close()
        
class rand_motif_network():
    def __init__(self, ID, netsize, rand_degree):
        self.net = network([],[])
        self.ID = ID
        self.netsize = netsize
        self.rand_degree = rand_degree
    def build_rand_motif_network(self, remove=True, percent=30):
        mot=motif(self.ID)
        ssize=self.netsize//len(mot.neulist)
        print('Number of motifs in network = '+str(ssize))
        for s in range(ssize):
            self.net.add_motif(mot)
        self.net.add_random_connection(iteration=self.net.size*(self.net.size-1)*self.rand_degree/100.0-self.net.connsize)
        if remove:
            self.net.remove_random_connection(iteration=len(self.net.pairlist)*percent/100.0)
    def build_concise(self):
        ssize=self.netsize//4
        concise_network = network([], [])
        concise_network.add_node(ssize)
        for pair in self.net.pairlist:
            pre, post = pair
            concise_pre = (int(pre.split('_')[1])+1)//2
            concise_post = (int(post.split('_')[1])+1)//2
            to_test = ('exc_'+str(concise_pre), 'exc_'+str(concise_post))
            if to_test not in concise_network.pairlist:
                concise_network.add_connection(to_test)
        return concise_network
    
class M_network():
    def __init__(self, nodes, pairlist):
        self.nodes = nodes
        self.neulist = sort_nodes(self.nodes)
        self.pairlist = pairlist
        self.size = len(self.nodes)
        self.connsize = len(self.pairlist)
        self.ut, self.bias = {}, {}
    def update_all(self):
        self.nodes = self.nodes
        self.neulist = sort_nodes(self.nodes)
        self.pairlist = self.pairlist
        self.size = len(self.nodes)
        self.connsize = len(self.pairlist)
        self.ut, self.bias = self.ut, self.bias
        self.redundant_warning(self.nodes)
        self.redundant_warning(self.pairlist)
    def update_ut_by_motif(self, motif, synwgt):
        #e1e2,e1i1,e2e1,e2i2,i1e1,i1i2,i2e2,i2i1
        c, reclist=0, ['AMPA','AMPA','AMPA','AMPA','GABA','GABA','GABA','GABA']
        for pair in motif.pairlist:
            pre, post = convert_neuname(pair[0], motif.name), convert_neuname(pair[1], motif.name)
            if pre not in self.ut.keys(): self.ut[pre]=[]
            self.ut[pre].append([post,reclist[c],1,synwgt[c],1])
            c+=1
    def update_ut_by_connection(self, pre, post, rec, synwgt):
        if pre not in self.ut.keys(): self.ut[pre]=[]
        self.ut[pre].append([post,rec,1,synwgt,1])
    def update_bias_by_motif(self, motif, bias):
        #e1,e2,i1,i2
        c=0
        for neu in motif.neulist:
            neuron = convert_neuname(neu, motif.name)
            if (0, neuron) in self.bias.keys():
                warn('Repeated bias current detected!') 
                print(self.bias.keys())
            self.bias[(0,neuron)]=['ramp',bias[c],0]
            c+=1
    def add_motif(self, motif): #Interactive
        for node in motif.neulist:
            self.nodes.append(convert_neuname(node, motif.name))
        for pair in motif.pairlist:
            self.pairlist.append((convert_neuname(pair[0], motif.name), convert_neuname(pair[1], motif.name)))
        self.update_all()
    def add_connection(self, pair): #Interactive
        for neuron in pair:
            if neuron not in self.nodes: self.nodes.append(neuron)
        self.pairlist.append(pair)
        self.update_all()
    def remove_connection(self, pair): #Interactive
        self.pairlist.remove(pair)
        self.update_all()
    def redundant_warning(self, nodes):
        if len(list(set(nodes)))!=len(nodes): raise bcs.AlgorithmError('Redundant nodes detected!')
    def draw(self): #Interactive
        G=nx.DiGraph()
        G.add_edges_from(self.pairlist)
        fig = plt.figure(figsize=(50,50))
        ax1 = fig.add_subplot(111)
        ax1.set_xticks([])
        ax1.set_yticks([])
        nx.draw_networkx(G, pos=nx.shell_layout(G), ax=ax1, node_size=2, with_labels=True, width=0.2)
        plt.savefig('network.png', dpi=200)
        plt.clf()
    
if __name__=='__main__':
    pass