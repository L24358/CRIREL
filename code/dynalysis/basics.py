import numpy as np
import warnings
import itertools

class InputError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return str(self.value)

class AlgorithmError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return str(self.value)

def readline(filename):
    ''' filename: string, name of file
    returns: the content of the file (in lines) as a list '''   
    with open(filename) as f:
        data = f.readlines()
    return data

def readcolumn(filename, to_float=False, to_int=False):
    ''' args: to_int, to_float, if else returns str '''
    data=readline(filename)
    col_num=len(data[0].split(' '))
    res=[[] for i in range(col_num)]
    for row in data:
        l = row.split('\n')[0].split(' ')
        for i in range(col_num):
            if not to_float and not to_int: res[i].append(l[i])
            elif to_float: res[i].append(float(l[i]))
            elif to_int: res[i].append(int(l[i]))
    return res

def readcolumn_uneven(filename, seperator='\t'):
    ''' args: seperator, '\t' or space '''
    res=[[]]
    data=readline(filename)
    for row in data:
        rowlist = [item for item in row.split('\n')[0].split(seperator) if item!='']
        if len(res) < len(rowlist): res += [[np.nan] *len(res[0]) for i in range(len(rowlist)-len(res))]
        for j in range(len(rowlist)):
            res[j].append(rowlist[j])
    return res

def fano(l, enable_warning=False):
    warnings.filterwarnings('error')
    try:
        return np.std(l)**2/np.mean(l)
    except (ZeroDivisionError, RuntimeWarning):
        if enable_warning: print('Error: function Fano: Zero Division Error.')
        return 0

def spiketrain(filename):
    ''' returns t and a list of spiketrains '''
    tdata, spikedata=readcolumn(filename)
    tdata=list(map(lambda x:float(x), tdata))
    spikedata=list(map(lambda x:int(x), spikedata))
    nlist=[[] for i in range(max(spikedata)+1)]
    
    for s in range(len(tdata)):
        t=tdata[s]
        spike=spikedata[s]
        nlist[spike].append(t)        
    return nlist

def filename_generator(fname, i):
    return fname if i==1 else fname+'_'+str(i)

def yieldSubList(l):
    for sublist in l:
        yield sublist
        
def output_clf(filename):
    output=open(filename, 'w')
    output.write('')
    output.close()
        
def output_line(filename, line, readwrite='a'):   
    output=open(filename, readwrite)
    output.write(line+'\n')
    output.close()
    return 0

def output_line2(output, line):   
    output.write(line+'\n')
        
def output_single(filename, l, readwrite='a'):
    output = open(filename, readwrite)
    for item in l[:-1]:
        output.write(str(item)+' ')
    output.write(str(l[-1])+'\n')
    output.close()
    return 0

def output_single_col(filename, l, readwrite='a'):
    output = open(filename, readwrite)
    for item in l:
        output.write(str(item)+'\n')
    output.close()
    return 0

def output_single_col2(output, l):
    for item in l: output.write(str(item)+'\n')
		
def output_double(filename, l, readwrite='a'):
    output = open(filename, readwrite)
    for item in l:
        for it in item[:-1]:
            output.write(str(it)+' ')
        output.write(str(item[-1])+'\n')
    output.close()
    return 0
			
def output_zero(filename, l, readwrite='a'):
    output = open(filename, readwrite)
    for item in l:
        output.write(item)
    output.close()
    return 0

def output_slashn(filename, readwrite='a'):
    output=open(filename, readwrite)
    output.write('\n')
    output.close()
    return 0

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
    
def to_string(l):
    try:
        return list(map(lambda x:str(x), l))
    except TypeError:
        print(l)
        raise InputError('func to_string')
        
def to_float(l):
    try:
        return list(map(lambda x:float(x), l))
    except TypeError:
        print(l)
        raise InputError('func to_float')
        
def to_int(l):
    try:
        return list(map(lambda x:int(x), l))
    except TypeError:
        print(l)
        raise InputError('func to_int')
		
def param_files(currents):
    l = list(itertools.product(*currents))
    return l

def read_difference(col1,col2,fname='Frate.txt'):
    all_cols=readcolumn(fname)
    time,list1,list2=all_cols[0],all_cols[col1],all_cols[col2]
    return to_float(time), np.array(to_float(list1))-np.array(to_float(list2))

def read_sum(col1,col2,fname='Frate.txt'):
    all_cols=readcolumn(fname)
    time,list1,list2=all_cols[0],all_cols[col1],all_cols[col2]
    return to_float(time), np.array(to_float(list1))+np.array(to_float(list2))

def zipp(nlist, plist, med='_'):
    res = []
    for i in range(len(nlist)):
        res.append(str(nlist[i])+'='+str(plist[i]))
    return med.join(res)