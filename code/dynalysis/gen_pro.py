import os
import numpy as np
import networkx as nx
from dynalysis.basics import readline, output_single_col, output_line, output_clf, InputError
from dynalysis.classes import motif

def membrane_current(filename, pop_name, event_t, gmean, gstd):
    l = ['EventTime '+str(event_t),\
         'Type=ChangeMembraneNoise', 'Label=#1#',\
         'Population:' + str(pop_name),\
         'GaussMean=' + str(gmean), 'GaussSTD=' + str(gstd),\
         'EndEvent', '%--------------------------------------', '']
    output_single_col(filename, l, readwrite='a')
    return 0

def external_stimulus(filename, pop_name, event_t, receptor, fire_r):
    l = ['EventTime ' + str(event_t),\
         'Type=ChangeExtFreq','Label=#1#',\
         'Population:' + str(pop_name),\
         'Receptor: ' + str(receptor),\
         'FreqExt=' + str(fire_r),\
         'EndEvent','%--------------------------------------','']
    output_single_col(filename, l, readwrite='a')
    return 0

def Vm_output (filename, output_filename):
    l = ['OutControl', 'FileName:' + str(output_filename),\
         'Type=MemPot','population:AllPopulation', 'EndOutputFile',\
         '%--------------------------------------']
    output_single_col(filename, l, readwrite='a')
    return 0

def firing_r_output (filename, output_filename):
    l = ['OutControl', 'FileName:' + str(output_filename),\
         'Type=FiringRate','FiringRateWindow=50','PrintStep=10',\
         'population:AllPopulation','EndOutputFile',\
         '%--------------------------------------']
    output_single_col(filename, l, readwrite='a')
    return 0

def spike_output (filename, output_filename):
    l = ['OutControl','FileName:' + str(output_filename),\
         'Type=Spike','population:AllPopulation',\
         'EndOutputFile','%--------------------------------------']
    output_single_col(filename, l, readwrite='a')
    return 0

def synwgt_output (filename, output_filename):
    l = ['OutControl','FileName:' + str(output_filename),\
         'Type=SynapticWeight','population:AllPopulation',\
         'EndOutputFile','%--------------------------------------']
    output_single_col(filename, l, readwrite='a')
    return 0

def limitation_output (filename, output_filename):
    l = ['OutControl','FileName:' + str(output_filename),\
         'Type=DeltaMemPotLimit','population:AllPopulation',\
         'EndOutputFile','%--------------------------------------']
    output_single_col(filename, l, readwrite='a')
    return 0

def group_define (filename, groupname, groupneuron):
    output = open(filename, 'a')
    print >> output, 'GroupName:' + str(groupname)
    print >> output, 'GroupMembers:',
    for i in groupneuron:
        print >> output, groupneuron[i] + ',',
    print >> output, 'EndGroupMembers'
    print >> output, '%--------------------------------------'
    output.close()
    
def stimulus(filename, key, l, trial_t):
    cmd=l[0]
    if cmd=='ramp':
        membrane_current(filename, key[1], key[0], l[1], l[2])
    elif cmd=='pulse':
        membrane_current(filename, key[1], key[0], l[2], l[3])
        membrane_current(filename, key[1], key[0]+float(l[1]), 0, l[3])
    elif cmd=='pulse_b':
        membrane_current(filename, key[1], key[0], l[2], l[3])
        membrane_current(filename, key[1], key[0]+float(l[1]), l[4], l[3])
    elif cmd=='ramp_b':
        membrane_current(filename, key[1], key[0], l[1]+l[3], l[2])
    elif cmd=='multi_ramp':
        neulist=l[1].split('/')
        curlist=l[2].split('/')
        for i in range(len(neulist)):
            membrane_current(filename, neulist[i], key[0], curlist[i], l[3])
    elif cmd=='ext':
        external_stimulus(filename, key[1], key[0], l[1], l[2])
    elif cmd=='sine':
        nums = trial_t/l[3]
        x = np.arange(nums)
        t = np.linspace(l[4], trial_t, nums)
        y = [np.sin(2*np.pi*l[1] * (float(i)/nums))*l[2] for i in x]
        for j in range(len(y)):
            membrane_current(filename, key[1], t[j], y[j], l[5])
        #density
        if l[6]:
            import matplotlib.pyplot as plt
            plt.plot(t, y, 'b.')
            plt.xlabel('time (ms)')
            plt.ylabel('Input Current (A)')
            plt.ylim(-1, 1)
            plt.savefig('Input Current.png')
    
def exec_pro(**kwargs):
    #kwargs
    kw = {'filename':'network.pro', 'ID':2447, 'event_dic':{},\
          'output_Frate':True, 'output_Spike':True, 'output_MemPot':False,\
          'gmean':0.25, 'gstd':1.5, 'trial_t':20000, 'Frate':'Frate.txt',\
          'graph':[], 'Mot':None, 'output_synwgt':False,'add_connections':[],\
          'baseline':False}
    kw.update(kwargs)
    filename, ID, event_dic = kw['filename'], kw['ID'], kw['event_dic']
    output_Frate, output_Spike, output_MemPot = kw['output_Frate'], \
                                                kw['output_Spike'], kw['output_MemPot']
    gmean, gstd, trial_t = kw['gmean'], kw['gstd'], kw['trial_t']
    Frate, Mot, baseline = kw['Frate'], kw['Mot'], kw['baseline']
    output_synwgt, add_connections=kw['output_synwgt'],kw['add_connections']
    
    #CM
    if Mot!=None:
        M=Mot
        graph = nx.DiGraph()
        graph.add_edges_from(M.pairlist+add_connections)
    elif ID!=False:
        M=motif(ID)
        graph = nx.DiGraph()
        try:
            graph.add_edges_from(M.pairlist+add_connections)
        except ValueError:
            print('Connectivity Matrix input incorrect, please try again.')
    else:
        graph=kw['graph']
        
    #data specification
    diary, neulist = event_dic.copy(), list(graph.nodes())
    if baseline: #baseline activity
        for node in neulist: diary[(1, node)]=['ramp',gmean,gstd] 
    if os.path.isfile('pro_info.txt'):
        events=readline('pro_info.txt')
        for event in events:
            time, neuron = float(event.split(' ')[0]), event.split(' ')[1]
            cmd, val = event.split(' ')[2].split('=')
            diary[(time, neuron)]=[cmd]+map(lambda x:float(x),val.split(','))
    
    #stimulation
    output_clf(filename)
    for key in sorted(diary.keys()):
        stimulus(filename, key, diary[key], trial_t)
        
    #end trial
    l = ['EventTime ' + str(trial_t),\
         'Type=EndTrial','Label=End_of_the_trial','EndEvent','',\
         '%--------------------------------------']
    output_single_col(filename, l, readwrite='a')

    #define group
    group_filenames = []
    output_line(filename, 'DefineMacro', readwrite='a')
    for i in group_filenames:
        group_define(filename, i, 'something')
    output_line(filename, 'EndDefineMacro', readwrite='a')
    output_line(filename, '%--------------------------------------', readwrite='a')

    #neuron output
    if output_Frate:
        firing_r_output(filename, Frate)
    if output_Spike:
        spike_output(filename, 'Spike.txt')
    if output_MemPot:
        Vm_output(filename, 'MemPot.dat')
    if output_synwgt:
        synwgt_output(filename, 'SynWgt.txt')

    #end protocal
    output_line(filename, 'EndOutControl', readwrite='a')
    output_line(filename, '%--------------------------------------', readwrite='a')
    return graph
    
if __name__=='__main__':
    pass