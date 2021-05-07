'''===========================================
Plotting parameter sweep for switch (bistability) and decision making.
This is the same code as Fig_1_new_3435.py, but with the parameters setup so that it runs for motif B instead of A.
==========================================='''

import os
import numpy as np
import random as rand
import seaborn as sns
import matplotlib.pyplot as plt
import dynalysis.basics as bcs
import dynalysis.classes as clss
import dynalysis.visualization as vis
from scipy import stats
from itertools import combinations, permutations
from matplotlib.colors import LinearSegmentedColormap
from dynalysis.gen_conf import exec_conf
from dynalysis.gen_pro import exec_pro

N = 2 #Number of excitatory neurons / inhibitory neurons
EElist = list(range(0,205,5)) #Sweep over g_ee
IElist = list(range(0,65,5)) #Sweep over g_ie
IIlist = [0] #Sweep over g_ii
duration, amp, noise = 100, 5, 0 #Duration, amplitude and noise std for pulses used in all trials
small2perc = {0.95:0.6, 0.9:0.7, 0.8:0.7, 1:0.7}
repeat = 30
motherpath = os.getcwd()

def bistable(mot, motherpath=motherpath):
    '''
    Makes file for bistable trials to be run.
    mot: class motif
    '''
	ed = {}
	for n in range(N):
		ed[(2+n*1000, 'exc_'+str(n+1))] = ['pulse',duration,amp,noise]
		for m in range(N): ed[(502+n*1000, 'exc_'+str(m+1))]=['pulse',duration,-amp,noise]
	b1 = clss.branch(str(mot.ID)+'_b', motherpath)
	for ee in EElist:
		for ie in IElist:
			for ii in IIlist:
				b2 = clss.branch('EE='+str(ee)+'_IE='+str(ie)+'_II='+str(ii), b1.pathlink)
				b2.mkdir()
                os.chdir(b2.pathlink)
				exec_conf([30, ee, ie, ii], Mot=mot, suppress=True)
				exec_pro(Mot=mot, trial_t=1502+(N-1)*1000, output_Spike=False, event_dic=ed, baseline=True,\
						gmean=0, gstd=noise)
	os.chdir(b1.pathlink)
	vis.edplot(ed)
	clss.copy_flysim(destination=b1.pathlink)
    return 0
	
def decision(mot, small, motherpath=motherpath):
	'''
    Makes file for decision making trials to be run.
    mot: class motif
    small: float, the percentage of the smaller input compared to the larger input
    '''
	ed = {}
	combs = sorted(permutations(mot.neulist[:N], 2))
	for c in range(len(combs)):
		ed[(2+c*1000, combs[c][0])] = ['pulse',duration,amp,noise]
		ed[(2+c*1000, combs[c][1])] = ['pulse',duration,amp*small,noise]
		for m in range(N): ed[(502+c*1000, 'exc_'+str(m+1))]=['pulse',duration,-amp,noise]
	b1 = clss.branch(str(mot.ID)+'_small='+str(small)+'_d', motherpath)
	for ee in EElist:
		for ie in IElist:
			for ii in IIlist:
				b2 = clss.branch('EE='+str(ee)+'_IE='+str(ie)+'_II='+str(ii), b1.pathlink)
				b2.mkdir()
				os.chdir(b2.pathlink)
				exec_conf([30, ee, ie, ii], Mot=mot, suppress=True)
				exec_pro(Mot=mot, trial_t=1502+(len(combs)-1)*1000, output_Spike=False, event_dic=ed, baseline=True,\
						gmean=0, gstd=noise)
	os.chdir(b1.pathlink)
	vis.edplot(ed, trial_t = 1502+(len(combs)-1)*1000)
	clss.copy_flysim(destination=b1.pathlink)
    return 0
	
def def_bistable(FRs):
    '''
    Determines whether the trial has bistability.
    FRs: array-like, 2-dimensional array containing firing rate of neurons.
    '''
	for m in range(0,N):
		FR = FRs[m]
		for n in range(N):
			if np.mean(FR[100*n+31:100*n+51]) > 10:return True	
	return False

def def_decision2(FRs, mot, tpe):
    '''
    Determines whether the trial makes decisions or not.
    FRs: array-like, 2-dimensional array containing firing rate of neurons.
    mot: class motif
    tpe: string, 'exc' or 'inh', representing whether to evaluate the excitatory or inhibitory neurons
    '''
	def significant(val1, val2):
		if abs(val1-val2) > 10: return True
		else: return False
	dic, neulist = {}, mot.neulist
	perms = sorted(permutations(neulist[:N], 2))
	if tpe=='exc': combs = sorted(combinations(neulist[:N], 2))
	elif tpe=='inh': combs = sorted(combinations(neulist[N:], 2))
	else: raise bcs.InputError('tpe can only be exc, inh')
	for p in range(len(perms)):
		for c in range(len(combs)):
			large, small = combs[c] #e.g. exc_1, exc_2
			prim = np.mean(FRs[neulist.index(large)][100*p+6:100*p+11])
			secd = np.mean(FRs[neulist.index(small)][100*p+6:100*p+11])
			if significant(prim, secd):
				#((e1, e2), (e3, e4)) = (neurons the inputs are given to, neurons that are judged)
				if prim > secd: dic[(perms[p], combs[c])] = combs[c][0]
				else: dic[(perms[p], combs[c])] = combs[c][1]
			else: dic[(perms[p], combs[c])] = 'neither'
	for key in sorted(dic.keys()):
		inps, decs = key
		res1 = dic[(inps, decs)]
		res2 = dic[((inps[1], inps[0]), decs)]
		if (res1!='neither') and (res2!='neither') and res1!=res2: return True
	return False
	
def heatmap(ID, tpe, small, xaxis, yaxis, outputfigname, motherpath=os.getcwd()):
    '''
    Plots the result of the parameter sweep.
    ID: int, ID of the motif
    tpe: string, 'exc' or 'inh', representing whether to evaluate the excitatory or inhibitory neurons
    small: float, the percentage of the smaller input compared to the larger input
    xaxis: array_like, the x axis of the parameter sweep
    yaxis: array_like, the y axis of the parameter sweep
    outputfigname: string, name of the output figure
    '''
	b1 = clss.branch(str(ID)+'_b', motherpath)
	b3 = clss.branch(str(ID)+'_small='+str(small)+'_d', motherpath)
	ylist = []
	yaxis.reverse()
	for y in yaxis:
		xlist = []
		for x in xaxis:
			b2 = clss.branch('EE='+str(y)+'_IE='+str(x)+'_II='+str(0), b1.pathlink) 
			b4 = clss.branch('EE='+str(y)+'_IE='+str(x)+'_II='+str(0), b3.pathlink) 
			FRs = bcs.readcolumn(os.path.join(b2.pathlink, 'Frate.txt'), to_float=True)[1:]
			b = def_bistable(FRs)
			all_d = []
			for rp in range(1,repeat+1):
				fname = os.path.join(b4.pathlink, bcs.filename_generator('Frate.txt', rp))
				FRs = bcs.readcolumn(os.path.join(b4.pathlink, fname), to_float=True)[1:]
				temp_d = def_decision2(FRs, clss.motif(ID, n=2*N, listtype='0'*N+'1'*N), tpe)
				all_d.append(temp_d)
			if all_d.count(True) >= repeat*small2perc[small]: d=True
			else: d=False
			if b and d: deter=3
			elif not b and d: deter=2
			elif b and not d: deter=1
			else: deter=0
			xlist.append(deter)
		ylist.append(xlist)
        
	#plot
	myColors = ('whitesmoke','royalblue','lightcoral')
	cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
	ax = sns.heatmap(ylist, xticklabels=xaxis, yticklabels=yaxis, cmap=cmap)
	new_yaxis = list(np.arange(0,200+10,20))
	new_yaxis.reverse()
	plt.xticks(np.arange(.5,12.5+3,3), np.arange(0,60+15,15))
	plt.yticks(np.arange(.5,40.5+4,4), new_yaxis)
	colorbar = ax.collections[0].colorbar
	colorbar.set_ticks([])
	plt.plot()
	plt.savefig(outputfigname, dpi=200)
	plt.clf()
	return ylist


if __name__ == '__main__':
    #=====Makes files to run=====#
    if True:
        decision(clss.motif(1791), 0.8)
        bistable(clss.motif(1791))
    #====Analyzes the results=====#
    if False:
        heatmap(1791, 'exc', 0.8, IElist, EElist, '1791_small=0.8_pass='+str(small2perc[0.8])+'.png')


