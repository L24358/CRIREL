import os
import warnings
import numpy as np
import random as rand
import matplotlib.pyplot as plt
import dynalysis.basics as bcs
import dynalysis.classes as clss
from itertools import combinations
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA, IncrementalPCA
from dynalysis.visualization import subplot
from sub_basics import *

combs=[p for p in combinations([0,1,2,3],1)]+[p for p in combinations([0,1,2,3],2)]+\
        [p for p in combinations([0,1,2,3],3)]+[p for p in combinations([0,1,2,3],4)]

def warn(*args, **kwargs):
    pass
warnings.warn = warn

def filter_data(l):
	count, res = 0, []
	for comb in combs:
		res += l[40+100*count:90+100*count]
		count += 1
	return res

def get_state_train(b, runnum, repeat=1):
	'''
	state_train is a list of coordinates.
	returns: 1) list, a list of all coordinates.
			 2) list, a list containing sub-lists of coordinates, with each sub-list being a different trial.
			 3) list, a list containing sub-lists of firing rates for each neuron
	'''
	os.chdir(b.pathlink)
	state_train_all, state_train_trucks, state_neurons = [], [], [[] for i in range(4)]
	for rp in range(1,repeat+1):
		#obtain firing rate trains
		x = bcs.readcolumn(bcs.filename_generator('Frate.txt', rp), to_float=True)
		x = [filter_data(x[i]) for i in range(1,len(x))]
		#state train
		state_train = list(zip(x[0],x[1],x[2],x[3])) #for 4 neurons only!
		state_train_trucks.append(state_train)
		state_train_all += state_train
		for i in range(4): state_neurons[i] += x[i]
	return state_train_all, state_train_trucks, state_neurons
	
def EP(state_train, nmax=4, plot=False):
	'''
	Performs the elbow test and returns relevant info about EP.
	returns: 1) int, number of EP, 2) list, coors of EP
			 3) list, state_train labeled into clusters, 4) float, ratio of SSE
	'''
	all_sse, all_mean, all_labels = [], [], [] #SSE represents how well the data performs
	for k in range(1,nmax+1):
		SSE=0
		labels = KMeans(n_clusters=k, max_iter=400, precompute_distances = 'auto').fit_predict(state_train)
		all_labels.append(labels)
		mean_train, temp = [[] for c in range(k)], []
		for s in range(len(labels)):
			mean_train[labels[s]].append(list(state_train[s]))
		for i in range(k):
			if mean_train[i]==[]: temp.append([0,0,0,0])
			else: temp.append(list(np.mean(np.array(mean_train[i]), axis=0)))
		all_mean.append(temp)
		for j in range(k):
			for ss in mean_train[j]:
				SSE += vect_len(ss, all_mean[k-1][j])**2
		all_sse.append(SSE)
	diff = [all_sse[i]-all_sse[i+1] for i in range(nmax-1)]
	ratios = []
	for i in range(nmax-2):
		if diff[i+1]!=0:ratios.append(diff[i]/float(diff[i+1]))
		else: ratios.append(np.inf)
	#plot
	if plot:
		plt.plot([i for i in range(1,nmax+1)], all_sse, 'k', alpha=0.5)
		plt.scatter([i for i in range(1,nmax+1)], all_sse, s=20, c='b')
		plt.xticks([i for i in range(1, nmax+1)])
		plt.xlabel('k (number of clusters)', fontsize=14)
		plt.ylabel('SSE (Hz squared)', fontsize=14)
		plt.title('Ratios: '+', '.join(bcs.to_string(ratios)))
		plt.savefig('elbow.png', dpi=100)
		plt.close('all')
	#determine if elbow exists
	for d in range(nmax-2):
		if much_larger_than(diff[d], diff[d+1], 7): return d+2, all_mean[d+1], all_labels[d+1], ratios[d]
	return 1, all_mean[0], all_labels[0], 0

def pca(fit_train, trans_train, dim=2):
	'''
	Performs pca on fit_train, and transforms trans_train using the results.
	returns: 1) list, new_trans_train, i.e. a list of coordinates, under pca transformation.
			 2) list, eigenvectors of the pca transformation.
	'''
	pca = IncrementalPCA(n_components=dim)
	try:
		pca.fit(fit_train)
	except:
		return [[0]*dim]*len(fit_train), 0
	comp = pca.components_
	new_trans_train = pca.transform(trans_train)
	return new_trans_train, comp 

def sort_clusters(nclusters, sorted_train, bin_size=0.01):
	'''
	returns: storage, which stores the attractors as garages
			 and len(time periods) that the system stays in said attractors items.
	'''	
	#count transitions and store
	STRG = storage([i for i in range(nclusters)], ['list' for i in range(nclusters)])
	previous_n, previous_name=0, sorted_train[0]
	for n in range(len(sorted_train)):
		name = sorted_train[n]
		if name != previous_name:
			STRG.store_list(previous_name, (n-previous_n)*bin_size)
			previous_name=name
			previous_n=n
	return STRG

def escape_rate(nclusters, STRG):
	'''
	returns: dict, with attractors as keys and its escape rate as values
	'''
	escape_dic={}
	for nc in range(nclusters):
		try:
			escape_dic[nc]=np.mean(STRG.retrieve_garage(nc))
		except RuntimeWarning:
			escape_dic[nc]=0
	return escape_dic

def get_fano(state_neurons):
	'''
	returns: the fano factor averaged the trials for each neuron.
	'''
	res=[]
	for neuron in state_neurons:
		res.append(bcs.fano(neuron))
	return res

def print_commands(plot, plot_fr, plot_pca, plot_elbow, trans, runnum, only):
	print('CAUTION: combs must be adjusted manually.')
	print('***Just to be clear, these are the commands for this run:')
	print('The trial you are running is: run'+runnum)
	print('Trans code is: '+trans+'\n')
	if plot and plot_fr: print('The firing rate graph will be plotted.')
	if plot and plot_pca: print('The pca graph will be plotted.')
	if plot and plot_elbow: print('The elbow graph will be plotted.')
	if not plot: print('Nothing will be plotted.')
	print('These actions will be done: '+', '.join(only))
		
def confidence(FFlist, ratio, esc_rate, harsh_criteria=10):
	deter1=lower_than(FFlist, harsh_criteria) #whether FF is VERY low
	#deter2=lower_than([float(ratio)], 15) #whether ratio is <15
	deter3=lower_than(bcs.to_float(esc_rate.split('_')),0.05) #whether all escape rate is < 0.05
	if deter1 and deter3: return '90_FP'
	elif deter1 and (not deter3): return '70_FP'
	else: return '30_FP'
	
def continuity(lists):
	res=[]
	for l in lists:
		previous_i, count, accumulate = 0, 1, []
		for j in range(len(l)):
			i=l[j]
			if j==(len(l)-2):
				accumulate.append(count)
			elif i !=0: count+=1
			elif previous_i!=0 and i==0:
				accumulate.append(count)
				count=1
			previous_i=i
		if accumulate==[]: res.append(0)
		else: res.append(np.mean(accumulate))
	return res

def correlation(lists):
	res=[]
	combs=combinations('0123',2)
	for comb in combs:
		A, B = lists[int(comb[0])], lists[int(comb[1])]
		try:
			corr, pval = pearsonr(A,B)
		except RuntimeWarning:
			corr, pval = 0, 0
		to_append='_'.join([comb[0]+comb[1],str(corr),str(pval)])
		res.append(to_append)
	return res

def determine_corr(num, FFlist, comb):
	FF1, FF2 = FFlist[int(comb[0])], FFlist[int(comb[1])]
	if FF1==0 and FF2==0: return 'none'
	elif num > 0.5: return 'pos'
	elif num < -0.5: return 'neg'
	elif FF1<5 and FF2<5 and num < 0.5: return 'pos'
	elif FF1<5 and FF2<5 and num > -0.5: return 'neg'
	else: return 'none'

def determine_states(runnum, motherpath=os.getcwd()):
	b_res=clss.branch('results_'+str(runnum), motherpath)
	corrlink=os.path.join(b_res.pathlink,'corr.txt')
	infolink=os.path.join(b_res.pathlink,'info.txt')
	info_entry, corr_entry=clss.entry([' 0', ' 0_6'], [' 1', ' 4']), clss.entry([' 0', ' 0_6'], [])
	infodata = info_entry.readdata_and_fix(infolink)
	corrdata = corr_entry.readdata_and_fix(corrlink, quick_format=True)
	motdic = {} #{motif: [[numEP, (comb, relation),..],..]}
	for key in corrdata.keys():
		#arrange each ps into [numEP, (comb, relation),..]
		numEP, motif, FF = infodata[key][0], key[1], bcs.to_float(infodata[key][1].split('_'))
		if motif not in motdic.keys(): motdic[motif]=[]
		temp = [numEP]
		for val in corrdata[key]:
			comb, crvalue, pvalue = val.split('_')
			relation = determine_corr(float(crvalue), FF, comb)
			temp.append((comb, relation))
		#Try to catch errors in states:
		relations = [combo[1] for combo in temp[1:]]
		if relations.count('none')>=3: temp=[numEP,('01','none'),('02','none'),('03','none'),\
											('12','none'),('13','none'),('23','none')]
		if relations.count('pos')>=3: temp=[numEP,('01','pos'),('02','pos'),('03','pos'),\
											('12','pos'),('13','pos'),('23','pos')]
		#Determine if there is already a qualitatively similar parameter set in the motif
		to_append = True
		for pms in motdic[motif]: #[[numEP, (comb, relation),..],..]
			exit = True
			for item in temp: #[numEP, (comb, relation),..]
				if item not in pms:
					exit = False
					break
			if exit:
				to_append = False
				break
		if to_append: motdic[motif].append(temp)
	return motdic
	
#determine_states('3')


def main(runnum, plot=False, outfile='info.txt', trans='ffffssff', motherpath=os.getcwd(), remedy=False, **kwargs):
	'''
	1) Plots the elbow and the frate graphs, moves all of them to a folder called 'graphs'.
	2) Performs kmeans on the data to get EP-relevant info, and:
		[1] Returns the pms along with their corresponding cooridnates of the EPs on 2-dimensions (determined by pca).
		[2] Plots the data of each param_set (PS) onto pca, labels the EPs, and moves it to folder 'graphs'.
		[3] Outputs PS, #EP, EP_coor, ratio, FF to file outfile.
	parameters:
	*plot: If any graph is plotted at all, it must be set to True.
	*only: Can select a few actions only.
	'''
	#kwargs
	kw={'plot_fr':False, 'plot_pca':False, 'plot_elbow':False, 'corrfile':'corr.txt', 'only':[], 'EPfile':'EPcoor.txt'}
	kw.update(kwargs)
	plot_fr, plot_pca, plot_elbow = kw['plot_fr'], kw['plot_pca'], kw['plot_elbow']
	corrfile, only, EPfile = kw['corrfile'], kw['only'], kw['EPfile']
	print_commands(plot, plot_fr, plot_pca, plot_elbow, trans, runnum, only)
	#dir
	runpath=os.path.join(os.getcwd(),'run'+runnum)
	os.chdir(runpath)
	alldirs=[dr for dr in os.listdir(runpath) if os.path.isdir(os.path.join(runpath,dr))]
	allpms=[]
	#deters
	deter_esc = (only==[] or ('esc' in only))
	deter_pca = (only==[] or ('pca' in only))
	deter_FF = (only==[] or ('FF' in only))
	deter_pw = (only==[] or ('pw' in only))
	deter_corr = (only==[] or ('corr' in only))
	deter_info = (only==[] or ('info' in only))
	deter_EP = (only==[] or ('EP' in only))
	#result files and outputs
	b_res=clss.branch('results_'+str(runnum), motherpath)
	b_graphs=clss.branch('graphs',b_res.pathlink)
	if os.path.exists(os.path.join(b_res.pathlink, outfile)):
		done=bcs.readcolumn(os.path.join(b_res.pathlink, outfile))[0]
	else:
		done=[]
	if plot: b_graphs.mkdir()
	else: b_res.mkdir()
	if deter_info and not remedy: bcs.output_clf(os.path.join(b_res.pathlink,outfile))
	if deter_corr and not remedy: bcs.output_clf(os.path.join(b_res.pathlink,corrfile))
	#analysis
	count=0
	for dr in alldirs:
		if dr not in done:
			count+=1
			print(str(count)+':'+dr)
			#specifications
			pm=parameter([],[])
			pm.assign_name(dr, trans_method=trans)
			b=clss.branch(dr,runpath)
			rp=10 if len(os.listdir(b.pathlink))>9 else 1 #due to flawed dataset
			#get EPs
			state_train_all, state_train_truck, state_neurons = get_state_train(b, runnum, repeat=rp)
			nclusters, EP_coors, label_train, ratio = EP(state_train_all, nmax=5, plot=(plot and plot_elbow))
			EP4 = ['_'.join(bcs.to_string(item)) for item in EP_coors]
			#calculate escape rate
			if deter_esc:
				accumulation=0
				STRG = storage([nc for nc in range(nclusters)], ['list' for i in range(nclusters)])
				for state_train in state_train_truck:
					lt=label_train[accumulation:len(state_train)+accumulation]
					accumulation+=len(state_train)
					new_STRG=sort_clusters(nclusters,lt)
					for nc in range(nclusters): STRG.massive[nc]+=new_STRG.massive[nc]
				ed = escape_rate(nclusters, STRG)
			#pcaPS
			if deter_pca:
				trans_train, comp = pca(state_train_all, state_train_all+EP_coors)
				x, y = [item[0] for item in trans_train], [item[1] for item in trans_train]
				pm.add_pair(('EP',EP_coors))
				allpms.append(pm)
				if plot and plot_pca:
					plt.plot(x[:-(nclusters)], y[:-(nclusters)], 'k.', alpha=0.5)
					plt.plot(x[-(nclusters):], y[-(nclusters):], 'b.')
					plt.xlabel('dim1', fontsize=14)
					plt.ylabel('dim2', fontsize=14)
					plt.savefig(dr+'_pcaPS.png', dpi=100)
					plt.close('all')
			#fano factor
			if deter_FF:
				FF=get_fano(state_neurons)
				FFall='_'.join(bcs.to_string(FF))
			#pulse width
			if deter_pw:
				pwidth=continuity(state_neurons)
				pall='_'.join(bcs.to_string(pwidth))
			#correlation
			if deter_corr:
				all_corrs=correlation(state_neurons)
			#move graphs and outputs
			if plot:
				if plot_elbow: os.rename('elbow.png',dr+'_elbow.png')
				if plot_fr: subplot(fname=os.path.join(b.pathlink, 'Frate.txt'), outputfigname=dr+'_frate.png', tstep=5,\
								   title='Fano Factors: '+FFall, tight=False, dpi=100)
				if plot_elbow: b_graphs.move_from(dr+'_elbow.png',b.pathlink)
				if plot_fr: b_graphs.move_from(dr+'_frate.png',b.pathlink)
				if plot_pca: b_graphs.move_from(dr+'_pcaPS.png',b.pathlink)
			if deter_info: 
				vals=[ed[key] for key in sorted(ed.keys())]
				numEP, esc_rate, ratio= str(len(ed.keys())), '_'.join(bcs.to_string(vals)), str(ratio)
				bcs.output_line(os.path.join(b_res.pathlink,outfile),\
								' '.join([dr,numEP, esc_rate,ratio,FFall,pall]))
			if deter_corr:
				bcs.output_line(os.path.join(b_res.pathlink,corrfile),\
								' '.join([dr]+all_corrs))
			if deter_EP:
				bcs.output_line(os.path.join(b_res.pathlink,EPfile),\
								' '.join([dr]+EP4))
	os.chdir(motherpath)
	return allpms
	
def plot_pcaEPs(runnum, plot=False, motherpath=os.getcwd(), trans='fffsfss', feature='i'):
	'''
	Plots the cooridnates of the EPs on 2-dimensions (determined by pca).
	'''
	b_res=clss.branch('results_'+str(runnum),motherpath)
	b_pcaEP=clss.branch('pcaEP',b_res.pathlink)
	b_pcaEP.mkdir()
	os.chdir(b_pcaEP.pathlink)
	#reevaluation
	newfile=os.path.join(b_res.pathlink, 'new-info.txt')
	EPfile=os.path.join(b_res.pathlink, 'EPcoor.txt')
	Ent=clss.entry(' 0', [' 6'])
	Ent2=clss.entry(' 0', [])
	data=Ent.readdata_and_fix(newfile)
	EPcoorss=Ent2.readdata_and_fix(EPfile, quick_format=True)
	#sort pms by motifs(ID)
	IDdic, coldic={}, {}
	for key in data.keys():
		pm=parameter([], [])
		pm.assign_name(key, trans_method='iiiiss')
		motID, col=pm.extract('ID'), pm.extract(feature)
		if motID not in IDdic.keys(): IDdic[motID]=[]
		if motID not in coldic.keys(): coldic[motID]=[]
		actual_coors=[bcs.to_float(cr.split('_')) for cr in EPcoorss[key]]
		#reevaluation
		if data[key][0]=='o':
			IDdic[motID]+=actual_coors
			coldic[motID]+=[int(col)]*len(actual_coors)
		else:
			IDdic[motID]+=[list(np.mean(actual_coors,axis=0))]
			coldic[motID]+=[int(col)]
	for motID in IDdic.keys():
		print(motID)
		EP_coors=IDdic[motID]
		#pca
		trans_train, vectors = pca(EP_coors, EP_coors)
		vec_strs = [str(round(vec[0],2))+' '+str(round(vec[1],2)) for vec in vectors]
		#elbow then pca
		nclusters, new_EP_coors, label_train, ratio = EP(EP_coors, nmax=6, plot=False)
		new_trans_train = pca(EP_coors, new_EP_coors)[0]
		if plot:
			plt.plot([item[0] for item in trans_train], [item[1] for item in trans_train], 'k.')
			plt.plot([item[0] for item in new_trans_train], [item[1] for item in new_trans_train], 'b.')
			plt.xlabel('dim1: '+vec_strs[0], fontsize=14)
			plt.ylabel('dim2: '+vec_strs[1], fontsize=14)
			plt.savefig('pcaEP_'+motID+'.png', dpi=200)
			plt.close('all')
			#try
			from dynalysis.data_visualization import plot_continuous_colorbar
			plot_continuous_colorbar([item[0] for item in trans_train], [item[1] for item in trans_train],\
									 coldic[motID], 'dim1: '+vec_strs[0], 'dim2: '+vec_strs[1], feature,\
									svf='pcaEP_'+motID+'_'+feature+'.png')
		bcs.output_clf('pcaEP_'+motID+'.txt')
		bcs.output_double('pcaEP_'+motID+'.txt', EP_coors)
	os.chdir(motherpath)
	return 0

def evaluation_by_FF(runnum, infile='info.txt', FF_criteria=30, cp=False):
	'''
	Evaluates the number of clusters based off ratio, Fano Factor, correlation and ON/OFF.
	parameters:
	*data: a dictionary with motifID as key and [numEP, esc_rate, ratio, FF1, FF2, FF3, FF4] as values
	'''
	#dir
	file=os.path.join(os.getcwd(),'results_'+runnum, infile)
	ofile=os.path.join(os.getcwd(),'results_'+runnum, 're-'+infile)
	nfile=os.path.join(os.getcwd(),'results_'+runnum, 'new-'+infile)
	b_res_path=os.path.join(os.getcwd(),'results_'+runnum)
	b_graphs=clss.branch('graphs', b_res_path)
	if cp:
		b_regraphs_s=clss.branch('re-graphs-s', b_res_path)
		b_regraphs_n=clss.branch('re-graphs-n', b_res_path)
		b_regraphs_s.mkdir()
		b_regraphs_n.mkdir()
	bcs.output_clf(ofile)
	bcs.output_clf(nfile)
	Ent=clss.entry(' 0', [' 1', ' 2', ' 3', ' 4_0', ' 4_1', ' 4_2', ' 4_3', ' 5_0', ' 5_1', ' 5_2', ' 5_3'])
	data=Ent.readdata_and_fix(file)
	#main
	new_data={}
	for key in data:
		numEP, esc_rate, ratio, FF1, FF2, FF3, FF4, pw1, pw2, pw3, pw4 = data[key]
		FFlist=bcs.to_float([FF1, FF2, FF3, FF4])
		FFstring='_'.join([FF1, FF2, FF3, FF4])
		pwlist=bcs.to_float([pw1, pw2, pw3, pw4])
		pwstring='_'.join([pw1, pw2, pw3, pw4])
		deter1=(int(numEP)>1) #can only have FPs if numEP>1 by definition
		deter2=lower_than(FFlist, FF_criteria) #if FF is too low
		deter3=lower_than(bcs.to_float(esc_rate.split('_')),0.1) #whether all escape rate is < 0.1
		deter4=lower_than(pwlist,5) #whether all pulse width < 5
		if deter1 and deter2 and deter3: #identify False Positives type-s (saturation)
			conf=confidence(FFlist, ratio, esc_rate)
			new_data[key]=['1', '0', '0', FFstring, pwstring, conf]
			bcs.output_line(ofile,' '.join([key]+new_data[key]))
			if cp: b_graphs.cp_to(key+'_frate.png', b_regraphs_s.pathlink)
		elif deter1 and deter4: #identify False Positives type-n (noise-induced firing)
			new_data[key]=['1', '0', '0', FFstring, pwstring, 'n']
			bcs.output_line(ofile,' '.join([key]+new_data[key]))
			if cp: b_graphs.cp_to(key+'_frate.png', b_regraphs_n.pathlink)
		else: #False Negatives not implemented
			new_data[key]=[numEP, esc_rate, ratio, FFstring, pwstring, 'o'] #correct trials with confidence='o'
		bcs.output_line(nfile,' '.join([key]+new_data[key]))
	return data

def evaluation_by_FF_only(runnum, infile='info.txt', FF_criteria=30, cp=False):
	'''
	Evaluates the number of clusters based off Fano Factor.
	parameters:
	*data: a dictionary with motifID as key and [numEP, esc_rate, ratio, FF1, FF2, FF3, FF4] as values
	'''
	#dir
	file=os.path.join(os.getcwd(),'results_'+runnum, infile)
	ofile=os.path.join(os.getcwd(),'results_'+runnum, 're-'+infile)
	nfile=os.path.join(os.getcwd(),'results_'+runnum, 'new-'+infile)
	b_res_path=os.path.join(os.getcwd(),'results_'+runnum)
	b_graphs=clss.branch('graphs', b_res_path)
	if cp:
		b_regraphs_s=clss.branch('re-graphs-s', b_res_path)
		b_regraphs_n=clss.branch('re-graphs-n', b_res_path)
		b_regraphs_s.mkdir()
		b_regraphs_n.mkdir()
	bcs.output_clf(ofile)
	bcs.output_clf(nfile)
	Ent=clss.entry(' 0', [' 1', ' 2', ' 3', ' 4_0', ' 4_1', ' 4_2', ' 4_3', ' 5_0', ' 5_1', ' 5_2', ' 5_3'])
	data=Ent.readdata_and_fix(file)
	#main
	new_data={}
	for key in data:
		numEP, esc_rate, ratio, FF1, FF2, FF3, FF4, pw1, pw2, pw3, pw4 = data[key]
		FFlist=bcs.to_float([FF1, FF2, FF3, FF4])
		FFstring='_'.join([FF1, FF2, FF3, FF4])
		pwlist=bcs.to_float([pw1, pw2, pw3, pw4])
		pwstring='_'.join([pw1, pw2, pw3, pw4])
		deter1=(int(numEP)>1) #can only have FPs if numEP>1 by definition
		deter2=lower_than(FFlist, FF_criteria) #if FF is too low
		if deter1 and deter2: #identify False Positives type-s (saturation)
			conf=confidence(FFlist, ratio, esc_rate)
			new_data[key]=['1', '0', '0', FFstring, pwstring, conf]
			bcs.output_line(ofile,' '.join([key]+new_data[key]))
			if cp: b_graphs.cp_to(key+'_frate.png', b_regraphs_s.pathlink)
		else: #False Negatives not implemented
			new_data[key]=[numEP, esc_rate, ratio, FFstring, pwstring, 'o'] #correct trials with confidence='o'
		bcs.output_line(nfile,' '.join([key]+new_data[key]))
	return data

def evaluation_by_mean(runnum, infile='new-info.txt', outfile='motif-numEP.txt'):
	#dir
	file=os.path.join(os.getcwd(),'results_'+runnum, infile)
	ofile=os.path.join(os.getcwd(),'results_'+runnum, outfile)
	bcs.output_clf(ofile)
	Ent=clss.entry(' 0_6', [' 1', ' 4_0', ' 4_1', ' 4_2', ' 4_3'])
	data=Ent.readdata(file)
	for key in data.keys():
		all_numEP = [int(num[0]) for num in data[key]]
		tot = np.sum(all_numEP)
		avg = np.mean(all_numEP)
		c1, c2, c3 = all_numEP.count(1), all_numEP.count(2), all_numEP.count(3)
		adj_avg_list, c0=[], 0
		for num in data[key]:
			if np.sum([float(FF) for FF in num[1:]])!=0: adj_avg_list.append(int(num[0]))
			else: c0+=1
		adj_avg = np.mean(adj_avg_list)
		to_output = [key, str(tot), str(avg), str(adj_avg), str(c0), str(c1-c0), str(c2), str(c3)]
		bcs.output_line(ofile, ' '.join(to_output))
	return 0

def evaluation_by_mean_set(setnum, infile='new-info.txt', ofile='motif-numEP', allnum=['1']):
	alls=np.load('set'+str(setnum)+'.npy')
	ofile=os.path.join(os.getcwd(), ofile+'-set'+str(setnum)+'.txt')
	bcs.output_clf(ofile)
	Ent, data=clss.entry(' 0', [' 1', ' 4_0', ' 4_1', ' 4_2', ' 4_3', ' 0_5']), {}
	for runnum in allnum:
		file=os.path.join(os.getcwd(),'results_'+runnum, infile)
		data.update(Ent.readdata(file))
	motifdic={}
	for sett in alls:
		e,i,n,comb,amp,bi,mot,ie,pb=sett ##
		pm=parameter([e,i,n,amp,bi,str(mot.typeID)+'-'+str(mot.ID),ie,pb],['e','i','noise','A','Bi','ID','gIE','pb']) ##
		if data[pm.name][0][5] not in motifdic: motifdic[data[pm.name][0][5]]=[]
		motifdic[data[pm.name][0][5]].append(data[pm.name][0][:5])
	for key in motifdic:
		all_numEP = [int(num[0]) for num in motifdic[key]]
		tot = np.sum(all_numEP)
		avg = np.mean(all_numEP)
		c1, c2, c3, c4 = all_numEP.count(1), all_numEP.count(2), all_numEP.count(3), all_numEP.count(4)
		adj_avg_list, c0=[], 0
		for num in motifdic[key]:
			if np.sum([float(FF) for FF in num[1:]])!=0: adj_avg_list.append(int(num[0]))
			else: c0+=1
		adj_avg = np.mean(adj_avg_list)
		to_output = [key, str(tot), str(avg), str(adj_avg), str(c0), str(c1-c0), str(c2), str(c3), str(c4)]
		bcs.output_line(ofile, ' '.join(to_output))
		
def evaluation_by_ratioSSE(runnum, criteria, infile='info.txt', cp=False):
	#dir
	file=os.path.join(os.getcwd(),'results_'+runnum, infile)
	nfile=os.path.join(os.getcwd(),'results_'+runnum, 'new2-'+infile)
	b_res_path=os.path.join(os.getcwd(),'results_'+runnum)
	b_graphs=clss.branch('graphs', b_res_path)
	if cp:
		b_regraphs_R=clss.branch('re-graphs-R', b_res_path)
		b_regraphs_R.mkdir()
	bcs.output_clf(nfile)
	Ent=clss.entry(' 0', [' 3', ' 1'])
	data=Ent.readdata_and_fix(file)
	Ent=clss.entry(' 0', [])
	alldata=Ent.readdata_and_fix(file, quick_format=True)
	#main
	for key in data:
		ratios = bcs.to_float(data[key][0].split('_'))
		passed=[]
		for r in range(len(ratios)):
			if ratios[r] > criteria: passed.append(r)
		if alldata[key][0] != str(r+1) and cp:
			b_graphs.cp_to(key+'_frate.png', b_regraphs_R.pathlink)
		if len(passed)==1: bcs.output_line(nfile,' '.join([key, str(r+1)]+alldata[key][1:]))
		else: bcs.output_line(nfile,' '.join([key, str(1)]+alldata[key][1:]))
	return data
	
if __name__=='__main__':
	pass
	runnum='31'
	pms=main(runnum, plot=True, plot_fr=True, plot_pca=True, plot_elbow=True, remedy=True)
	evaluation_by_ratioSSE(runnum, 7, infile='info.txt', cp=True)
	data=evaluation_by_FF_only(runnum, cp=True, FF_criteria=37)
	plot_pcaEPs(runnum, plot=True, feature='i')
	setnum='31'
	evaluation_by_mean_set(setnum, allnum=[setnum])
	
	