'''===========================================
1. Adjust the parameters
2. Adjust rnum, which will become the name of the file
==========================================='''

import os
import numpy as np
import dynalysis.basics as bcs
import dynalysis.classes as clss
from itertools import combinations
from dynalysis.gen_pro import exec_pro
from dynalysis.gen_conf import exec_conf
from dynalysis.visualization import edplot
from sub_basics import parameter

#parameters
elist, ilist=[90,105,120,135,150], [30,40,50,60,70,80,90]
noise=[0] #since it's a brief current pulse
repeat=1
amplitude=[1,3,5,7]
bistable=['F']
motif=[clss.motif(4086,listtype='0011')]
combs=[p for p in combinations([0,1,2,3],1)]+[p for p in combinations([0,1,2,3],2)]+\
        [p for p in combinations([0,1,2,3],3)]+[p for p in combinations([0,1,2,3],4)]
gIE = [15,30,45]
pblist = [0]

def record_runs(runnum):
	b=clss.branch('run'+str(runnum),os.getcwd())
	alldirs=[dr for dr in os.listdir(b.pathlink) if os.path.isdir(os.path.join(b.pathlink, dr))]
	for dr in alldirs:
		bcs.output_line('parameters.txt',dr)
		
def in_record(name, motherpath):
	path=os.path.join(motherpath, 'parameters.txt')
	data=[row.split('\n')[0] for row in bcs.readline(path)]
	if name in data: return True
	return False

def export_set(num):
	alls=bcs.param_files([elist,ilist,noise,combs,amplitude,bistable,motif,gIE,pblist]) ##
	np.save('set'+str(num)+'.npy', alls)

def mkfiles(runnum, partial=False, exclude=[]):
	'''Bi is not implemented yet'''
	simpath=os.getcwd()
	motherpath=clss.branch('run'+str(runnum),os.getcwd()).pathlink
	done=clss.branch('done'+str(runnum),os.getcwd())
	if not partial:
		clss.branch('run'+str(runnum),os.getcwd()).mkdir()
		clss.copy_flysim(motherpath)
	else:
		print(done.pathlink)
		done.mkdir()
	#exclude
	allsets=[elist,ilist,noise,amplitude,bistable,motif,gIE,pblist]
	param_all=bcs.param_files(allsets)
	print(len(param_all))
	#main
	plot_ed = True
	for sett in param_all:
		e,i,n,amp,bi,mot,ie,pb=sett
		ed, pm={}, parameter([e,i,n,amp,bi,str(mot.typeID)+'-'+str(mot.ID),ie,pb],['e','i','noise','A','Bi','ID','gIE','pb'])
		if not in_record(pm.name, simpath):
			count = 0
			for neu in mot.neulist: ed[(2, neu)]=['ramp',pb,0] #pulse_b
			for comb in combs:
				for c in comb: 
					ed[(3+1000*count,mot.neulist[c])]=['pulse_b',100,amp,n,pb] #pulse_b
					for neu in mot.neulist: ed[(900+1000*count,neu)]=['pulse_b',100,-3,0,pb] #pulse_b
				count += 1
			b=clss.branch(pm.name,motherpath)
			if not partial: b.mkdir()
			if (partial==False) or ('Frate.txt_'+str(repeat) not in os.listdir(b.pathlink)):
				os.chdir(b.pathlink)
				exec_conf([30,e,ie,i], Mot=mot, suppress=True)
				exec_pro(Mot=mot,trial_t=1000*count, output_Spike=False, event_dic=ed,\
						 gmean=0, gstd=n, baseline=True)
				if plot_ed: edplot(ed); plot_ed = False; print(pm.name)
			elif partial:
				print('skipped: '+pm.name)
				b.move_self(done.pathlink)

#run
rnum = '31'
export_set(rnum)
mkfiles(rnum)

