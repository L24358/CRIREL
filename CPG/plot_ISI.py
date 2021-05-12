import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
import dynalysis.basics as bcs
import dynalysis.classes as clss

def myround(x, base=0.01): return base * round(x/base)

def ISI(train):
	res = []
	for t in range(len(train)-1): res.append(myround(train[t+1]-train[t]))
	return res

mother = os.getcwd()
for ID in [1131,3426,2347,3401,274,1290,3435,1791,4086]:
	Freyja = os.path.join(mother,'all_scan','scan_'+str(ID))
	b_graph = clss.branch('scan_'+str(ID)+'_mempot', os.getcwd()); b_graph.mkdir()

	count = 0
	for gei in range(0,60,5): #65
		for gie in range(0,200,10): #200
			b = clss.branch('gei='+str(gei)+'_gie='+str(gie), Freyja)
			fname = os.path.join(b.pathlink, 'Spike.txt')
			nlist = bcs.spiketrain(fname)
			spktrain = [t for t in nlist[0] if t >= 0.5]
			ISItrain = ISI(spktrain)
			
			#fig = plt.figure(figsize=(10,5))
			#ax1 = fig.add_subplot(121)
			#fname = os.path.join(b.pathlink, 'MemPot.dat')
			#data = bcs.readcolumn(fname)
			#for neu in range(1,5):
			#	y = bcs.to_float(data[neu])[5:1005]
			#	ax1.plot(y, label='neu = '+str(neu+1))
			#plt.legend()
			
			#ax2 = fig.add_subplot(122)
			#plt.hist(ISItrain)
			
			#plt.savefig(os.path.join(b_graph.pathlink, 'gei='+str(gei)+'_gie='+str(gie)+'.png'), dpi=300)
			#plt.close('all')
			
			if len(list(set(ISItrain))) > 1: count += 1
	
	bcs.output_line(os.path.join(mother,'counts.txt'), 'ID='+str(ID)+': '+str(count))
			
			