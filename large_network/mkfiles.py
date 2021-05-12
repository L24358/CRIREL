import os
import numpy as np
import dynalysis.classes as clss
from dynalysis.gen_pro import exec_pro
from dynalysis.gen_conf import exec_conf
from dynalysis.classes2 import Lokis_lair
from dynalysis.visualization import edplot

N = 100 #number of neurons
perc = [0.5,0.5,0.5,0.5] #percentage of gie, gee, gie, gii
Lokes = Lokis_lair(N,0.25,perc) #0.25: percentage of inhibitory neurons
Freyja = os.path.join(os.getcwd(),'memory_v3 (pii=0.5)')

ed={}; tnum=100; bias=-5 #tnum=100
for neu in Lokes.nodes: ed[(1,neu)]=['ramp',bias,0]
for tn in range(tnum):
	for neu in Lokes.nodes: #random pulse to each neuron
		amp = np.random.normal(3,3)
		ed[(1002+tn*3000,neu)]=['pulse_b',500,bias+amp*np.sign(amp),0,bias]
	for neu in Lokes.nodes[:75]: #only exc are reset
		ed[(1002+tn*3000+1500,neu)]=['pulse_b',500,-30,0,bias]
edplot(ed, keys=['exc_1','exc_2','inh_1'],trial_t=3000*tnum+1002)

for gie in range(0,30,2):
	for gii in range(0,30,1):
		b = clss.branch('gie='+str(gie)+'_gii='+str(gii), Freyja); b.mkdir(); os.chdir(b.pathlink)
		exec_conf([6,10,gie,gii], Mot=Lokes, suppress=True)
		exec_pro(Mot=Lokes,trial_t=3000*tnum+1002, output_Spike=False, event_dic=ed, baseline=True,\
				gmean=bias,gstd=0)
		clss.copy_flysim(b.pathlink)

