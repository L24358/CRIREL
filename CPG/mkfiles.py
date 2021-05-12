import os
import numpy as np
import dynalysis.classes as clss
from gen_pro import exec_pro
from dynalysis.gen_conf import exec_conf

mother = os.getcwd()
for ID in [1131,3426,2347,3401,274,1290,3435,1791,4086]: #1131,3426,2347,3401,274,1290,3435,1791,4086
	Freyja = os.path.join(mother,'scan_'+str(ID))
	e1b=0.5; e2b=0.5
	ed={(0, 'exc_1'):['ramp',e1b+3,0],\
		(0, 'exc_2'):['ramp',e2b+3,0],\
		(0, 'inh_1'):['ramp',-0.5,0],\
		(0, 'inh_2'):['ramp',-0.5,0]}
	
	for gei in range(0,60,5): #65
		for gie in range(0,200,10): #200
			b = clss.branch('gei='+str(gei)+'_gie='+str(gie), Freyja); b.mkdir(); os.chdir(b.pathlink)
			exec_conf([gei,65,gie,5], ID=ID)
			exec_pro(ID=ID,output_Spike=True, trial_t=2000, event_dic=ed,\
				output_Frate = True, output_MemPot = True)
			clss.copy_flysim(b.pathlink)

