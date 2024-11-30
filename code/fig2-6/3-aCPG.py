import subprocess
import dynalysis.basics as bcs
import dynalysis.classes as clss
from dynalysis.gen_conf import exec_conf
from dynalysis.gen_pro import exec_pro
from dynalysis.visualization import subplot, edplot

ut={}
ed={(0, 'inh_1'):['pulse',15,3,0],\
    (2, 'exc_1'):['ramp',2.4,0],\
    (2, 'exc_2'):['ramp',2.4,0],\
    (2, 'inh_1'):['ramp',0,0],\
    (2, 'inh_2'):['ramp',0,0]}
exec_conf([70,10,70,25], ID=3435, GABA=['GABA',40,-90,0,0,0])
exec_pro(ID=3435,output_Spike=False, trial_t=1000, output_MemPot=True, event_dic=ed, MemPot='3-aCPG.dat',\
		output_Frate = False)
subprocess.call(['./ss_flysim_iter=0.sh'])

#edplot(ed) #plots the protocol
#subplot(fname='3-aCPG.dat') #plots the results