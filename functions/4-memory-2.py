import subprocess
import dynalysis.basics as bcs
import dynalysis.classes as clss
from dynalysis.gen_conf import exec_conf
from dynalysis.gen_pro import exec_pro
from dynalysis.visualization import subplot, edplot

ut={}
ed={(0, 'exc_1'):['ramp',0.5,0],\
    (0, 'exc_2'):['ramp',0.5,0],\
    (0, 'inh_1'):['ramp',0.5,0],\
    (0, 'inh_2'):['ramp',0.5,0],\
    (1000, 'exc_1'):['pulse_b',30,0.9,0,0.5],\
    (1000, 'exc_2'):['pulse_b',30,0.9,0,0.5],\
    (1000, 'inh_1'):['pulse_b',30,0.9,0,0.5],\
    (1000, 'inh_2'):['pulse_b',30,0.9,0,0.5],\
    (2000, 'exc_1'):['pulse_b',30,-1,0,0.5],\
    (2000, 'exc_2'):['pulse_b',30,-1,0,0.5],\
    (2000, 'inh_1'):['pulse_b',30,-1,0,0.5],\
    (2000, 'inh_2'):['pulse_b',30,-1,0,0.5]} 
exec_conf([45,75,2,30], ID=3435)

exec_pro(ID=3435,output_Spike=False, trial_t=3000, output_MemPot=True, event_dic=ed, MemPot='4-memory-2.dat',\
		output_Frate = False)
subprocess.call(['./ss_flysim_iter=0.sh'])

#edplot(ed) #plots the protocol
#subplot(fname='4-memory-2.dat') #plots the results