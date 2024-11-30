import subprocess
import dynalysis.basics as bcs
import dynalysis.classes as clss
from dynalysis.gen_conf import exec_conf
from dynalysis.gen_pro import exec_pro
from dynalysis.visualization import subplot, edplot

p=-5
dur=100
i=2
e=0
ut={}
ed={(2, 'inh_1'):['pulse_b',15,5,0,i],\
    (200,'inh_1'):['pulse_b',+dur,i+p,0,+i],\
    (200,'inh_2'):['pulse_b',+dur,i+p,0,i],\
    (500,'inh_1'):['pulse_b',+dur,i+p,0,+i],\
    (500,'inh_2'):['pulse_b',+dur,i+p,0,i],\
    (800,'inh_2'):['pulse_b',+dur,i+p,0,+i],\
    (800,'inh_1'):['pulse_b',+dur,i+p,0,i],\
    (0, 'exc_1'):['ramp',2.4+e,0],\
    (0, 'exc_2'):['ramp',2.4+e,0],\
    (0, 'inh_1'):['ramp',i,0],\
    (0, 'inh_2'):['ramp',i,0]}
exec_conf([70,10,70,25], ID=3435, GABA=['GABA',40,-90,0,0,0])
exec_pro(ID=3435,output_Spike=False, trial_t=1300, output_MemPot=True, event_dic=ed, MemPot='9-toggle-decide.dat',\
		output_Frate = False)
subprocess.call(['./ss_flysim_iter=0.sh'])

#edplot(ed) #plots the protocol
#subplot(fname='9-toggle-decide.dat') #plots the results