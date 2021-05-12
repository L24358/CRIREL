import subprocess
import dynalysis.basics as bcs
import dynalysis.classes as clss
from dynalysis.gen_conf import exec_conf
from dynalysis.gen_pro import exec_pro
from dynalysis.visualization import subplot, edplot

e1b=0.5
e2b=0.5

ut={}
ed={(0, 'exc_1'):['ramp',e1b+3,0],\
    (0, 'exc_2'):['ramp',e2b+3,0],\
    (0, 'inh_1'):['ramp',-0.5,0],\
    (0, 'inh_2'):['ramp',-0.5,0]}
exec_conf([30,65,200,5], ID=3435, update_tar=ut)
exec_pro(ID=3435,output_Spike=False, trial_t=1000, output_MemPot=True, event_dic=ed, MemPot='2-sCPG.dat',\
		output_Frate = False)
subprocess.call(['./ss_flysim_iter=0.sh'])

#edplot(ed) #plots the protocol
#subplot(fname='2-sCPG.dat') #plots the results