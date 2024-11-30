import subprocess
import dynalysis.basics as bcs
import dynalysis.classes as clss
from dynalysis.gen_conf import exec_conf
from dynalysis.gen_pro import exec_pro
from dynalysis.visualization import subplot, edplot

ut={}
ed={(100, 'exc_1'):['pulse',50,0.5,0],\
    (100, 'exc_2'):['pulse',50,0.5,0],\
    (400, 'exc_1'):['pulse',50,-0.5,0],\
    (400, 'exc_2'):['pulse',50,-0.5,0]}
exec_conf([0,110,0,0], ID=3435, update_tar=ut)
exec_pro(ID=3435,output_Spike=False, trial_t=1000, output_MemPot=True, event_dic=ed, MemPot='1-bistability.dat',\
		output_Frate = False)
subprocess.call(['./ss_flysim_iter=0.sh'])

#edplot(ed) #plots the protocol
#subplot(fname='1-bistability.dat') #plots the results