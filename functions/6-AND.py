import os
import subprocess
from dynalysis import basics as bcs
from dynalysis import classes as clss
from dynalysis.gen_conf import exec_conf
from dynalysis.gen_pro import exec_pro
from dynalysis.visualization import subplot, edplot

#=====parameters=====#
e1i1=30
e1i2=0
e2i1=0
e2i2=30
i1e1=60
i1e2=0
i1i2=30
i2e1=0
i2e2=60
i2i1=30
ee5=50
input3=45
input4=45
duration=1000
bias_current=-.65
eee= -.2
iii= -.6

exec_conf([0,30,0,0], ID=4095, update_tar={\
          'exc_1':[['inh_1','AMPA',e1i1,1,1], ['inh_2','AMPA',e1i2,1,1], ['exc_5','AMPA',ee5,1,1]],\
          'exc_2':[['inh_1','AMPA',e2i1,1,1], ['inh_2','AMPA',e2i2,1,1], ['exc_5','AMPA',ee5,1,1]],\
          'inh_1':[['exc_1','GABA',i1e1,1,1], ['exc_2','GABA',i1e2,1,1], ['inh_2','GABA',i1i2,1,1]],\
          'inh_2':[['exc_1','GABA',i2e1,1,1], ['exc_2','GABA',i2e2,1,1], ['inh_1','GABA',i2i1,1,1]],\
          'exc_3':[['exc_1','AMPA',input3,1,1], ['inh_1','AMPA',input3,1,1]],\
          'exc_4':[['exc_2','AMPA',input4,1,1], ['inh_2','AMPA',input4,1,1]]},\
          event_dic={'exc_1':{'Taum':20}, 'exc_2':{'Taum':20},\
                     'inh_1':{'Taum':20}, 'inh_2':{'Taum':20}},\
          add_connections=[('exc_3','exc_1'),('exc_3', 'inh_1'),\
                           ('exc_4','exc_2'),('exc_4', 'inh_2'),\
                           ('exc_1','exc_5'),('exc_2','exc_5')])

exec_pro(ID=4095, trial_t=5500, output_Spike=False, gmean=0, gstd=0,\
         event_dic={(100, 'exc_1'):['pulse_b',duration,1+eee,0,eee],\
                    (100, 'exc_2'):['pulse_b',duration,1+eee,0,eee],\
                    (2105, 'exc_1'):['pulse_b',duration,1+eee,0,eee],\
                    (2100, 'exc_2'):['pulse_b',duration,1+eee,0,eee],\
                    (4100, 'exc_1'):['pulse_b',duration,1+eee,0,eee],\
                    (4105, 'exc_2'):['pulse_b',duration,1+eee,0,eee],\
                    (100, 'inh_1'):['pulse_b',duration,1+iii,0,eee],\
                    (100, 'inh_2'):['pulse_b',duration,1+iii,0,eee],\
                    (2105, 'inh_1'):['pulse_b',duration,1+iii,0,eee],\
                    (2100, 'inh_2'):['pulse_b',duration,1+iii,0,eee],\
                    (4100, 'inh_1'):['pulse_b',duration,1+iii,0,eee],\
                    (4105, 'inh_2'):['pulse_b',duration,1+iii,0,eee],\
                    (2, 'exc_1'):['ramp',eee,0],\
                    (2, 'exc_2'):['ramp',eee,0],\
                    (2, 'inh_1'):['ramp',iii,0],\
                    (2, 'inh_2'):['ramp',iii,0],\
                    (2, 'exc_5'):['ramp',bias_current,0]},\
        output_MemPot=True, output_Frate=False, MemPot='6-AND.dat')

subprocess.call(['./ss_flysim_iter=0.sh'])

#edplot(ed) #plots the protocol
#subplot(fname='6-AND.dat', cols=[5,6,7], col_labels=['output','inh_1','inh_2']) #plots the results

