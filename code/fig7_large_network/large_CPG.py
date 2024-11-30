import os
import numpy as np
import dynalysis.classes as clss
from dynalysis.gen_pro import exec_pro
from dynalysis.gen_conf import exec_conf
from dynalysis.classes2 import Lokis_lair
from dynalysis.visualization import edplot

N = 100 #number of neurons
perc = [0.5,0.5,0.5,1] #percentage of gie, gee, gie, gii
Lokes = Lokis_lair(N,0.25,perc)
Freyja = os.getcwd()

ed={}; biasi=0; biase=1
for neu in Lokes.nodes[:75]: ed[(1,neu)]=['ramp',biase,0]
for neu in Lokes.nodes[75:]: ed[(1,neu)]=['ramp',biasi,0]
edplot(ed, keys=['exc_1','exc_2','inh_1'])

gie = 32
gii = 40
gee = 20
gei = 32
exec_conf([gei,gee,gie,gii], Mot=Lokes, suppress=True)
exec_pro(Mot=Lokes,trial_t=10000, output_Spike=False, event_dic=ed, baseline=False, output_MemPot=True)
clss.copy_flysim(Freyja)

