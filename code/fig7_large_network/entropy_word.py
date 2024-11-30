import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import dynalysis.basics as bcs
from PRetina.MIv3 import mi_quick

folder = 'gee=10_pii=0.5 (tnum=100)(memory v2)'

matrix = []
for gii in range(0,30,1):
    temp = []
    for gie in range(0,30,2):
        file = os.path.join(folder,'gie='+str(gie)+'_gii='+str(gii),'Frate.txt')
        data = bcs.readcolumn(file); print('gie='+str(gie)+'_gii='+str(gii))
        dic = {}; c = 0
        for tn in range(100):
            word = []
            for neu in range(1,len(data)): #100-150: s, 150-250: m
                mean = np.mean(bcs.to_float(data[neu][150+tn*300:251+tn*300]))
                state = int(round(mean//40,0))
                word.append(state)
            word = tuple(word)
            if word not in dic.keys(): dic[word] = c; c+= 1
        MI = mi_quick(sorted(dic.values()),sorted(dic.values()),0,(list(range(0,c)),list(range(0,c))))
        temp.append(MI) ##
    matrix.append(temp)
print(matrix)
    
ax = sns.heatmap(matrix, cmap='YlGnBu')
ax.set_xticks(np.arange(0,15+5,5))
ax.set_yticks(list(range(0,30+5,5)))
ax.set_xticklabels(list(range(0,30+10,10)))
ax.set_yticklabels(list(range(0,30+5,5)))
ax.invert_yaxis()
plt.xlabel('gie'); plt.ylabel('gii')
plt.savefig(folder+'_MI3 (bi, quo=40).png', dpi=300)
