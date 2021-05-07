'''===========================================
Plotting bar chart for number of EPs for different motifs.
==========================================='''

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

x=np.arange(7)
c1=[1290,11760+6300,6450+6300,2475,6015,18900,18900] #c0+c1
c2=[3960,840,0,2660,6645,0,0]
c3=[12645,0,6150,10035,6000,0,0]
c4=[1005,0,0,3690,240,0,0]
b1=plt.bar(x, c1, width=0.5, color='steelblue')
b2=plt.bar(x, c2, width=0.5, bottom=c1, color='mediumseagreen')
b3=plt.bar(x, c3, width=0.5, bottom=np.array(c1)+np.array(c2), color='goldenrod')
b4=plt.bar(x, c4, width=0.5, bottom=np.array(c1)+np.array(c2)+np.array(c3), color='indianred')
plt.xticks(x, ['A','B','C','D','E','F','G'])
plt.xlabel('ID')
plt.ylabel('counts')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend((b1[0], b2[0], b3[0], b4[0]), ('EP=1', 'EP=2', 'EP=3', 'EP=4'), loc='upper right')
plt.savefig('bar.png',dpi=300)