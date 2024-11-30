import os
import numpy as np
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
import dynalysis.basics as bcs
import dynalysis.classes as clss
from dynalysis.scalebar import add_scalebar

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def palette(name,num):
	return sns.color_palette(name,num)

c1='#1f77b4' #blue
c4='white'
c2='white'
c3 = 'firebrick'
scale = (1,50)

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
data = bcs.readcolumn('Frate.txt')
time = bcs.to_float(data[0])
for neu in range(1,75+1): plt.plot(time,bcs.to_float(data[neu]),color=palette("BrBG",75)[neu-1],linewidth=1.5)
ax.axis('off')
scalebar_h = add_scalebar(ax, scale[0], '', 'upper right', matchx=False, matchy=False)
scalebar_v = add_scalebar(ax, scale[1], '', 'lower left', matchx=False, matchy=False, horizontal=False)
ax.add_artist(scalebar_v)
ax.add_artist(scalebar_h)
plt.savefig('exc.png',dpi=300); plt.clf() #colorFader(c4,c1,neu/75)

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
for neu in range(76,101): plt.plot(time,bcs.to_float(data[neu]),color=palette("BrBG",25)[neu-76],linewidth=1.5)
ax.axis('off')
scalebar_h = add_scalebar(ax, scale[0], '', 'upper right', matchx=False, matchy=False)
scalebar_v = add_scalebar(ax, scale[1], '', 'lower left', matchx=False, matchy=False, horizontal=False)
ax.add_artist(scalebar_v)
ax.add_artist(scalebar_h)
plt.savefig('inh.png',dpi=300); plt.clf() #colorFader(c2,c3,(neu-76)/25)

	