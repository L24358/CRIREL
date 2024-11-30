import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import dynalysis.basics as bcs
import dynalysis.classes as clss
from matplotlib import colors, gridspec

#=====subfunctions=====#
def read_event_dic(ed, trial_t=None):
    neu_dic, max_time={},0
    for key in sorted(ed.keys()):
        start_time, neuron = key
        #pulse or ramp
        if ed[key][0]=='pulse':
            duration, mag=ed[key][1], ed[key][2]
            to_append=[(start_time,mag),(start_time+duration,0)]
            end_time=start_time+duration
        elif ed[key][0]=='ramp':
            to_append=[(start_time,ed[key][1])]
            end_time=start_time
        elif ed[key][0]=='pulse_b':
            duration, mag, bias=ed[key][1], ed[key][2], ed[key][4]
            to_append=[(start_time,mag),(start_time+duration,bias)]
            end_time=start_time+duration
        elif ed[key][0]=='ramp_b':
            to_append=[(start_time,ed[key][1]+ed[key][3])]
            end_time=start_time
        else: raise bcs.InputError('func read_event_dic')
        #maxtime
        if end_time > max_time: max_time=end_time
        #append neuron to neu_dic key list
        if neuron not in neu_dic.keys(): neu_dic[neuron]={'x':[0],'y':[0]}
        #append the magnitudes
        for t in to_append:
            neu_dic[neuron]['x'].append(t[0])
            neu_dic[neuron]['x'].append(t[0])
            neu_dic[neuron]['y'].append(neu_dic[neuron]['y'][-1])
            neu_dic[neuron]['y'].append(t[1])
    if trial_t!=None: max_time=trial_t
    for neuron in neu_dic.keys():
        neu_dic[neuron]['x'].append(max_time)
        neu_dic[neuron]['y'].append(neu_dic[neuron]['y'][-1])
    return neu_dic

def read_data(col, fname='Frate.txt'):
    all_cols=bcs.readcolumn(fname)
    time, to_plot=bcs.to_float(all_cols[0]), bcs.to_float(all_cols[col])
    return time, to_plot

#=====main=====#
def edplot(ed, trial_t=None, **kwargs):
    #kwargs
    kw={'label_fontsize':16, 'tick_fontsize':12, 'outputfigname':'protocol.png',\
        'xlabel':'time (s)', 'title':None, 'figsize':(12,8), 'keys':[]}
    kw.update(kwargs)
    label_fontsize, tick_fontsize=kw['label_fontsize'], kw['tick_fontsize']
    outputfigname, keys=kw['outputfigname'], kw['keys']
    xlabel, title, figsize=kw['xlabel'], kw['title'], kw['figsize']
    #main
    neu_dic=read_event_dic(ed)
    if keys==[]:
        all_cols=[neu_dic[neuron] for neuron in clss.sort_nodes(neu_dic.keys())]
        col_labels=[neuron for neuron in clss.sort_nodes(neu_dic.keys())]
    else:
        all_cols=[neu_dic[neuron] for neuron in clss.sort_nodes(keys)]
        col_labels=[neuron for neuron in clss.sort_nodes(keys)]
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i in range(len(col_labels)):
        ax = fig.add_subplot(len(col_labels), 1, i+1)
        ax.set_ylabel(col_labels[i], fontsize=label_fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.plot(all_cols[i]['x'], all_cols[i]['y'], 'k')
    plt.xlabel(xlabel, fontsize=label_fontsize)
    if title!=None: fig.suptitle(title, fontsize=20)
    plt.tight_layout()
    plt.savefig(outputfigname, dpi=150)
    plt.close('all')
    return all_cols

def subplot(cols=[1,2,3,4], fname='Frate.txt', **kwargs):
    #kwargs
    kw={'label_fontsize':16, 'tick_fontsize':12, 'outputfigname':'Frate.png',\
        'col_labels':['exc1','exc2','inh1','inh2'], 'xlabel':'time (s)',\
        'title':None, 'figsize':(12,8)}
    kw.update(kwargs)
    label_fontsize, tick_fontsize=kw['label_fontsize'], kw['tick_fontsize']
    outputfigname, col_labels=kw['outputfigname'], kw['col_labels']
    xlabel, title, figsize=kw['xlabel'], kw['title'], kw['figsize']
    #main
    all_cols=bcs.readcolumn(fname)
    time, to_plot=bcs.to_float(all_cols[0]), [bcs.to_float(all_cols[c]) for c in cols]
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i in range(len(cols)):
        ax = fig.add_subplot(len(cols), 1, i+1)
        ax.set_xticks(list(np.arange(0,time[-1],0.5)))
        ax.set_ylabel(col_labels[i], fontsize=label_fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.plot(time, to_plot[i], 'k')
    plt.xlabel(xlabel, fontsize=label_fontsize)
    if title!=None: fig.suptitle(title, fontsize=20)
    plt.tight_layout()
    plt.savefig(outputfigname, dpi=150)
    plt.close('all')
    return 0

def subplot2(xlist, ylist, **kwargs):
    #kwargs
    kw={'label_fontsize':16, 'tick_fontsize':12, 'outputfigname':'any.png',\
        'title':None, 'figsize':(12,8), 'xcol_labels':None, 'ycol_labels':None}
    kw.update(kwargs)
    label_fontsize, tick_fontsize=kw['label_fontsize'], kw['tick_fontsize']
    outputfigname=kw['outputfigname']
    title, figsize=kw['title'], kw['figsize']
    xcol_labels=['time (s)']*len(xlist) if kw['xcol_labels']==None else kw['xcol_labels']
    ycol_labels =['fr']*len(ylist) if kw['ycol_labels']==None else kw['ycol_labels']
    #main
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i in range(len(xlist)):
        ax = fig.add_subplot(len(xlist), 1, i+1)
        ax.set_xlabel(xcol_labels[i], fontsize=label_fontsize)
        ax.set_ylabel(ycol_labels[i], fontsize=label_fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.plot(xlist[i], ylist[i], 'k')
    if title!=None: fig.suptitle(title, fontsize=20)
    plt.tight_layout()
    plt.savefig(outputfigname, dpi=150)
    plt.close('all')
    return 0

def subplot3(xlist, ylist, indexlist, **kwargs):
    #kwargs
    kw={'label_fontsize':16, 'tick_fontsize':12, 'outputfigname':'any.png',\
        'title':None, 'figsize':(12,8), 'xcol_labels':None, 'ycol_labels':None,\
        'titlelist':None, 'collist':None, 'title_fontsize':20, 'wspace':None,\
        'hspace':None,'offaxis':[]}
    kw.update(kwargs)
    label_fontsize, tick_fontsize=kw['label_fontsize'], kw['tick_fontsize']
    outputfigname,titlelist=kw['outputfigname'],kw['titlelist']
    title, figsize,collist=kw['title'], kw['figsize'],kw['collist']
    xcol_labels=['time (s)']*len(xlist) if kw['xcol_labels']==None else kw['xcol_labels']
    ycol_labels =['fr']*len(ylist) if kw['ycol_labels']==None else kw['ycol_labels']
    title_fontsize,wspace,hspace=kw['title_fontsize'],kw['wspace'],kw['hspace']
    offaxis=kw['offaxis']
    #main
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i in range(len(xlist)):
        ax = fig.add_subplot(int(indexlist[i][1]), int(indexlist[i][0]), int(indexlist[i][2:]))
        ax.set_xlabel(xcol_labels[i], fontsize=label_fontsize)
        ax.set_ylabel(ycol_labels[i], fontsize=label_fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i+1 in offaxis: ax.set_axis_off()
        if collist!=None and collist[i]!=None: color=collist[i]
        else: color='k'
        if titlelist!=None and titlelist[i]!=None:
            ax.set_title(titlelist[i],fontsize=title_fontsize)
        ax.plot(xlist[i], ylist[i], color)
    if title!=None: fig.suptitle(title, fontsize=title_fontsize)
    if wspace==None and hspace==None: plt.tight_layout()
    if wspace!=None: plt.subplots_adjust(wspace=wspace)
    if hspace!=None: plt.subplots_adjust(hspace=hspace)
    plt.savefig(outputfigname, dpi=150)
    plt.close('all')
    return 0

def linesplot(xlist, ylist, **kwargs):
    #kwargs
    kw={'label_fontsize':16, 'tick_fontsize':12, 'outputfigname':'line.png',\
        'title':None, 'figsize':(12,8), 'xlabel':'time (s)', 'ylabel':'Fr (Hz)',\
        'collist':None,'ticksoff':False}
    kw.update(kwargs)
    label_fontsize, tick_fontsize=kw['label_fontsize'], kw['tick_fontsize']
    outputfigname,ticksoff=kw['outputfigname'],kw['ticksoff']
    title, figsize=kw['title'], kw['figsize']
    xlabel,ylabel,collist=kw['xlabel'],kw['ylabel'],kw['collist']
    #main
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    ax = fig.add_subplot(111)
    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if ticksoff:
        plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
        plt.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    else:
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if collist==None:
        for i in range(len(xlist)): ax.plot(xlist[i], ylist[i])
    else:
        for i in range(len(xlist)): ax.plot(xlist[i], ylist[i], collist[i])
    if title!=None: fig.suptitle(title, fontsize=20)
    plt.tight_layout()
    plt.savefig(outputfigname, dpi=150)
    plt.close('all')
    return 0

def sublinesplot(xlist, ylist, indexlist, **kwargs):
    #kwargs
    kw={'label_fontsize':16, 'tick_fontsize':12, 'outputfigname':'any.png',\
        'title':None, 'figsize':(12,8), 'xcol_labels':None, 'ycol_labels':None,\
        'titlelist':None, 'collist':None, 'title_fontsize':20, 'wspace':None,\
        'hspace':None,'offaxis':[], 'dpi':150, 'points':[],'emptypoints':[]}
    kw.update(kwargs)
    label_fontsize, tick_fontsize=kw['label_fontsize'], kw['tick_fontsize']
    outputfigname,titlelist=kw['outputfigname'],kw['titlelist']
    title, figsize,collist=kw['title'], kw['figsize'],kw['collist']
    xcol_labels=['time (s)']*len(xlist) if kw['xcol_labels']==None else kw['xcol_labels']
    ycol_labels =['fr']*len(ylist) if kw['ycol_labels']==None else kw['ycol_labels']
    title_fontsize,wspace,hspace=kw['title_fontsize'],kw['wspace'],kw['hspace']
    offaxis,dpi,points,emptypoints=kw['offaxis'],kw['dpi'],kw['points'],kw['emptypoints']
    #main
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i in range(len(xlist)):
        ax = fig.add_subplot(int(indexlist[i][1]), int(indexlist[i][0]), int(indexlist[i][2:]))
        ax.set_xlabel(xcol_labels[i], fontsize=label_fontsize)
        ax.set_ylabel(ycol_labels[i], fontsize=label_fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i+1 in offaxis: ax.set_axis_off()
        if titlelist!=None and titlelist[i]!=None:
            ax.set_title(titlelist[i],fontsize=title_fontsize)
        for j in range(len(xlist[i])):
            if (i,j) in points: lstyle, marker,face='None', 'o','original'
            elif (i,j) in emptypoints: lstyle, marker,face='None', 'o','none'
            else: lstyle, marker,face='-', ',','original'
            if collist==None:
                mfc='k' if face=='original' else face
                ax.plot(xlist[i][j], ylist[i][j], linestyle=lstyle, marker=marker, mfc=mfc)
            else:
                mfc=collist[i][j] if face=='original' else face
                ax.plot(xlist[i][j], ylist[i][j], color=collist[i][j], linestyle=lstyle, marker=marker, mfc=mfc)
    if title!=None: fig.suptitle(title, fontsize=title_fontsize)
    if wspace==None and hspace==None: plt.tight_layout()
    if wspace!=None: plt.subplots_adjust(wspace=wspace)
    if hspace!=None: plt.subplots_adjust(hspace=hspace)
    plt.savefig(outputfigname, dpi=dpi)
    plt.close('all')
    return 0

def stateplot(statedic, **kwargs):
    '''Plots the evolution of the state of a variable through time'''
    #kwargs
    kw={'state2color':{}, 'objectorder':None, 'halfspan':0.1, 'figsize':None,\
        'outputfigname':'stateplot.png', 'colorpanel':None, 'stateorder':None}
    kw.update(kwargs)
    state2color, objectorder=kw['state2color'], kw['objectorder']
    halfspan, figsize, outputfigname=kw['halfspan'], kw['figsize'], kw['outputfigname']
    colorpanel, stateorder=kw['colorpanel'], kw['stateorder']
    #declare
    allstates, maxlength=[], 0
    objectlist=objectorder if type(objectorder)==list else sorted(statedic.keys())
    #get all states and maxlength
    for obj in objectlist:
        for state in statedic[obj].keys():
            if state not in allstates: allstates.append(state)
            for instance in statedic[obj][state]:
                maxlength=instance[1] if instance[1] > maxlength else maxlength
    #set state2color, panel and stateorder
    if state2color=={}:
        panel=sns.color_palette('Set2',len(allstates))
        for s in range(len(allstates)): state2color[allstates[s]]=panel[s]
    else: panel=colorpanel
    stateorder=sorted(allstates) if type(stateorder)!=list else stateorder
    #cmap
    fig = plt.figure() if figsize==None else plt.figure(figsize)
    gs = gridspec.GridSpec(1, 2, width_ratios=[20, 1]) 
    ax2 = plt.subplot(gs[1])
    col_map = colors.ListedColormap(panel)
    cbar = mpl.colorbar.ColorbarBase(ax2, cmap=col_map, orientation='vertical',\
                                     boundaries=[0,1,2,3,4],\
                                     ticks=np.arange(0.5,0.5+len(allstates),1))
    cbar.set_ticklabels(stateorder)
    #plot lines
    ax1 = plt.subplot(gs[0])
    for obj in objectlist:
        y_at=objectlist.index(obj)+1
        for state in statedic[obj].keys():
            for instance in statedic[obj][state]:
                mn,mx=float(instance[0])/maxlength,float(instance[1])/maxlength
                plt.axhspan(y_at-halfspan, y_at+halfspan,xmin=mn,xmax=mx,\
                            facecolor=state2color[state])
    #settings
    #ax1.spines['top'].set_visible(False)
    #ax1.spines['right'].set_visible(False)
    plt.yticks([i for i in range(1,len(objectlist)+1)],objectlist)
    ax1.set_xlim([0,maxlength])
    ax1.set_ylim([0,len(objectlist)+1])
    plt.savefig(outputfigname, dpi=150)
    
