import os
import numpy as np
import dynalysis.basics as bcs
import dynalysis.classes as clss

#=====basics=====#
def vect_len(coor1, coor2):
    '''returns: length of vector'''
    if len(coor1)!=len(coor2): raise bcs.AlgorithmError('func vect_len')
    vect = [coor1[i]-coor2[i] for i in range(len(coor1))]
    return np.linalg.norm(np.array(vect))

def much_larger_than(a, b, criteria):
    '''returns: T/F(a >> b*criteria)'''
    if a > b*criteria and b!=0: return True
    else: return False
	
def deter_order(coors):
    if coors[0][2]>=coors[1][2]: return coors
    else: return list(reversed(coors))
	
def lower_than(l, criteria):
	'''returns: T/F(all items in l is smaller than criteria)'''
	for item in l:
		if item>criteria: return False
	return True
	
#=====classes=====#
class attractor():
    '''Defines an attractor.'''
    def __init__(self, name, coor, basin):
        self.name = name
        self.coor = coor
        self.basin = basin
        self.in_basin = False
    def deter_in_basin(self, state):
        if vect_len(state, self.coor) <= self.basin: return True
        return False
    def update_status(self, new_status):
        self.in_basin = new_status
    def check_overlap(self, other):
        if vect_len(self.coor, other.coor) < (self.basin+other.basin): return True
        return False
    
class storage():
    def __init__(self, init_name=[], init_type=[]):
        self.massive = {}
        for t in range(len(init_type)):
            tpe = init_type[t]
            if tpe is 'dict': self.massive[init_name[t]]={}
            elif tpe is 'list': self.massive[init_name[t]]=[]
    def add_garage(self, new_name, new_type):
        if new_type is 'dict': self.massive[new_name]={}
        elif new_type is 'list': self.massive[new_name]=[]
    def store_list(self, garage, item):
        self.massive[garage].append(item)
    def store_dict(self):
        pass
    def retrieve_garage(self, garage):
        return self.massive[garage]

class parameter():
	'''Defines a set of variables [pset] with corr. names [pname].
	For example, [coor, radius]=[(1,3), 5]'''
	def __init__(self, pset, pname):
		if len(pset)!=len(pname): raise bcs.InputError('class parameter')
		self.pset = pset
		self.pname = pname
		self.parampair = self.get_ordered_parampair()
		self.name = self.get_ordered_name()
	def get_ordered_parampair(self):
		parampair=[]
		for i in range(len(self.pname)): parampair.append((self.pname[i],self.pset[i]))
		return parampair
	def get_ordered_name(self):
		org_str=''
		for pair in self.parampair: org_str=org_str+pair[0]+'='+str(pair[1])+'_'
		return org_str[:-1]
	def update_all(self):
		self.parampair = self.get_ordered_parampair()
		self.name = self.get_ordered_name()
	def add_pair(self, pair): #interactive
		if type(pair)==tuple: pair=[pair]
		for p in pair:
			self.pset.append(p[1])
			self.pname.append(p[0])
		self.update_all()
	def remove_by_name(self, name):
		sett=self.extract(name)
		self.pname.remove(name)
		self.pset.remove(sett)
		self.update_all()
	def assign_parampair(self, parampair): #interactive
		self.pset = [pair[1] for pair in parampair]
		self.pname = [pair[0] for pair in parampair]
		self.update_all()
	def assign_name(self, name, trans_method=None): #interactive
		strpairlist = name.split('_')
		self.pset=[]
		for sp in range(len(strpairlist)):
			subject=strpairlist[sp].split('=')[1]
			if trans_method==None: #if trans_method is not specified
				try:
					if float(subject)==int(subject): self.pset.append(int(subject))
					else: self.pset.append(float(subject))
				except:
					self.pset.append(subject)
			else:
				if trans_method[sp]=='f': trans_subject=float(subject)
				elif trans_method[sp]=='i': trans_subject=int(float(subject))
				else: trans_subject = subject
				self.pset.append(trans_subject)
		self.pname = [strpair.split('=')[0] for strpair in strpairlist]
		self.update_all()
	def equivalence_check_by_name(self, othername): #interactive
		otherpm = parameter([],[])
		otherpm.assign_name(othername)
		for ps in self.pset:
			if ps not in otherpm.pset: return False
		for pn in self.pname:
			if pn not in otherpm.pname:return False
		return True
	def extract(self, pname): #interactive
		for pair in self.parampair:
			if pname==pair[0]: return pair[1]
		raise bcs.InputError('class parameter')
    
class bunch():
    '''Defines a variable with all of its possible values.
    For example, radius:[1,2,3,4,5]'''
    def __init__(self, name, sets):
        self.name = name
        self.sets = sets
        self.namebunch = self.get_ordered_namebunch()
        self.pairbunch = self.get_ordered_pairbunch()
    def get_ordered_namebunch(self):
        namebunch=[]
        for s in self.sets: namebunch.append(parameter([s],[self.name]).name) 
        return namebunch
    def get_ordered_pairbunch(self):
        pairbunch=[]
        for s in self.sets: pairbunch.append((self.name, s))
        return pairbunch            

#=====Others=====#
def add_continuity(runnum, outfile='info.txt', motherpath=os.getcwd()):
	#dir
	runpath=os.path.join(os.getcwd(),'run'+runnum)
	os.chdir(runpath)
	b_res=clss.branch('results_'+str(runnum), motherpath)
	alldirs=[dr for dr in os.listdir(runpath) if os.path.isdir(os.path.join(runpath,dr))]
	ofile=os.path.join(b_res.pathlink,outfile)
	Ent=clss.entry(' 0', [' 1', ' 2', ' 3', ' 4'])
	data=Ent.readdata_and_fix(ofile)
	bcs.output_clf(ofile)
	#analysis
	count=0
	for dr in alldirs:
		count+=1
		print(count)
		b=clss.branch(dr,runpath)
		rp=10 if len(os.listdir(b.pathlink))>9 else 1 #due to flawed dataset
		state_train_all, state_train_truck, state_neurons = get_state_train(b, repeat=rp)
		pwidth=continuity(state_neurons)
		pall='_'.join(bcs.to_string(pwidth))
		data[dr]=list(data[dr][:-1])+[data[dr][-1].split('\n')[0], pall]
	for key in data.keys():
		bcs.output_line(ofile,' '.join([key]+list(data[key])))
	os.chdir(motherpath)
	
def add_correlation(runnum, outfile='corr.txt', motherpath=os.getcwd()):
	#dir
	runpath=os.path.join(os.getcwd(),'run'+runnum)
	os.chdir(runpath)
	alldirs=[dr for dr in os.listdir(runpath) if os.path.isdir(os.path.join(runpath,dr))]
	#result files and outputs
	b_res=clss.branch('results_'+str(runnum), motherpath)
	bcs.output_clf(os.path.join(b_res.pathlink,outfile))
	#analysis
	for dr in alldirs:
		print(dr)
		#specifications
		b=clss.branch(dr,runpath)
		rp=10 if len(os.listdir(b.pathlink))>9 else 1 #due to flawed dataset
		#get EPs
		state_train_all, state_train_truck, state_neurons = get_state_train(b, repeat=rp)
		all_corrs=correlation(state_neurons)
		bcs.output_line(os.path.join(b_res.pathlink,outfile),\
						' '.join([dr]+all_corrs))
	os.chdir(motherpath)
	return 0

def get_parampairs(paramsets, bunches, iternum):
    if iternum <=0: return paramsets
    else:
        new_paramsets, bunch = [], bunches[len(bunches)-iternum].pairbunch
        for pair in bunch:
            for param in paramsets: new_paramsets.append(param+[pair])
        return get_parampairs(new_paramsets, bunches, iternum-1)

def analyze_bunch(bunches, motherpath=os.getcwd()):
    parampairs = get_parampairs([], bunches, len(bunches))
    for parampair in parampairs:
        pm = parameter([],[])
        pm.update_by_parampair(parampair)
        fname = os.path.join(motherpath, 'Frate_'+pm.name+'.txt')
        FF = get_Fano_Factor(fname)
		
def fix():
	import shutil
	b_run=clss.branch('run8',os.getcwd())
	b_fix=clss.branch('run8-3',os.getcwd())
	b_fix.mkdir()
	alldirs=[dr for dr in os.listdir(b_run.pathlink) if os.path.isdir(os.path.join(b_run.pathlink,dr))]
	for dr in alldirs:
		target=os.path.join(b_run.pathlink,dr)
		if 'Frate.txt' not in os.listdir(target):
			shutil.move(target, b_fix.pathlink)
	return 0


def fix2():
	import shutil
	b_run=clss.branch('run6',os.getcwd())
	b_fix=clss.branch('run6-3',os.getcwd())
	b_fix.mkdir()
	alldirs=[dr for dr in os.listdir(b_run.pathlink) if os.path.isdir(os.path.join(b_run.pathlink,dr))]
	for dr in alldirs:
		target=os.path.join(b_run.pathlink,dr,'Frate.txt')
		try:
			data=bcs.readcolumn(target)[1]
		except:
			print(dr)
			shutil.move(target, b_fix.pathlink)
	return 0
	
def detect_loop(pairlist, loopnum):
	tree=[list(item) for item in pairlist]
	for pair in pairlist:
		first=True
		for branch in tree:
			if pair[0]==branch[-1]:
				tree.append(branch+[pair[1]])
				first=False
	for branch in tree:
		#if the first and last element of the branch is the same, and
		#if the length of the branch is loopnum+1 (to prevent 1-2-1-3-1 counts for loopnum=3)
		#if the number of unique neurons in branch is loopnum  (to prevent 1-2-1-3-1 counts for loopnum=4)
		if branch[0]==branch[-1] and len(branch)==loopnum+1 and len(set(branch))==loopnum: return True
	return False
					