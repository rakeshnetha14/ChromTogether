import os
import sys
import numpy as np
import random
import copy
import csv
import networkx as nx
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

os.system('mkdir randomizations')

def bipartite_randomizations(file1,file2,nb,peaks,dimension):
	"""
	This function calculates the number of spatial chromatin interactions for each of the TF pair on which the TF pair are co-occurring.
	"""
	f=open(file1,'r');
	s=[line.split('\t') for line in f];
	f.close()
	#####
	realGedges=[(int(s[i][5]),int(s[i][6])) for i in range(len(s))];
	realG = nx.Graph()
	realG.add_edges_from(realGedges)
	#print realG.number_of_nodes()
	realdegree=list(realG.degree)
	realdegreeseq=[realdegree[i][1] for i in range(len(realdegree))]
	#realdegreeseq=list(realG.degree().values())
	#realdegreeseq=[realG.degree()[i] for i in range(realG.number_of_nodes())]
	ddg={}
	for i in range(len(realdegreeseq)):
		ddg[i]=realdegreeseq[i]
	ddg_uni=np.sort(np.unique(list(ddg.values())))
	tr=0;
	while ddg_uni[tr]<=10:
		tr+=1;
	l=[];
	for i in range(max(realdegreeseq)+1):
		l.append([])
	for i in range(len(realdegreeseq)):
		l[ddg[i]].append(i)
	#####
	f1=open(file2,'r');
	s1=[line.split('\t') for line in f1];
	f1.close()
	ntfs=len(s1[0])-3;
	if dimension=='3D':
		d2={}
		realBedges=[];
		BB={};
		k=0
		for i in range(len(s1)):
			TFlist=[];
			for j in range(ntfs):
				if(int(float(s1[i][3+j]))==1):
					TFlist.append(j)
					realBedges.append([i,len(s1)+j,k]);
					BB[k]=i;
					k+=1;
			d2[i]=TFlist
	elif dimension=='1D':
		d2={}
		realBedges=[];
		BB={};
		k=0
		for i in range(len(s1)):
			TFlist=[];
			for j in range(ntfs):
				if(int(float(s1[i][3+j]))==2):
					TFlist.append(2*j)
					TFlist.append((2*j)+1)
					realBedges.append([i,len(s1)+(2*j),k]);
					realBedges.append([i,len(s1)+(2*j)+1,k]);
					BB[k]=i;
					BB[k+1]=i;
					k+=2;
				elif(int(float(s1[i][3+j]))==1):
					TFlist.append(2*j)
					realBedges.append([i,len(s1)+(2*j),k]);
					BB[k]=i;
					k+=1;
			d2[i]=TFlist
	lB=[];
	for i in range(max(realdegreeseq)+1):
		lB.append([])
	for i in range(len(realBedges)):
		if realBedges[i][0] in ddg:
			lB[ddg[realBedges[i][0]]].append(i)
	d1={}
	j1=0
	for i in range(len(realBedges)):
		ll=realBedges[i][0:2];
		if tuple(ll) in d1:
			d1[tuple(ll)]+=1
			j1+=1
		else:
			d1[tuple(ll)]=1
	if dimension=='3D':
		rxy1=np.zeros([ntfs,ntfs],dtype=int)
		for i in range(len(realGedges)):
			for z in d2[int(realGedges[i][0])]:
				for z1 in d2[int(realGedges[i][1])]:
					if ((int(realGedges[i][0]),len(s1)+z) in d1 and (int(realGedges[i][1]),len(s1)+z1) in d1):
						rxy1[z][z1]+=1;
		rxy2=rxy1+np.transpose(rxy1)
	elif dimension=='1D':
		rxy=np.zeros([ntfs,ntfs],dtype=int)
		for i in range(len(s1)):
			for z in d2[int(i)]:
				for z1 in d2[int(i)]:
					if z!=z1:
						rxy[z//2][z1//2]+=1;	
	#return rxy2,rxy
	"""
	The following part of function creates random biparite graphs from the real bipartite network between chromatin fragments and TFs. 
	For each random network, co-occurring spatial interactions are measured for each TF pair and is compared with the corresponding number in the
	real network.
	"""
	countn=np.zeros([ntfs,ntfs],dtype=int);
	countn1=np.zeros([ntfs,ntfs],dtype=int);
	#countn2=np.zeros([ntfs,ntfs],dtype=int);
	#countn3=np.zeros([ntfs,ntfs],dtype=int);
	name1='randomizations/rand_'+peaks+'_'+dimension+'_'
	#name2='randomizations/rand_'+peaks+'_'+dimension+'_'
	cfrl1=[-1,0,1]
	cfrl3=[0,1]
	cfrl2=[-2,-1,0,1,2]
	d3={}
	d3=copy.deepcopy(d1)
	d4=copy.deepcopy(d2)
	randomBedges=copy.deepcopy(realBedges)
	for ii in range(nb):
		print(ii)
		xy1=np.zeros([ntfs,ntfs],dtype=int);
		xy=np.zeros([ntfs,ntfs],dtype=int);
		jj=0
		for i in range(int(len(realBedges)*0.3)):
			ren=np.random.choice(len(realBedges),1)
			rer1=randomBedges[ren[0]][0]
			ret1=randomBedges[ren[0]][1]
			retp1=randomBedges[ren[0]][2]
			wd=ddg[BB[retp1]]
			if wd<=5:
				wdf=wd
				ren1=np.random.choice(len(lB[int(wdf)]),1)
				rer2=randomBedges[lB[int(wdf)][ren1[0]]][0]
				ret2=randomBedges[lB[int(wdf)][ren1[0]]][1]
				retp2=randomBedges[lB[int(wdf)][ren1[0]]][2]
			elif wd<=10:
				cfr=np.random.choice(5,1)
				wdf=wd+cfrl2[int(cfr)]
				while len(lB[int(wdf)])==0:
					cfr=np.random.choice(5,1)
					wdf=wd+cfrl2[int(cfr)]
				ren1=np.random.choice(len(lB[int(wdf)]),1)
				rer2=randomBedges[lB[int(wdf)][ren1[0]]][0]
				ret2=randomBedges[lB[int(wdf)][ren1[0]]][1]
				retp2=randomBedges[lB[int(wdf)][ren1[0]]][2]
			else:
				wdf=random.choice(ddg_uni[tr:])
				while len(lB[int(wdf)])==0:
					wdf=random.choice(ddg_uni[tr:])
				ren1=np.random.choice(len(lB[int(wdf)]),1)
				rer2=randomBedges[lB[int(wdf)][ren1[0]]][0]
				ret2=randomBedges[lB[int(wdf)][ren1[0]]][1]
				retp2=randomBedges[lB[int(wdf)][ren1[0]]][2]
			if(tuple([rer1,ret2]) in d3 or tuple([rer2,ret1]) in d3):
				pass
				#jj+=1
			else:
				randomBedges[ren[0]][1]=ret2
				randomBedges[ren[0]][2]=retp2
				randomBedges[lB[int(wdf)][ren1[0]]][1]=ret1
				randomBedges[lB[int(wdf)][ren1[0]]][2]=retp1
				del d3[tuple([rer1,ret1])]
				del d3[tuple([rer2,ret2])]
				d3[tuple([rer1,ret2])]=1
				d3[tuple([rer2,ret1])]=1
				d4[rer1].remove(ret1-len(s1))
				d4[rer1].append(ret2-len(s1))
				d4[rer2].remove(ret2-len(s1))
				d4[rer2].append(ret1-len(s1))
		if dimension=='3D':
			for i in range(len(realGedges)):
				for z in d4[int(realGedges[i][0])]:
					for z1 in d4[int(realGedges[i][1])]:
						xy1[z][z1]+=1;
			xy2=xy1+np.transpose(xy1)
			with open(name1+str(ii)+'.csv', 'w') as f:
    				writer = csv.writer(f)
    				writer.writerows(xy2)
			for z in range(ntfs):
				for z1 in range(ntfs):
					if(xy2[z][z1]>=rxy2[z][z1]):
						countn[z][z1]+=1;
					if(xy2[z][z1]<=rxy2[z][z1]):
						countn1[z][z1]+=1;
		elif dimension=='1D':
			for i in range(len(s1)):
				for z in d4[int(i)]:
					for z1 in d4[int(i)]:
						if z!=z1:
							xy[z//2][z1//2]+=1;	
			with open(name1+str(ii)+'.csv', 'w') as f:
    				writer = csv.writer(f)
    				writer.writerows(xy)
			for z in range(ntfs):
				for z1 in range(ntfs):
					if(xy[z][z1]>=rxy[z][z1]):
						countn[z][z1]+=1;
					if(xy[z][z1]<=rxy[z][z1]):
						countn1[z][z1]+=1;
	with open('randomizations/ranp_'+peaks+'_'+dimension+'.csv', 'w') as f:
    		writer = csv.writer(f)
    		writer.writerows(countn)
	with open('randomizations/rann_'+peaks+'_'+dimension+'.csv', 'w') as f:
    		writer = csv.writer(f)
    		writer.writerows(countn1)
	return ntfs


def fdrestimate(mylist,m):
	pi0=1
	t=np.linspace(0.001,1,1000)
	FDR=[]
	for j in range(len(t)):
		count=0
		for i in range(len(mylist)):
			if mylist[i]<=t[j]:
				count=count+1
		value=pi0*m*t[j]/float(count)
		FDR.append(value)
	
	return [FDR,t] 


def qvalueCalculate(mylist,fdr,t):
	qvalue=[]
	for i in range(len(mylist)):
		temp=[]
		for j in range(len(fdr)):
			if t[j]>=mylist[i]:
				temp.append(fdr[j])
		qvalue.append(min(temp))
	
	return qvalue 

def qvalues(file1,file2,nb,ntfs,peaks,dimension):
	"""
	This function calculates the empirical q-values for both attracting and repelling TF pairs.
	"""
	pv=[];
	with open(file1, 'r') as f2:
		reader = csv.reader(f2)
		for row in reader:
			for z1 in range(ntfs):
				pv.append((int(row[z1])+1)/float(nb+1))
	pv1=[]
	with open(file2, 'r') as f3:
		reader = csv.reader(f3)
		for row in reader:
			for z1 in range(ntfs):
				pv1.append((int(row[z1])+1)/float(nb+1))
	fdr,t=fdrestimate(pv,ntfs*ntfs)
	q=qvalueCalculate(pv,fdr,t)
	fdr,t=fdrestimate(pv1,ntfs*ntfs)
	qdash=qvalueCalculate(pv1,fdr,t)
	f4=open('qv_'+peaks+'_'+dimension+'.dat','w')
	count=0
	for i in range(ntfs):
		for j in range(ntfs):
			if pv[count]>=0.8:
				value=1-qdash[count]
			else:
				if q[count]>=0.95:
					value = 0.50
				else:
					value=q[count] 
			f4.write(str(value)+'\t')
			count+=1
		f4.write('\n')

def heatmap_generation(file1,file2,ntfs):
	"""
	This function generates heatmap showing both attracting and repelling TF pairs.
	"""
	cm=np.zeros([1000,4],dtype=float)
	for i in range(50):
		cm[i,1]=(1-((0.75/50)*i));
	for j in range(950,1000):
		cm[j,0]=0.265+((0.75/50)*(j-950));
	for k in range(1000):
		cm[k,3]=1
	cm1=colors.LinearSegmentedColormap.from_list('my_colormap', cm)
	f5=open(file1,'r');
	s5=[line.split('\t') for line in f5]
	f5.close()
	M=np.zeros([ntfs,ntfs],dtype=float);
	for i in range(ntfs):
		for j in range(ntfs):
			M[i][j]=float(s5[i][j])
	[file1_name,file1_extension]=file1.split('.')
	f6=open(file2,'r');
	s6=f6.readlines();
	f6.close();
	fig, ax = plt.subplots()
	ax1=ax.imshow(M,cmap=cm1, vmin=0, vmax=1)
	ax.grid(False)
	ax.set_xticks(range(ntfs))
	ax.set_yticks(range(ntfs))
	ax.tick_params(axis=u'both', which=u'both',length=3)
	ax.set_xticklabels(s6, rotation='vertical', fontsize=6, fontweight = 'bold')
	ax.xaxis.set_ticks_position('top')
	ax.set_yticklabels(s6, fontsize=6, fontweight = 'bold')
	fig.colorbar(ax1)
	plt.tight_layout(rect=[0,0,1,1])
	plt.savefig(file1_name+'.png',dpi=600,bbox_inches='tight')
	plt.savefig(file1_name+'.pdf',dpi=600,bbox_inches='tight')
	plt.close()



chromatinfile=sys.argv[1]
nr1=int(sys.argv[2])
com=int(sys.argv[3])
[chromatinfile_name,chromatinfile_extension]=chromatinfile.split('.')
chromatinfile19=chromatinfile_name+'_intra_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_asnodes.'+chromatinfile_extension
chromatinfile20=chromatinfile_name+'_intra_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge_oc.'+chromatinfile_extension
chromatinfile21=chromatinfile_name+'_intra_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge_ocd.'+chromatinfile_extension
chromatinfile22=chromatinfile_name+'_intra_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge_mc.'+chromatinfile_extension
chromatinfile23=chromatinfile_name+'_intra_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge_mcd.'+chromatinfile_extension

if com==2:
	ntfs1=bipartite_randomizations(chromatinfile19,chromatinfile20,nr1,'chip','3D')
	ntfs1=bipartite_randomizations(chromatinfile19,chromatinfile21,nr1,'chip','1D')
	qvalues('randomizations/ranp_chip_3D.csv','randomizations/rann_chip_3D.csv',nr1,ntfs1,'chip','3D')
	qvalues('randomizations/ranp_chip_1D.csv','randomizations/rann_chip_1D.csv',nr1,ntfs1,'chip','1D')
	heatmap_generation('qv_chip_3D.dat','tf_chipname.txt',ntfs1)
	heatmap_generation('qv_chip_1D.dat','tf_chipname.txt',ntfs1)
	ntfs2=bipartite_randomizations(chromatinfile19,chromatinfile22,nr1,'pwm','3D')
	ntfs2=bipartite_randomizations(chromatinfile19,chromatinfile23,nr1,'pwm','1D')
	qvalues('randomizations/ranp_pwm_3D.csv','randomizations/rann_pwm_3D.csv',nr1,ntfs2,'pwm','3D')
	qvalues('randomizations/ranp_pwm_1D.csv','randomizations/rann_pwm_1D.csv',nr1,ntfs2,'pwm','1D')
	heatmap_generation('qv_pwm_3D.dat','tf_pwmname.txt',ntfs2)
	heatmap_generation('qv_pwm_1D.dat','tf_pwmname.txt',ntfs2)
elif com==1:
	ntfs2=bipartite_randomizations(chromatinfile19,chromatinfile22,nr1,'pwm','3D')
	ntfs2=bipartite_randomizations(chromatinfile19,chromatinfile23,nr1,'pwm','1D')
	qvalues('randomizations/ranp_pwm_3D.csv','randomizations/rann_pwm_3D.csv',nr1,ntfs2,'pwm','3D')
	qvalues('randomizations/ranp_pwm_1D.csv','randomizations/rann_pwm_1D.csv',nr1,ntfs2,'pwm','1D')
	heatmap_generation('qv_pwm_3D.dat','tf_pwmname.txt',ntfs2)
	heatmap_generation('qv_pwm_1D.dat','tf_pwmname.txt',ntfs2)
elif com==0:
	ntfs1=bipartite_randomizations(chromatinfile19,chromatinfile20,nr1,'chip','3D')
	ntfs1=bipartite_randomizations(chromatinfile19,chromatinfile21,nr1,'chip','1D')
	qvalues('randomizations/ranp_chip_3D.csv','randomizations/rann_chip_3D.csv',nr1,ntfs1,'chip','3D')
	qvalues('randomizations/ranp_chip_1D.csv','randomizations/rann_chip_1D.csv',nr1,ntfs1,'chip','1D')
	heatmap_generation('qv_chip_3D.dat','tf_chipname.txt',ntfs1)
	heatmap_generation('qv_chip_1D.dat','tf_chipname.txt',ntfs1)



