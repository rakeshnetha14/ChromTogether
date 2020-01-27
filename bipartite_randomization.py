import os
import sys
import numpy as np
import random
import copy
import csv
import networkx as nx

def bipartite_randomizations(file1,file2,nb):
	"""
	This function calculates for each TF pairs the number of spatial chromatin interactions on which they are co-occurring.
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
	lB=[];
	for i in range(max(realdegreeseq)+1):
		lB.append([])
	for i in range(len(realBedges)):
		if ddg.has_key(realBedges[i][0]):
			lB[ddg[realBedges[i][0]]].append(i)
	d1={}
	j1=0
	for i in range(len(realBedges)):
		ll=realBedges[i][0:2];
		if d1.has_key(tuple(ll)):
			d1[tuple(ll)]+=1
			j1+=1
		else:
			d1[tuple(ll)]=1
	rxy1=np.zeros([ntfs,ntfs],dtype=int)
	for i in range(len(realGedges)):
		for z in d2[int(realGedges[i][0])]:
			for z1 in d2[int(realGedges[i][1])]:
				if(d1.has_key((int(realGedges[i][0]),len(s1)+z)) and d1.has_key((int(realGedges[i][1]),len(s1)+z1))):
					rxy1[z][z1]+=1;
	rxy2=rxy1+np.transpose(rxy1)
	rxy=np.zeros([ntfs,ntfs],dtype=int)
	for i in range(len(s1)):
		for z in d2[int(i)]:
			for z1 in d2[int(i)]:
				rxy[z][z1]+=1;	
	#return rxy2,rxy
	"""
	This function creates random biparite graphs from real bipartite network between chromatin fragments and TFs. For each random network co-occurrence of TF pairs in spatial interactions is calculated and is compared with real network co-occurrence for each TF pair.
	"""
	countn=np.zeros([ntfs,ntfs],dtype=int);
	countn1=np.zeros([ntfs,ntfs],dtype=int);
	os.system('mkdir randomizations')
	name1='randomizations/rand_'
	cfrl1=[-1,0,1]
	cfrl3=[0,1]
	cfrl2=[-2,-1,0,1,2]
	d3={}
	d3=copy.deepcopy(d1)
	d4=copy.deepcopy(d2)
	randomBedges=copy.deepcopy(realBedges)
	for ii in range(nb):
		print ii
		xy1=np.zeros([ntfs,ntfs],dtype=int);
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
				ren1=np.random.choice(len(lB[int(wdf)]),1)
				rer2=randomBedges[lB[int(wdf)][ren1[0]]][0]
				ret2=randomBedges[lB[int(wdf)][ren1[0]]][1]
				retp2=randomBedges[lB[int(wdf)][ren1[0]]][2]
			else:
				wdf=random.choice(np.unique(ddg.values())[10:])
				ren1=np.random.choice(len(lB[int(wdf)]),1)
				rer2=randomBedges[lB[int(wdf)][ren1[0]]][0]
				ret2=randomBedges[lB[int(wdf)][ren1[0]]][1]
				retp2=randomBedges[lB[int(wdf)][ren1[0]]][2]
			if(d3.has_key(tuple([rer1,ret2])) or d3.has_key(tuple([rer2,ret1]))):
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
		for i in range(len(realGedges)):
			for z in d4[int(realGedges[i][0])]:
				for z1 in d4[int(realGedges[i][1])]:
					#if(d3.has_key((int(realGedges[i][0]),len(s1)+z)) and d3.has_key((int(realGedges[i][1]),len(s1)+z1))):
					xy1[z][z1]+=1;
		xy2=xy1+np.transpose(xy1)	
		with open(name1+str(ii)+'.csv', "wb") as f:
    			writer = csv.writer(f)
    			writer.writerows(xy2)
		for z in range(ntfs):
			for z1 in range(ntfs):
				if(xy2[z][z1]>=rxy2[z][z1]):
					countn[z][z1]+=1;
				if(xy2[z][z1]<=rxy2[z][z1]):
					countn1[z][z1]+=1;
	with open('randomizations/ranp.csv', "wb") as f:
    		writer = csv.writer(f)
    		writer.writerows(countn)
	with open('randomizations/rann.csv', "wb") as f:
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

def qvalues(file1,file2,nb,ntfs):
	"""
	This function calculates the empirical q-values for both attracting and repelling TF pairs.
	"""
	pv=[];
	with open(file1, 'rb') as f2:
		reader = csv.reader(f2)
		for row in reader:
			for z1 in range(ntfs):
				pv.append((int(row[z1])+1)/float(nb+1))
	pv1=[]
	with open(file2, 'rb') as f3:
		reader = csv.reader(f3)
		for row in reader:
			for z1 in range(ntfs):
				pv1.append((int(row[z1])+1)/float(nb+1))
	fdr,t=fdrestimate(pv,ntfs)
	q=qvalueCalculate(pv,fdr,t)
	fdr,t=fdrestimate(pv1,ntfs)
	qdash=qvalueCalculate(pv1,fdr,t)
	f4=open('qvalue_matrix.dat','w')
	count=0
	for i in range(ntfs):
		for j in range(ntfs):
			if pv[count]>=0.8:
				value=1-qdash[count]
			else:
				if q[count]>=0.95:
					value = 0.935
				else:
					value=q[count] 
			f4.write(str(value)+'\t')
			count+=1
		f4.write('\n')



chromatinfile=sys.argv[1]
nr1=sys.argv[2]
[chromatinfile_name,chromatinfile_extension]=chromatinfile.split('.')
chromatinfile19=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_asnodes.'+chromatinfile_extension
chromatinfile20=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge_oc.'+chromatinfile_extension

ntfs1=bipartite_randomizations(chromatinfile19,chromatinfile20,nr1)
qvalues('randomizations/ranp.csv','randomizations/rann.csv',nr1,ntfs1)



