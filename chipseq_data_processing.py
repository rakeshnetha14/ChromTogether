import os
import sys
import numpy as np

def bipartite_info(file1,chipseq_folder):
	"""
	This function creates a data file which contain the binding information for each TF in columns on each chromatin 
	fragment in the rows.
	"""
	chipf=os.listdir(chipseq_folder)
	chipf.sort()
	#####	
	f=open(file1,'r');
	s=[line.split('\t') for line in f];
	f.close();
	s_array=np.array(s)
	chID=np.unique(s_array[:,0])
	chID_d={}
	chID_d={}
	chID_d[s[0][0]]=0;
	#####
	ra=[0];
	chid=s[0][0];
	k=0;
	i1=0;
	for l in range(len(chID)-1):
		while(chid==s[k][0]):
			k+=1;
		ra.append(k);
		chID_d[s[k][0]]=i1+1;
		chid=s[k][0];
		i1+=1;
	ra.append(len(s));
	#####
	s=np.array(s);
	b1=s[:,1];
	b2=s[:,2];
	b1=list(map(int,b1));
	b2=list(map(int,b2));
	#####
	chipf_name=[]
	for z in range(len(chipf)):
		#print z;
		chipf_name.append(chipf[z].split('.')[0]);
		with open(chipseq_folder+'/'+chipf[z],'r') as f1:
			s1=[line.split('\t')[0:3] for line in f1];
		f1.close();
		s1_array=np.array(s1)
		chID1=np.unique(s1_array[:,0])
		chID1_d={}
		chID1_d[s1[0][0]]=0;
		ra1=[0];
		chid=s1[0][0];
		k=0;
		i1=0;
		for l in range(len(chID1)-1):
			while(chid==s1[k][0]):
				k+=1;
			ra1.append(k);
			chID1_d[s1[k][0]]=i1+1;
			chid=s1[k][0];
			i1+=1;
		ra1.append(len(s1));
		s1=np.array(s1);
		a1=s1[:,1];
		a2=s1[:,2];
		a1=list(map(int,a1));
		a2=list(map(int,a2));
		[file1_name,file1_extension]=file1.split('.');
		f2=open(file1_name+'_oc.'+file1_extension,'r');
		data1=f2.readlines();
		f2.close();
		ind=[];
		ind1=[];
		for i in range(len(chID)):
			if chID[i] in chID1_d:
				for l in range(ra1[chID1_d[chID[i]]],ra1[chID1_d[chID[i]]+1]):
					for k in range(ra[chID_d[chID[i]]],ra[chID_d[chID[i]]+1]):
						if((a1[l]<=b1[k]<=a2[l]) or (a1[l]<=b2[k]<=a2[l]) or (b1[k]<=a1[l]<=b2[k]) or (b1[k]<=a2[l]<=b2[k])):
							ind.append(k);
		q2=np.unique(ind);
		c3=np.zeros(len(s));
		j2=0;
		j3=0;
		for i in range(len(data1)):
			if(j2<len(q2)):
				if(q2[j2]==i):
					c3[i]=1
					j2+=1
		
			tx='\t'+str(c3[i])+'\n'
			data1[i]=data1[i].replace('\n',tx);
		f3=open(file1_name+'_oc.'+file1_extension,'w');
		data2=''.join(data1);
		f3.write(data2);
		f3.close();
		f4=open(file1_name+'_ocd.'+file1_extension,'r');
		data3=f4.readlines();
		f4.close();
		ind=[];
		ind1=[];
		for i in range(len(chID)):
			if chID[i] in chID1_d:
				for l in range(ra1[chID1_d[chID[i]]],ra1[chID1_d[chID[i]]+1]):
					for k in range(ra[chID_d[chID[i]]],ra[chID_d[chID[i]]+1]):
						if((a1[l]<=b1[k]<=a2[l]) or (a1[l]<=b2[k]<=a2[l]) or (b1[k]<=a1[l]<=b2[k]) or (b1[k]<=a2[l]<=b2[k])):
							ind.append(k);
		q2, q2c = np.unique(ind, return_counts=True)
		c3=np.zeros(len(s));
		j2=0;
		j3=0;
		for i in range(len(data3)):
			if(j2<len(q2)):
				if(q2[j2]==i):
					if q2c[j2]>=2:
						c3[i]=2
					else:
						c3[i]=1
					j2+=1
			tx='\t'+str(c3[i])+'\n'
			data3[i]=data3[i].replace('\n',tx);
		f5=open(file1_name+'_ocd.'+file1_extension,'w');
		data4=''.join(data3);
		f5.write(data4);
		f5.close();
	f6=open('tf_chipname.txt','w');
	data5='\n'.join(chipf_name);
	f6.write(data5);
	f6.close();
	s=[];s1=[];data1=[];data2=[];data3=[];data4=[];


print "Creating bipartite graph between genomic fragments and TFs ..."
chromatinfile=sys.argv[1]
chipseq_folder=sys.argv[2]
[chromatinfile_name,chromatinfile_extension]=chromatinfile.split('.')
chromatinfile18=chromatinfile_name+'_intra_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge.'+chromatinfile_extension

[chromatinfile18_name,chromatinfile18_extension]=chromatinfile18.split('.')
command1= 'cp '+chromatinfile18+' '+chromatinfile18_name+'_oc.'+chromatinfile18_extension
command2= 'cp '+chromatinfile18+' '+chromatinfile18_name+'_ocd.'+chromatinfile18_extension

os.system(command1)
os.system(command2)

bipartite_info(chromatinfile18,chipseq_folder)



