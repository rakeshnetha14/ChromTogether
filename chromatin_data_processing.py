import os
import sys
import numpy as np
import networkx as nx

def individual_fragments(IntData):
	"""
	This function is to create indidual fragment data file format from chromatin interaction data file, which is used for sorting and merging.
	"""
	f=open(IntData,'r');
	s=[line.split('\t') for line in f]
	f.close()
	
	tx=[];
	for i in range(len(s)):
		ss=str(s[i][0])+'\t'+str(s[i][1])+'\t'+str(s[i][2])+'\n'+str(s[i][0])+'\t'+str(s[i][3])+'\t'+str(s[i][4]).strip('\n')+'\n';
		tx.append(ss);
	
	[IntData_name,IntData_extension]=IntData.split('.')
	f1=open(IntData_name+'_individual_fragments.'+IntData_extension,'w');
	data=''.join(tx);
	f1.write(data);
	f1.close();
	s=[];tx=[];data=[];

	#return f1


def merging_fragments(file1,file2):
	"""
	This function corrects the coordinates in chromatin interaction data file to the corresponding overlapping regions after merging the overlapping fragments into one fragment
	"""
	f2=open(file1,'r');
	s2=[line.split('\t') for line in f2];
	f2.close()
	ra=[0];
	chid=s2[0][0];
	k=0;
	for l in range(22):
		while(chid==s2[k][0]):
			k+=1;
		ra.append(k);
		chid=s2[k][0];
	ra.append(len(s2));
	####
	f3=open(file2,'r');
	s3=[line.split('\t') for line in f3];
	f3.close()
	#####
	ra1=[0];
	chid=s3[0][0];
	k=0;
	for l in range(22):
		while(chid==s3[k][0]):
			k+=1;
		ra1.append(k);
		chid=s3[k][0];
	ra1.append(len(s3));
	#####
	for i in range(23):
		for k in range(ra[i],ra[i+1]):
			l1=ra1[i]
			for l in range(l1,ra1[i+1]):
				if (int(s3[l][1])<=int(s2[k][1]) and int(s2[k][2])<=int(s3[l][2])):
					s2[k][1]=str(int(s3[l][1]))
					s2[k][2]=str(int(s3[l][2]))
					break
				else:
					l1+=1
	for i in range(23):
		for k in range(ra[i],ra[i+1]):
			for l in range(ra1[i],ra1[i+1]):
				if (int(s3[l][1])<=int(s2[k][3]) and int(s2[k][4])<=int(s3[l][2])):
					s2[k][3]=str(int(s3[l][1]))
					s2[k][4]=s3[l][2]
					break
	#####
	tx=[];
	for i in range(len(s2)):
		ss=s2[i][0]+'\t'+s2[i][1]+'\t'+s2[i][2]+'\t'+s2[i][3]+'\t'+s2[i][4];
		tx.append(ss);
	[file1_name,file1_extension]=file1.split('.')
	f4=open(file1_name+'_aftermerging.'+file1_extension,'w');
	data1=''.join(tx);
	f4.write(data1);
	f4.close();
	s2=[];s3=[];tx=[];data1=[];

def duplicate_removal(file3):
	"""
	This function is to remove duplicate interactions if any present in the chromatin interaction data
	"""
	f5=open(file3,'r');
	s5=[line.split('\t') for line in f5];
	f5.close();
	######
	ra=[0];
	chid=s5[0][0];
	k=0;
	for l in range(22):
		while(chid==s5[k][0]):
			k+=1;
		ra.append(k);
		chid=s5[k][0];
	ra.append(len(s5));
	######
	ind=[];
	for k in range(23):
		for i in range(ra[k],ra[k+1]):
			for j in range(i+1,ra[k+1]):
				if(int(s5[i][1])==int(s5[j][1]) and int(s5[i][2])==int(s5[j][2]) and int(s5[i][3])==int(s5[j][3]) and int(s5[i][4])==int(s5[j][4])):
					ind.append(j);
	ind1=np.unique(ind);
	ind2=list(set(range(len(s5))).difference(set(ind1)));
	ind3=np.sort(ind2);
	#####
	f6=open(file3,'r');
	data2=f6.readlines();
	f6.close();
	#####
	data2=np.array(data2);
	data3=[];
	for i in range(len(ind3)):
		data3.append(data2[ind3[i]]);
	#####
	[file3_name,file3_extension]=file3.split('.')
	f7=open(file3_name+'_rrm.'+file3_extension,'w');
	data4=''.join(data3);
	f7.write(data4);
	f7.close();
	s5=[];data2=[];data3=[];data4=[];

def selfinteraction_removal(file4):
	"""
	This function is to remove self interactions if any present.
	"""
	f8=open(file4,'r');
	s8=[line.split('\t') for line in f8];
	f8.close();
	#####
	ind=[]
	for i in range(len(s8)):
		if int(s8[i][1])==int(s8[i][3]) and int(s8[i][2])==int(s8[i][4]):
			ind.append(i)
	ind2=list(set(range(len(s8))).difference(set(ind)));
	ind3=np.sort(ind2);
	#####
	f9=open(file4,'r');
	data5=f9.readlines();
	f9.close();
	data5=np.array(data5);
	data6=[];
	for i in range(len(ind3)):
		data6.append(data5[ind3[i]]);
	#####
	[file4_name,file4_extension]=file4.split('.')
	f10=open(file4_name+'_srm.'+file4_extension,'w');
	data7=''.join(data6);
	f10.write(data7);
	f10.close();
	s8=[];data5=[];data6=[];data7=[];

def shortrangeinteraction_removal(file5):
	"""
	This function is to filter out any interaction present within some cut-off distance of linear genome
	"""
	f11=open(file5,'r');
	s11=[line.split('\t') for line in f11];
	f11.close();
	f12=open(file5,'r');
	data8=f12.readlines();
	f12.close();
	data8=np.array(data8);
	data9=[];
	#####
	for i in range(len(s11)):
		so=np.sort([int(s11[i][1]),int(s11[i][2]),int(s11[i][3]),int(s11[i][4])])
		if (so[2]-so[1]>=2000):
			data9.append(data8[i])
	#####
	[file5_name,file5_extension]=file5.split('.')	
	f13=open(file5_name+'_2kbfil.'+file5_extension,'w');
	data10=''.join(data9);
	f13.write(data10);
	f13.close();
	s11=[];data8=[];data9=[];data10=[];

def makingnodefile(file6,file7):
	"""
	This function assigns each fragment in the chromatin interaction data file with the node number
	"""
	f14=open(file6,'r');
	d={}
	i=0;
	for line in f14:
		l=line.split('\t')
		name=l[0]+'-'+l[1]+'-'+str(int(l[2]))
		d[name]=i
		i+=1;
	f15=open(file7,'r');
	data11=f15.readlines();
	f15.close();
	f16=open(file7,'r');
	s16=[line.split('\t') for line in f16];
	f16.close();
	####
	for i in range(len(s16)):
		name1=s16[i][0]+'-'+s16[i][1]+'-'+str(int(s16[i][2]))
		name2=s16[i][0]+'-'+s16[i][3]+'-'+str(int(s16[i][4]))
		tx='\t'+str(d[name1])+'\t'+str(d[name2])+'\n'
		data11[i]=data11[i].replace('\n',tx);
	#####
	[file7_name,file7_extension]=file7.split('.')
	f17=open(file7_name+'_asnodes.'+file7_extension,'w');
	data12=''.join(data11);
	f17.write(data12);
	f17.close();
	data11=[];data12=[];

def highdegree_removal(file8,file8_1):
	"""
	This function is to filter out interactions containing fragment with node degree higher than some cut-off.
	"""
	f18=open(file8_1,'r');
	s18=[line.split('\t') for line in f18];
	f18.close()
	####
	realGedges=[(int(s18[i][5]),int(s18[i][6])) for i in range(len(s18))];
	realG = nx.Graph()
	realG.add_edges_from(realGedges)
	#print realG.number_of_nodes()
	####
	realdegree=list(realG.degree)
	realdegreeseq=[realdegree[i][1] for i in range(len(realdegree))]
	ddg={}
	for i in range(len(realdegreeseq)):
		ddg[i]=realdegreeseq[i]
	fil20=[];
	for i in range(len(s18)):
		if (ddg[realGedges[i][0]]>20 or ddg[realGedges[i][1]]>20):
			fil20.append(i)
	ind2=list(set(range(len(s18))).difference(set(fil20)));
	ind3=np.sort(ind2);
	#####
	f19=open(file8,'r');
	data13=f19.readlines();
	f19.close();
	data13=np.array(data13);
	data14=[];
	for i in range(len(ind3)):
		data14.append(data13[ind3[i]]);
	#####
	[file8_name,file8_extension]=file8.split('.')
	f20=open(file8_name+'_d20fil.'+file8_extension,'w');
	data15=''.join(data14);
	f20.write(data15);
	f20.close();
	s18=[];data13=[];data14=[];data15=[];

def highlength_removal(file9):
	"""
	This function is to remove interactions containing fragments whose length is higher than some cut-off.
	"""
	f21=open(file9,'r');
	s21=[line.split('\t') for line in f21];
	f21.close();
	f22=open(file9,'r');
	data16=f22.readlines();
	f22.close();
	data16=np.array(data16);
	data17=[];
	for i in range(len(s21)):
		if (int(s21[i][2])-int(s21[i][1])<30000 and int(s21[i][4])-int(s21[i][3])<30000):
			data17.append(data16[i])
	[file9_name,file9_extension]=file9.split('.')	
	f23=open(file9_name+'_30kblenfil.'+file9_extension,'w');
	data18=''.join(data17);
	f23.write(data18);
	f23.close();
	s21=[];data16=[];data17=[];data18=[];



chromatinfile=sys.argv[1]
[chromatinfile_name,chromatinfile_extension]=chromatinfile.split('.')
chromatinfile1=chromatinfile_name+'_sort.'+chromatinfile_extension
chromatinfile2=chromatinfile_name+'_individual_fragments.'+chromatinfile_extension
chromatinfile3=chromatinfile_name+'_individual_fragments_sort.'+chromatinfile_extension
chromatinfile4=chromatinfile_name+'_individual_fragments_sort_merge.'+chromatinfile_extension
chromatinfile5=chromatinfile_name+'_sort_aftermerging.'+chromatinfile_extension
chromatinfile6=chromatinfile_name+'_sort_aftermerging_sort.'+chromatinfile_extension
chromatinfile7=chromatinfile_name+'_sort_aftermerging_sort_rrm.'+chromatinfile_extension
chromatinfile8=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm.'+chromatinfile_extension
chromatinfile9=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil.'+chromatinfile_extension
chromatinfile10=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_individual_fragments.'+chromatinfile_extension
chromatinfile11=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_individual_fragments_sort.'+chromatinfile_extension
chromatinfile12=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_individual_fragments_sort_merge.'+chromatinfile_extension
chromatinfile13=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_asnodes.'+chromatinfile_extension
chromatinfile14=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil.'+chromatinfile_extension
chromatinfile15=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil.'+chromatinfile_extension
chromatinfile16=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments.'+chromatinfile_extension
chromatinfile17=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort.'+chromatinfile_extension
chromatinfile18=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge.'+chromatinfile_extension


command1 = 'bedtools sort -i'+' '+chromatinfile+' > '+chromatinfile1

os.system(command1)

individual_fragments(chromatinfile)

command2= 'bedtools sort -i'+' '+chromatinfile2+' > '+chromatinfile3
command3= 'bedtools merge -i'+' '+chromatinfile3+' > '+chromatinfile4
command4= 'bedtools sort -i'+' '+chromatinfile5+' > '+chromatinfile6
command5= 'bedtools sort -i'+' '+chromatinfile10+' > '+chromatinfile11
command6= 'bedtools merge -i'+' '+chromatinfile11+' > '+chromatinfile12
command7= 'bedtools sort -i'+' '+chromatinfile16+' > '+chromatinfile17
command8= 'bedtools merge -i'+' '+chromatinfile17+' > '+chromatinfile18

os.system(command2)
os.system(command3)
merging_fragments(chromatinfile1,chromatinfile4)
os.system(command4)
duplicate_removal(chromatinfile6)
selfinteraction_removal(chromatinfile7)
shortrangeinteraction_removal(chromatinfile8)
individual_fragments(chromatinfile9)
os.system(command5)
os.system(command6)
makingnodefile(chromatinfile12,chromatinfile9)
highdegree_removal(chromatinfile9,chromatinfile13)
highlength_removal(chromatinfile14)
individual_fragments(chromatinfile15)
os.system(command7)
os.system(command8)
makingnodefile(chromatinfile18,chromatinfile15)



