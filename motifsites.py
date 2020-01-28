import os
import sys
import re
import numpy as np

def motif_scanning(folder1,file1,file2):
	"""
	This function scans the fragments regions for the presence of motif sites for PWM matrix of various TFs.
	"""
	motf=os.listdir(folder1)
	motf.sort()
	tf_name=[]
	for i in range(len(motf)):
		[name,extension]=motf[i].split('.');
		tf_name.append(name);
		command= 'fimo -o '+name+' '+folder1+'/'+motf[i]+' '+file1
		os.system(command)
		f=open(name+'/fimo.tsv','r')
		s=[line.split('\t')[2] for line in f if len(line.split('\t'))>=2];
		s1=[re.split(',|:|-|\..',s[j]) for j in range(1,len(s))];
		f.close();
		mind=[]
		for k in range(len(s)-1):
			mind.append(int(s1[k][3])-1)
		mind1=np.sort(np.unique(mind))
		[file2_name,file2_extension]=file2.split('.');
		f1=open(file2_name+'_mc.'+file2_extension,'r')
		data1=f1.readlines();
		f1.close();
		mind2=list(set(range(len(data1)))-set(mind1))
		for l in range(len(mind1)):
			tx='\t'+'1'+'\n'
			data1[mind1[l]]=data1[mind1[l]].replace('\n',tx);
		for m in range(len(mind2)):
			tx='\t'+'0'+'\n'
			data1[mind2[m]]=data1[mind2[m]].replace('\n',tx);
		f2=open(file2_name+'_mc.'+file2_extension,'w')
		data2=''.join(data1);
		f2.write(data2);
		f2.close();
		mind3,mind3c=np.unique(mind, return_counts=True)
		f3=open(file2_name+'_mcd.'+file2_extension,'r')
		data3=f3.readlines();
		f3.close();
		mind4=list(set(range(len(data3)))-set(mind3))
		for l in range(len(mind3)):
			if mind3c[l]>=2:
				tx='\t'+'2'+'\n'
				data3[mind3[l]]=data3[mind3[l]].replace('\n',tx);
			else:
				tx='\t'+'1'+'\n'
				data3[mind3[l]]=data3[mind3[l]].replace('\n',tx);
		for m in range(len(mind4)):
			tx='\t'+'0'+'\n'
			data3[mind4[m]]=data3[mind4[m]].replace('\n',tx);
		f4=open(file2_name+'_mcd.'+file2_extension,'w')
		data4=''.join(data3);
		f4.write(data4);
		f4.close();
	f6=open('tf_pwmname.txt','w');
	data5='\n'.join(tf_name);
	f6.write(data5);
	f6.close();
	s=[];s1=[];data1=[];data2=[];data3=[];data4=[];

def rename_fa(file1):
	f=open(file1)
	ff=open(file1+'addname','w')
	i=0
	count=1
	for line in f:
		if i%2==0:
			l=line.split('\n')
			ff.write(l[0]+':'+str(count)+'\n') 
			count=count+1 
		else:
			ff.write(line)
		i=i+1


chromatinfile=sys.argv[1]
genomefasta=sys.argv[2]
pwm_folder=sys.argv[3]
[chromatinfile_name,chromatinfile_extension]=chromatinfile.split('.')
chromatinfile18=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge.'+chromatinfile_extension
chromatinfile24=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge.fa'
chromatinfile25=chromatinfile_name+'_sort_aftermerging_sort_rrm_srm_2kbfil_d20fil_30kblenfil_individual_fragments_sort_merge.faaddname'


[chromatinfile18_name,chromatinfile18_extension]=chromatinfile18.split('.')
command1= 'cp '+chromatinfile18+' '+chromatinfile18_name+'_mc.'+chromatinfile18_extension
command2= 'cp '+chromatinfile18+' '+chromatinfile18_name+'_mcd.'+chromatinfile18_extension
command3= 'bedtools getfasta -fi '+genomefasta+' -fo '+chromatinfile24+' -bed '+chromatinfile18
os.system(command1)
os.system(command2)
os.system(command3)

rename_fa(chromatinfile24)
motif_scanning(pwm_folder,chromatinfile25,chromatinfile18)


