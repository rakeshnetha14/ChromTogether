#!/bin/sh

# This is a pipeline to find the significantly attracting and repelling TF pairs in the spatial proximal regions of the chromatin data. Also, it contain also find the above TF pairs in sequentailly in the above considered regions of chromatin.

#Input: The analysis requires the following data
#1)Spatial proximal regions from any of the CHIA-PET, Capture Hi-C or HI-C experiments and the data should be in tab delimited with the following columns i.e. chromosome of fragment1, start position of fragment1, stop position of fragment1, chromosome of fragment2, start position of fragment2, stop position of fragment2.  
#2)Various transcription factor binding data from ChIP-Seq experiment for a given cell line. The chip-seq files can be in narrowPeak format with the first columns as chromosome ID, statt position of the peak, end position of the peak. Name all the chip-seq files as the TF name with narrowPeak as the extension. For example name of chip-seq file for ATF1 is "ATF1.narrowPeak". Store all the renamed chip-seq files in a directory "CHIP_SEQ_FILES".

read -p "Enter Chromatin interaction data file name : " file1
read -p "Enter the directory containing chip-seq files" directory1
read -p "Enter the number of randomizations" nr

#file1_name=$(echo $name | cut -f 1 -d '.')
#file1_extension=$(echo $name | cut -f 2 -d '.')
#file2=

#echo $file2

#bedtools sort -i $file1 > $file1_name'_sort.'$file1_extension

python chromatin_data_processing.py $file1
python chipseq_data_processing.py $file1 $directory1
python bipartite_randomization.py $file1 $nr


