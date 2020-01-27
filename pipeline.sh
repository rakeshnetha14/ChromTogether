#!/bin/sh

# This is a pipeline to find the significantly attracting and repelling TF pairs in the spatial proximal regions of the chromatin data. Also, it contain also find the above TF pairs in sequentailly in the above considered regions of chromatin.

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


