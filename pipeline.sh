#!/bin/sh

# This is a pipeline to find the significantly attracting and repelling TF pairs in the spatial proximal regions of the chromatin data. Also, it contain also find the above TF pairs in sequentailly in the above considered regions of chromatin.

read -p "Enter Chromatin interaction data file name : " file1
read -p "Enter the genome sequence fasta file" file2
read -p "Enter the directory containing chip-seq files" directory1
read -p "Enter the directory containing pwms" directory2
read -p "Enter the number of randomizations" nr

python chromatin_data_processing.py $file1
python chipseq_data_processing.py $file1 $directory1
python motifsites.py $file1 $file2 $directory2
python bipartite_randomization.py $file1 $nr


