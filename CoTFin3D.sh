#!/bin/sh

# This is a pipeline to find the significantly attracting and repelling TF pairs in the spatial proximal regions of 
# the chromatin data. It can also find co-occurrence of TF pairs in sequentailly contiguous regions of the genome.

helpFunction()
{
   echo ""
   echo "Usage: $0 -a file1 -b genomefastafile -c chipseq_directory -d pwm_directory -e number_of_randomizations"
   echo -e "\t-a chromatin fragment interaction file"
   echo -e "\t-b genome sequence fasta file"
   echo -e "\t-c directory containing chip-seq files and is optional if you want analyse binding sites"
   echo -e "\t-d directory containing pwm files of TFs and is optional if you want analyse motif sites"
   echo -e "\t-e number of randomization to perform"
   exit 1
}

while getopts "a:b:c:d:e:" opt
do
   case "$opt" in
      a ) file1="$OPTARG" ;;
      b ) genomefastafile="$OPTARG" ;;
      c ) chipseq_directory="$OPTARG" ;;
      d ) pwm_directory="$OPTARG" ;;
      e ) number_of_randomizations="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

if [ -z "$file1" ] || [ -z "$genomefastafile" ] || [ -z "$number_of_randomizations" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

if [ -z "$chipseq_directory" ] && [ -z "$pwm_directory" ]
then
   echo "both directories are missing";
   helpFunction
fi

python chromatin_data_processing.py $file1

if [ -z "$chipseq_directory" ] || [ -z "$pwm_directory" ]
then
   if [ ! -z "$chipseq_directory" ]
   then
      python chipseq_data_processing.py $file1 $chipseq_directory
      python bipartite_randomization.py $file1 $number_of_randomizations 0
   else 
      python motifsites.py $file1 $genomefastafile $pwm_directory
      python bipartite_randomization.py $file1 $number_of_randomizations 1
   fi
else
   python chipseq_data_processing.py $file1 $chipseq_directory
   python motifsites.py $file1 $genomefastafile $pwm_directory
   python bipartite_randomization.py $file1 $number_of_randomizations 2
fi


