# TF_colocalization_3Dchromatin
Pipeline to find colocalizing TF pairs in 3D chromatin

This is a pipeline to find transcription factor(TF) pairs co-localizing in spatially proximal chromatin regions. Also, the pipeline finds the same co-localizing TF pairs in sequential contiguous i.e. on linear genome.

### Getting Started

The pipelines also requires packages bedtools, MEME Suite and various python libraries i.e. numpy, networkx, matplotlib etc. preinstalled in the system.

### Running the pipeline

	Usage: ./pipeline.sh -a file1 -b genomefastafile -c chipseq_directory -d pwm_directory -e number_of_randomizations
		-a chromatin fragment interaction file
		-b genome sequence fasta file
		-c directory containing chip-seq files and is optional if you want analyse binding sites
		-d directory containing pwm files of TFs and is optional if you want analyse motif sites
		-e number of randomization to perform
  
  
### Input files  

The pipeline requires the following input files:

1)Spatial proximal regions from any of the CHIA-PET, Capture Hi-C or HI-C experiments and the data should be in tab delimited with the columns i.e. chromosome of fragment1, start position of fragment1, stop position of fragment1, chromosome of fragment2, start position of fragment2, stop position of fragment2.  

2)Chip-seq peak files for TFs you want analyze the co-occurrence for the same cell line/tissue used for above chromatin interaction data. The chip-seq files should be in narrowPeak format with the first three columns chromosome ID, start position of the peak, end position of the peak. Name all the chip-seq files as the TF name with narrowPeak as the extension. For example name of chip-seq file for ATF1 is "ATF1.narrowPeak". Store all the renamed chip-seq files in a single directory.

3)PWMs files in meme format for the TFs you want analyze co-occurrence for motif sites. Name the files with TF name i.e. for "ATF1" the PWM file is "ATF1.meme" and store all the files in a single directory.

4) Reference genome sequence fasta file. For example data set analysis, you can download the human reference genome(hg19) from the following link https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

### Output files

The pipeline generates multiple files, but some important files are
1) qv_chip_3D.dat, qv_chip_3D.png gives the q-values for each pair of TFs for co-occurrence of binding peaks in 3D spatial proximal regions in data and heatmap respectively.
2) qv_chip_1D.dat, qv_chip_1D.png gives the q-values for each pair of TFs for co-occurrence of binding peaks in sequential contiguous regions in data and heatmap respectively.
3) qv_pwm_3D.dat, qv_pwm_3D.png gives the q-values for each pair of TFs for co-occurrence of motif sites in 3D spatial proximal regions in data and heatmap respectively.
4) qv_pwm_1D.dat, qv_pwm_1D.png gives the q-values for each pair of TFs for co-occurrence of motif sites in sequential contiguous regions in data and heatmap respectively.

### Running the example

The directory "example" contains example_interactions.bed which contains 3D spatial interactions between regions, chipseq folder with few chip-seq files, and pwm folder with few PWM files.

Execute the following command:

	./pipeline.sh -a example/example_interactions.bed -b "path of the genome sequence file" -c chipseq -d pwm -e 500





