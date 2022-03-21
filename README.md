%# TF_colocalization_3Dchromatin
# Co-occurring TF pair Finder in 3D genome (CoTFin3D)
Pipeline to find colocalizing TF pairs in 3D chromatin

This is a pipeline to find transcription factor(TF) pairs co-localizing in spatially proximal chromatin regions. Also, the pipeline can find the co-localizing TF pairs in sequential contiguous i.e. on linear genome too.

### Getting Started

The pipeline requires packages such as bedtools, MEME Suite 5.1.0 and Python 3.6 and python libraries i.e. numpy, networkx, matplotlib etc. preinstalled in the system.

### Running the pipeline

	Usage: ./CoTFin3D.sh -a <filename> -b <genomefastafile> -c <chipseq_directory> -d <pwm_directory> -e <number_of_randomizations>
		-a chromatin fragment interaction file
		-b genome sequence fasta file
		-c directory containing chip-seq files and is optional if you want analyse binding sites
		-d directory containing pwm files of TFs and is optional if you want analyse motif sites
		-e number of randomization to perform
  
  
### Input files  

The pipeline requires the following input files:

1)Spatially proximal chromatin region data from any of the CHIA-PET, Capture Hi-C or Hi-C experiments and the data should be in tab delimited with the 6 columns i.e. chromosome of fragment1, start position of fragment1, stop position of fragment1, chromosome of fragment2, start position of fragment2, stop position of fragment2.  

2)Chip-seq peak files of TFs that one want analyze for co-occurrence of the same cell line/tissue of the chromatin interaction dataset. The chip-seq files should be in narrowPeak format which contains chromosome ID, start position of the peak, end position of the peak as first three columns. Filename of the chip-seq file should be named with TF narrowPeak as the extension. For example name of chip-seq file for ATF1 should be "ATF1.narrowPeak". Store all the renamed chip-seq files in a single directory.

3)PWMs files in meme format for the TFs you want analyze for co-occurrence of motif sites. Name the files with TF name i.e. for "ATF1"  the PWM file is "ATF1.meme" and store all the files in a single directory.

4)Reference genome sequence fasta file. For example data set provided with the tool, please download the human reference genome(hg19) from the following link https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz and then unzip and place it in the example folder.

### Output files

The pipeline generates multiple files and some important files among them are
1) qv_chip_3D.dat, qv_chip_3D.png gives the q-values for each pair of TFs for co-occurrence of binding peaks in 3D spatial proximal regions in data and heatmap respectively.
2) qv_chip_1D.dat, qv_chip_1D.png gives the q-values for each pair of TFs for co-occurrence of binding peaks in sequential contiguous regions in data and heatmap respectively.
3) qv_pwm_3D.dat, qv_pwm_3D.png gives the q-values for each pair of TFs for co-occurrence of motif sites in 3D spatial proximal regions in data and heatmap respectively.
4) qv_pwm_1D.dat, qv_pwm_1D.png gives the q-values for each pair of TFs for co-occurrence of motif sites in sequential contiguous regions in data and heatmap respectively.

### Running the example

The directory "example" provided here contains example_interactions.bed which contains 3D spatial interactions between regions, chipseq folder with few chip-seq files, and pwm folder with few PWM files.

Execute the following command:

	./pipeline.sh -a example/example_interactions.bed -b "path of the genome sequence file" -c chipseq -d pwm -e 500

