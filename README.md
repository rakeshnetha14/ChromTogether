# TF_colocalization_3Dchromatin
Pipeline to find colocalizing TF pairs in 3D chromatin

This is a pipeline to find transcription factor(TF) pairs co-localizing in spatially proximal chromatin regions. Also, the pipeline finds the same co-localizing TF pairs in sequential contiguous i.e. on linear genome.

The pipeline requires the following input files:

1)Spatial proximal regions from any of the CHIA-PET, Capture Hi-C or HI-C experiments and the data should be in tab delimited with the columns i.e. chromosome of fragment1, start position of fragment1, stop position of fragment1, chromosome of fragment2, start position of fragment2, stop position of fragment2.  
2)Chip-seq peaks files available for various transcription factors for the same cell line/tissue used for above chromatin interaction data. The chip-seq files should be in narrowPeak format with the first three columns chromosome ID, start position of the peak, end position of the peak. Name all the chip-seq files as the TF name with narrowPeak as the extension. For example name of chip-seq file for ATF1 is "ATF1.narrowPeak". Store all the renamed chip-seq files in a single directory.

The pipelines also requires packages bedtools, MEME Suite and various python libraries i.e. numpy, networkx, matplotlib etc. preinstalled in the system.  
