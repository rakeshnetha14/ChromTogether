# TF_colocalization_3Dchromatin
Pipeline to find colocalizing TF pairs in 3D chromatin

This is a pipeline to find transcription factor(TF) pairs co-localizing in spatially proximal chromatin regions. Also, the pipeline finds the same co-localizing TF pairs in sequential contiguous i.e. on linear genome.

#Input: The analysis requires the following data
#1)Spatial proximal regions from any of the CHIA-PET, Capture Hi-C or HI-C experiments and the data should be in tab delimited with the following columns i.e. chromosome of fragment1, start position of fragment1, stop position of fragment1, chromosome of fragment2, start position of fragment2, stop position of fragment2.  
#2)Various transcription factor binding data from ChIP-Seq experiment for a given cell line. The chip-seq files can be in narrowPeak format with the first columns as chromosome ID, statt position of the peak, end position of the peak. Name all the chip-seq files as the TF name with narrowPeak as the extension. For example name of chip-seq file for ATF1 is "ATF1.narrowPeak". Store all the renamed chip-seq files in a directory "CHIP_SEQ_FILES".
