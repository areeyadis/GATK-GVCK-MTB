# GATK-GVCK-MTB

Purpose:
Workflow for calling variants from whole genome sequencing data of Mycobacterium tuberculosis with GATK4.

Tool:
The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool from GATK4 in GVCF mode on a single sample according to GATK Best Practices. When executed the workflow scatters the HaplotypeCaller tool over a sample using an intervals list file. The output file produced will be a single gvcf file which can be used by the joint-discovery workflow.

Inputs:
    FastQ files

Outputs:
    One GVCF file and its index file

Intermediate files:
    SAM files
    BAM files and their index files
    Intervals list files

Exculding SNPs list:
    ExCluster.bed

Software version requirements:
    GATK version 4.1.2.0 
    BWA version 0.7.15
    SAMtools version 1.4
    SnpEff version 4.2
