# GATK-GVCK-MTB

Purpose:
Workflow for calling variants from whole genome sequencing data of Mycobacterium tuberculosis.

Tool:
The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool from GATK4 in GVCF mode on a single sample according to GATK Best Practices. When executed the workflow scatters the HaplotypeCaller tool over a sample using an intervals list file. The output file produced will be a single gvcf file which can be used by the joint-discovery workflow.
