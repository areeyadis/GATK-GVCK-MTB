sample_name=$1

########## Software path ###############

prinseq_path1="../pipeline_software/prinseq-lite-0.20.4/prinseq-lite.pl"
prinseq_path2="../pipeline_software/prinseq-lite-0.20.4/prinseq-graphs-noPCA.pl"
bwa_path="../bwa-0.7.15/bwa mem"
GenomeAnalysisTK_path="../GenomeAnalysisTK.jar"
snpEff_path="../snpEff/snpEff/snpEff.jar"
GATK4_path="../gatk-4.1.0.0/gatk"

########## Reference ###############

ref_fasta="../index2/NC_000962.3_1.fasta"
ExClus_Bed="../index2/ExCluster.bed"

########## Output path ###############

output_path1="result/qc_trimmed/gd/"
output_path2="result/qc_trimmed/html/"
output_prefix1="result/"
out_QC=$output_prefix1"fastQC/"
out_map=$output_prefix1"bwa_map/"
out_VC_samtool=$output_prefix1"samtools/"
out_VC_GATK=$output_prefix1"GATK/"
out_trim=$output_prefix1"trimmed/"
out_annotate=$output_prefix1"annotate/"
out_stat=$output_prefix1"stat/"
out_ExClus=$output_prefix1"ExClus/"
out_ExClus_indel=$output_prefix1"ExClus_indel/"
output_prefix2="result/snp_indel/"
out_snpcall=$output_prefix2"snp/"
out_indelcall=$output_prefix2"indel/"
output_prefix3="result/filtered/"
out_filtersnp=$output_prefix3"snp/"
out_filterindel=$output_prefix3"indel/"
output_prefix4="result/pass_selected/"
out_passsnp=$output_prefix4"snp/"
out_passindel=$output_prefix4"indel/"


################ step 1 Trim ####################

	${prinseq_path1} -verbose -fastq ${sample_name}_1.fastq -fastq2 ${sample_name}_2.fastq -out_format 3 -ns_max_n 0 -derep 12345 -trim_qual_left 25 -trim_qual_right 25 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 3 -trim_qual_step 1 -min_qual_mean 25 -min_len 75 -lc_method entropy -lc_threshold 70 -out_good ${out_trim}${sample_name}_good -out_bad ${out_trim}${sample_name}_bad

#### QCplot in loop ###############

	${prinseq_path1} -verbose -fastq ${out_trim}${sample_name}_good_1.fastq -fastq2 ${out_trim}${sample_name}_good_2.fastq -graph_data ${output_path1}${sample_name}.gd -out_good null -out_bad null

###### html generation ############

	${prinseq_path2} -verbose -i ${output_path1}${sample_name}.gd -html_all -o ${output_path2}${sample_name}


############# Codes for mapping, realignment and variants calling ###########


#1.Mapping the raw read to H37Rv using BWA MEM
	${bwa_path} $ref_fasta ${out_trim}${sample_name}_good_1.fastq ${out_trim}${sample_name}_good_2.fastq > ${out_map}${sample_name}.sam -R '@RG\tID:'Illumina'\tSM:'${sample_name}'\tSW:bwa'

#2.Re-organizing the mapped read using samtools

#2.1 Convert sam file to bam file [compress file to smaller size]
	samtools view -bS ${out_map}${sample_name}.sam -o ${out_map}${sample_name}.bam

#2.2 Sorting the bam file [facilitate retrieval: faster resding time]
	samtools sort ${out_map}${sample_name}.bam ${out_map}${sample_name}.sort

#2.3 Indexing the sorted bam file
	samtools index ${out_map}${sample_name}.sort.bam

#2.4 Remove duplicate reads
	java -jar ../picard.jar MarkDuplicates I=${out_map}${sample_name}.sort.bam O=${out_map}${sample_name}_RmDup.bam M=${out_map}${sample_name}_RmDup.txt

#2.5 Indexing the RmDup bam file
	samtools index ${out_map}${sample_name}_RmDup.bam

#3.Realignment using GATK

#3.1 Create the realigner target
	java -Xmx64g -jar ${GenomeAnalysisTK_path} -T RealignerTargetCreator -R $ref_fasta -I ${out_map}${sample_name}_RmDup.bam -o ${out_map}${sample_name}_RmDup.intervals

#3.2 Do realignment
	java -Xmx64g -jar ${GenomeAnalysisTK_path} -T IndelRealigner -R $ref_fasta -I ${out_map}${sample_name}_RmDup.bam -targetIntervals ${out_map}${sample_name}_RmDup.intervals -o ${out_map}${sample_name}_RmDup.realn.bam -rf DuplicateRead -rf MappingQualityUnavailable

############# OPTIONAL STAT ANALYSIS ##################

	java -Xmx64g -jar ${GenomeAnalysisTK_path} -T DepthOfCoverage -R $ref_fasta -I ${out_map}${sample_name}.sort.bam.realn.bam -o ${out_stat}${sample_name}.report

########################################################

#4.Variants calling

#4.1 Variants calling using GATK:GVCF
	${GATK4_path} --java-options "-Xmx4g" HaplotypeCaller -R $ref_fasta -I ${out_map}${sample_name}_RmDup.realn.bam -O ${out_VC_GATK}${sample_name}.raw.g.vcf -ERC GVCF -ploidy 1 -mbq 20

#5.UnionGVCFs using GenomicsDBImport
	${GATK4_path} --java-options "-Xmx4g -Xms4g" GenomicsDBImport -V ${out_VC_GATK}${sample_name}.raw.g.vcf --genomicsdb-workspace-path my_database/ -L NC_000962

#6.Convert GVCF to vcf file (Union): GenotpyeGVCFs
	${GATK4_path} --java-options "-Xmx4g" GenotypeGVCFs -R $ref_fasta -V gendb://my_database -O ${out_VC_GATK}${sample_name}.vcf

#7.Variants annontation using snpEff
	java -jar ${snpEff_path} -ud 0 -formatEff -v m_tuberculosis_NC ${out_VC_GATK}${sample_name}.vcf > ${out_annotate}${sample_name}.ann.vcf -s ${out_annotate}${sample_name}.ann.html

#8. SNPs indel Select 
	#8.1 SNPs
		java -Xmx64g -jar ${GenomeAnalysisTK_path} -T SelectVariants -R $ref_fasta -V ${out_annotate}${sample_name}.ann.vcf -selectType SNP -o ${out_snpcall}${sample_name}_snp.vcf
	#8.2 Indel
		java -Xmx64g -jar ${GenomeAnalysisTK_path} -T SelectVariants -R $ref_fasta -V ${out_annotate}${sample_name}.ann.vcf -selectType INDEL -o ${out_indelcall}${sample_name}_indel.vcf

#9. SNPs indel filter
	#9.1 SNPs
		java -Xmx64g -jar ${GenomeAnalysisTK_path} -T VariantFiltration -R $ref_fasta -V ${out_snpcall}${sample_name}_snp.vcf --filterExpression "QD < 2.0 || DP < 10 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -cluster 3 -window 50 --filterName "Failfilter" -o ${out_filtersnp}${sample_name}_filteredsnp.vcf
	#9.2 Indel
		java -Xmx64g -jar ${GenomeAnalysisTK_path} -T VariantFiltration -R $ref_fasta -V ${out_indelcall}${sample_name}_indel.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "Failfilter" -o ${out_filterindel}${sample_name}_filteredindel.vcf

#10. Passed_call
	#10.1 SNPs
		java -Xmx64g -jar ${GenomeAnalysisTK_path} -T SelectVariants -R $ref_fasta -V ${out_filtersnp}${sample_name}_filteredsnp.vcf -o ${out_passsnp}${sample_name}_filteredsnp.pass.vcf -env -ef
	#10.2 Indels
		java -Xmx64g -jar ${GenomeAnalysisTK_path} -T SelectVariants -R -R $ref_fasta -V ${out_filterindel}${sample_name}_filteredindel.vcf -o ${out_passindel}${sample_name}_filteredindel.pass.vcf -env -ef

#11. Exclude PE/PPE regions, repeat regions, mobile elements and drug resistance genes
	vcftools --vcf ${out_passsnp}${sample_name}_filteredsnp.pass.vcf --exclude-bed ${ExClust_Bed} --recode --recode-INFO-all --out ${out_ExClus}${sample_name}
	vcftools --vcf ${out_passindel}${sample_name}_filteredindel.pass.vcf --exclude-bed ${ExClus_Bed} --recode --recode-INFO-all --out ${out_ExClus_indel}${sample_name}_indel
