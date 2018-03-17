hg38 = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/Homo_sapiens_assembly38
hg38Index = /data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/somatics/hg38-index
STAR = /data/WHRI-EndocrinePituitaryGroup/software/STAR/bin/Linux_x86_64
PICARD = /data/WHRI-EndocrinePituitaryGroup/software/picard
GATK= /data/WHRI-EndocrinePituitaryGroup/software/gatk
KNOWNINDELS = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/Homo_sapiens_assembly38.known_indels.vcf
DBSNP = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/Homo_sapiens_assembly38.dbsnp138.vcf
GOLDINDELS = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf
WHERE = /data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/somatics/sp1452

CURRENT = sp1452
s1_1 = SP1452_S3_R1_001_val_1
s1_2 = SP1452_S3_R2_001_val_2

star1: 
	${STAR}/STAR --runMode genomeGenerate --genomeDir ${hg38Index} --genomeFastaFiles ${hg38}.fasta  --runThreadN 4

sp1452Aligned.out.sam:
	gunzip -c ${WHERE}/${s1_1}.fq.gz > ${WHERE}/${s1_1}.fq 
	gunzip -c ${WHERE}/${s1_2}.fq.gz > ${WHERE}/${s1_2}.fq
	${STAR}/STAR --genomeDir ${hg38Index} --readFilesIn ${WHERE}/${s1_1}.fq ${WHERE}/${s1_2}.fq --runThreadN 2 --outFileNamePrefix ${WHERE}/${CURRENT}

sp1452sjdbList.out.tab:
	${STAR}/STAR --runMode genomeGenerate --genomeDir ${WHERE} --genomeFastaFiles  ${hg38}.fasta \
	--sjdbFileChrStartEnd ${WHERE}/${CURRENT}SJ.out.tab --sjdbOverhang 75 --runThreadN 2 --outFileNamePrefix ${WHERE} 


Aligned.out.sam:
	${STAR}/STAR --genomeDir ${WHERE} --readFilesIn ${WHERE}/${s1_1}.fq ${WHERE}/${s1_2}.fq --runThreadN 2 


sp1452.dedupped.bam: 
	java -Dlog4j.configurationFile="log4j2.xml" -jar ${PICARD}/picard.jar AddOrReplaceReadGroups I=${WHERE}/Aligned.out.sam O=${WHERE}/rg_added_sorted.bam \
	SO=coordinate RGID=@E00461_116.SP1452_S3 RGLB=SP1452_S3 RGPL=ILLUMINA RGPU=@E00461_116_GW170602261_1.SP1452_S3 RGSM=SP1452
	java -Dlog4j.configurationFile="log4j2.xml" -jar ${PICARD}/picard.jar MarkDuplicates I=${WHERE}/rg_added_sorted.bam O=${WHERE}/${CURRENT}.dedupped.bam  CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT M=${WHERE}/output.metrics 
	
sp1452.split.bam: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${hg38}.fasta -I ${WHERE}/${CURRENT}.dedupped.bam -o ${WHERE}/${CURRENT}.split.bam \
	 -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

sp1452realignertargetcreator.intervals:
	java -jar ${GATK}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${hg38}.fasta \
	-known ${KNOWNINDELS} -I ${WHERE}/${CURRENT}.split.bam -o ${WHERE}/${CURRENT}.realignertargetcreator.intervals  

sp1452-realigned.bam: 
	java -Xmx50G -Djava.io.tmpdir=${WHERE}/tmp -jar ${GATK}/GenomeAnalysisTK.jar -T IndelRealigner -R ${hg38}.fasta \
	-targetIntervals ${WHERE}/${CURRENT}.realignertargetcreator.intervals -known ${KNOWNINDELS} -I ${WHERE}/${CURRENT}.split.bam \
	-o ${WHERE}/${CURRENT}-realigned.bam 

#Changed to skip indel realignment 
sp1452.recal_data.table: 
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-R ${hg38}.fasta \
	-I ${WHERE}/${CURRENT}.split.bam \
	-knownSites ${DBSNP} \
	-knownSites ${GOLDINDELS} \
	-o ${WHERE}/${CURRENT}.recal_data.table 

sp1452.post_recal_data.table: 
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-R ${hg38}.fasta \
	-I ${WHERE}/${CURRENT}.split.bam \
	-knownSites ${DBSNP}  \
	-knownSites ${GOLDINDELS} \
	-BQSR ${WHERE}/${CURRENT}.recal_data.table \
	-o ${WHERE}/${CURRENT}.post_recal_data.table 

sp1452.recabliration_plots.pdf: 
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T AnalyzeCovariates \
	-R ${hg38}.fasta \
	-before ${WHERE}/${CURRENT}.recal_data.table \
	-after ${WHERE}/${CURRENT}.post_recal_data.table \
	-plots ${WHERE}/${CURRENT}.recalibration_plots.pdf

#Assuming post recabliration is better than before recabliration 

sp1452.bam: 
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T PrintReads \
	-R ${hg38}.fasta  \
	-I  ${WHERE}/${CURRENT}.split.bam \
	-BQSR ${WHERE}/${CURRENT}.post_recal_data.table \
	-o ${WHERE}/${CURRENT}.bam
	#This is recablirated bam 


sp1452.vcf: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${hg38}.fasta -I ${WHERE}/${CURRENT}-recablirated.bam \
	-dontUseSoftClippedBases -stand_call_conf 20.0 -o ${CURRENT}.vcf


 
