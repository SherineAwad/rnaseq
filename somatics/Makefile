hg38 = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/Homo_sapiens_assembly38
hg38Index = /data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/somatics/hg38-index
STAR = /data/WHRI-EndocrinePituitaryGroup/software/STAR/bin/Linux_x86_64
PICARD = /data/WHRI-EndocrinePituitaryGroup/software/picard
WHERE = /data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/somatics
GATK4 = /data/WHRI-EndocrinePituitaryGroup/software/gatk-4.0.0.0
GATK = /data/WHRI-EndocrinePituitaryGroup/software/gatk
KNOWNINDELS = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/Homo_sapiens_assembly38.known_indels.vcf
DBSNP = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/Homo_sapiens_assembly38.dbsnp138.vcf
GOLDINDELS = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf
COSMIC = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/COSMIC/cosmic.vcf 
AFGNOME = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/hg38GATK/af-only-gnomad.hg38.vcf.gz
ANNOVAR = /data/WHRI-EndocrinePituitaryGroup/software/annovar
BBMAP = /data/WHRI-EndocrinePituitaryGroup/software/bbmap
#Using split.bam files temporary 
sp1522.bam: 
	ln -fs  ${WHERE}/sp1508/sp1508.split.bam sp1508.bam
	ln -fs  ${WHERE}/sp1510/sp1510.split.bam sp1510.bam
	ln -fs  ${WHERE}/sp1514/sp1514.split.bam sp1514.bam
	ln -fs  ${WHERE}/sp1515/sp1515.split.bam sp1515.bam
	ln -fs  ${WHERE}/sp1516/sp1516.split.bam sp1516.bam
	ln -fs  ${WHERE}/sp1517/sp1517.split.bam sp1517.bam
	ln -fs  ${WHERE}/sp1519/sp1519.split.bam sp1519.bam
	ln -fs  ${WHERE}/sp1520/sp1520.split.bam sp1520.bam
	ln -fs  ${WHERE}/sp1521/sp1521.split.bam sp1521.bam
	ln -fs  ${WHERE}/sp1522/sp1522.split.bam sp1522.bam

sp1508.vcf: sp1508.bam 
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1508.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1508.vcf

sp1510.vcf: sp1510.bam
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1510.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1510.vcf

sp1514.vcf: sp1514.bam
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1514.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1514.vcf

sp1515.vcf: sp1515.bam
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1515.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1515.vcf

sp1516.vcf: sp1516.bam
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1516.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1516.vcf

sp1517.vcf: sp1517.bam
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1517.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1517.vcf

sp1519.vcf: sp1519.bam
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1519.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1519.vcf

sp1520.vcf: sp1520.bam
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1520.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1520.vcf

sp1521.vcf: sp1521.bam	
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1521.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1521.vcf

sp1522.vcf: sp1522.bam
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${hg38}.fasta \
	-I:tumor sp1522.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--artifact_detection_mode \
	-o sp1522.vcf
pon.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T CombineVariants \
	--arg_file NORMALS.list \
	-minN 2 \
	--setKey "null" \
	--filteredAreUncalled \
	--filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
	-o pon.vcf.gz \
	-R ${hg38}.fasta

pon_sitesonly.vcf.gz:
	java -Dlog4j.configurationFile="log4j2.xml" -jar ${PICARD}/picard.jar MakeSitesOnlyVcf \
	I= pon.vcf.gz \
	O= pon_siteonly.vcf.gz	

sp1452mutect2.vcf.gz:  
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1452.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1452mutect2.vcf.gz \
	-R ${hg38}.fasta

sp1466mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1466.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1466mutect2.vcf.gz \
	-R ${hg38}.fasta

sp1470mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1470.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1470mutect2.vcf.gz \
	-R ${hg38}.fasta

sp1476mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1476.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1476mutect2.vcf.gz \
	-R ${hg38}.fasta

sp1477mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1477.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1477mutect2.vcf.gz \
	-R ${hg38}.fasta

sp1479mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1479.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1479mutect2.vcf.gz \
	-R ${hg38}.fasta


sp1481mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1481.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1481mutect2.vcf.gz \
	-R ${hg38}.fasta

sp1483mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1483.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1483mutect2.vcf.gz \
	-R ${hg38}.fasta



sp1485mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1485.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1485mutect2.vcf.gz \
	-R ${hg38}.fasta

sp1487mutect2.vcf.gz: 
	java -jar ${GATK}/GenomeAnalysisTK.jar -T MuTect2 \
	-I:tumor sp1487.bam \
	--dbsnp ${DBSNP} \
	--cosmic ${COSMIC} \
	--normal_panel pon_siteonly.vcf.gz \
	--output_mode EMIT_VARIANTS_ONLY \
	-o sp1487mutect2.vcf.gz \
	-R ${hg38}.fasta

sp1452_indels.vcf.gz: sp1452mutect2.vcf.gz
	java -jar ${GATK}/GenomeAnalysisTK.jar -T SelectVariants \
	-V sp1452mutect2.vcf.gz \
	-o mutect2_indels.vcf.gz \
	-selectType INDEL \
	-R ${hg38}.fastai \
	--list

sp1452.pass.recode.vcf: 
	vcftools --gzvcf sp1452mutect2.vcf.gz --remove-filtered-all --recode --out sp1452.pass

sp1466.pass.recode.vcf: 
	vcftools --gzvcf sp1466mutect2.vcf.gz --remove-filtered-all --recode --out sp1466.pass
sp1470.pass.recode.vcf: 
	vcftools --gzvcf sp1470mutect2.vcf.gz --remove-filtered-all --recode --out sp1470.pass 

sp1476.pass.recode.vcf:
	vcftools --gzvcf sp1476mutect2.vcf.gz --remove-filtered-all --recode --out sp1476.pass 

sp1477.pass.recode.vcf: 
	vcftools --gzvcf sp1477mutect2.vcf.gz --remove-filtered-all --recode --out sp1477.pass 

sp1479.pass.recode.vcf: 
	vcftools --gzvcf sp1479mutect2.vcf.gz --remove-filtered-all --recode --out sp1479.pass 

sp1481.pass.recode.vcf: 
	vcftools --gzvcf sp1481mutect2.vcf.gz --remove-filtered-all --recode --out sp1481.pass 

sp1483.pass.recode.vcf: 
	vcftools --gzvcf sp1483mutect2.vcf.gz --remove-filtered-all --recode --out sp1483.pass 

sp1485.pass.recode.vcf: 
	vcftools --gzvcf sp1485mutect2.vcf.gz --remove-filtered-all --recode --out sp1485.pass 

sp1487.pass.recode.vcf: 
	vcftools --gzvcf sp1487mutect2.vcf.gz --remove-filtered-all --recode --out sp1487.pass 

sp1452.annovar.vcf: 
	${ANNOVAR}/table_annovar.pl sp1452.pass.recode.vcf ${ANNOVAR}/humandb/ -buildver hg38 -out sp1452anno -remove -protocol knownGene,exac03,gnomad_genome,\
	cosmic70,nci60,clinvar_20170501,avsnp144,dbnsfp33a -operation g,f,f,f,f,f,f,f -nastring . -vcfinput

sp1466.annovar.vcf: 
	${ANNOVAR}/table_annovar.pl sp1466.pass.recode.vcf ${ANNOVAR}/humandb/ -buildver hg38 -out sp1466anno -remove -protocol knownGene,exac03,gnomad_genome,\
	cosmic70,nci60,clinvar_20170501,avsnp144,dbnsfp33a -operation g,f,f,f,f,f,f,f -nastring . -vcfinput


sp1487.pass.vcf:
	gunzip -c sp1476mutect2.vcf.gz |grep 'PASS' >sp1452.pass.vcf
	

AWKOUT:
	#awk '{num=split($8,a,";");val=a[num];sub(/[^0-9]*/,"",a[num]);if(a[num]+0>10){print a[num];print $8 > "out_file"}}'   Input_file
	awk '{num=split($8,a,";");val=a[num];sub(/[^0-9]*/,"",a[num]);if(a[num]+0>10){print a[num];print  > "sp1452.pass.vcf"}}'  sp1452.pass.txt
 

sp1452.coverage.out:
	python coverage-profile.py ${hg38}.fasta ${WHERE}/sp1452/sp1452Aligned.out.sam ${WHERE}/sp1452.coverage.out


sp1452clumpify:
	${BBMAP}/clumpify.sh in=${WHERE}/SP1452_S3_R1_001.fastq.gz in2=${WHERE}/SP1452_S3_R2_001.fastq.gz  out=sp1452_1.fq out2=sp1452_2.fq \
	optical=f 
