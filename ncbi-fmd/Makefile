WHERE = /data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/ncbi-fmd
TRIMMOMATIC = /data/WHRI-EndocrinePituitaryGroup/software/Trimmomatic-0.36
HTRANSCRIPTS = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/GRCh37_latest_rna.fna.gz
SALMON = /data/WHRI-EndocrinePituitaryGroup/software/Salmon-0.8.2_linux_x86_64/bin
GTF = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/Homo_sapiens/UCSC/hg19/Annotation/Genes
WHOLEGENOME = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta
BOWTIEINDEX = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex
BOWTIE2INDEX = /data/WHRI-EndocrinePituitaryGroup/sherine/databases/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index
TOPHAT2 = /data/home/hhx779/software/tophat-2.1.1.Linux_x86_64


fm-human-index:
	${SALMON}/salmon index -t ${HTRANSCRIPTS} -i fmhuman-index --type fmd

jf: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/JF-cells_S12_R1_001.fastq.gz -2 ${WHERE}/JF-cells_S12_R2_001.fastq.gz \
	-o JF --seqBias --gcBias --posBias -p 3 --writeUnmappedNames 
tc: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/TC-cells_S11_R1_001.fastq.gz -2 ${WHERE}/TC-cells_S11_R2_001.fastq.gz \
	-o TC --seqBias --gcBias --posBias -p 3 --writeUnmappedNames 

	

sp1452: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1452_S3_R1_001_val_1.fq.gz -2 ${WHERE}/SP1452_S3_R2_001_val_2.fq.gz \
	-o sp1452 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames 

sp1466: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1466_S2_R1_001_val_1.fq.gz -2 ${WHERE}/SP1466_S2_R2_001_val_2.fq.gz \
	-o sp1466 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1470:
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1470_S1_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1470_S1_R2_001_val_2.fq.gz \
	-o sp1470 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames
sp1476: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1476_S4_R1_001_val_1.fq.gz -2 ${WHERE}/SP1476_S4_R2_001_val_2.fq.gz \
	-o sp1476 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames
sp1477: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1477_S5_R1_001_val_1.fq.gz -2 ${WHERE}/SP1477_S5_R2_001_val_2.fq.gz \
	-o sp1477 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames 

sp1479: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1479_S6_R1_001_val_1.fq.gz -2 ${WHERE}/SP1479_S6_R2_001_val_2.fq.gz \
	-o sp1479 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames 

sp1481: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1481_S7_R1_001_val_1.fq.gz -2 ${WHERE}/SP1481_S7_R2_001_val_2.fq.gz \
	-o sp1481 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1483: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1483_S8_R1_001_val_1.fq.gz -2 ${WHERE}/SP1483_S8_R2_001_val_2.fq.gz \
	-o sp1483 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1485: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1485_S9_R1_001_val_1.fq.gz -2 ${WHERE}/SP1485_S9_R2_001_val_2.fq.gz \
        -o sp1485 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames



sp1487: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/SP1487_S10_R1_001_val_1.fq.gz -2 ${WHERE}/SP1487_S10_R2_001_val_2.fq.gz \
	-o sp1487 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames 
	


sp1521: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1521_S21_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1521_S21_R2_001_val_2.fq.gz \
	-o sp1521 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1508:
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1508_S13_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1508_S13_R2_001_val_2.fq.gz \
	-o sp1508 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1510:
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1510_S14_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1510_S14_R2_001_val_2.fq.gz \
	-o sp1510 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1514:
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1514_S15_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1514_S15_R2_001_val_2.fq.gz \
	-o sp1514 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1515: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1515_S16_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1515_S16_R2_001_val_2.fq.gz \
	-o sp1515 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames 

sp1516: 
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1516_S17_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1516_S17_R2_001_val_2.fq.gz \
	-o sp1516 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames


sp1517:
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1517_S18_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1517_S18_R2_001_val_2.fq.gz \
	-o sp1517 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1519:
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1519_S19_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1519_S19_R2_001_val_2.fq.gz \
	-o sp1519 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames


sp1520:
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1520_S20_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1520_S20_R2_001_val_2.fq.gz \
	-o sp1520 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames

sp1522:
	${SALMON}/salmon quant -i fmhuman-index --libType IU -1 ${WHERE}/Sp1522_S22_R1_001_val_1.fq.gz -2 ${WHERE}/Sp1522_S22_R2_001_val_2.fq.gz \
	-o sp1522 --seqBias --gcBias --posBias -p 3 --writeUnmappedNames
