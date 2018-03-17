library(tximport)
library(rjson)
library(readr)

samples <- read.table(file.path("/data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/salmon", "samples.txt"), header = FALSE)

files <- file.path("/data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/", "salmon", samples$V1, "quant.sf")

names(files) <- paste0("sample", 1:20)

tx2gene <- read.csv(file.path("/data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/salmon", "tx2gene.csv"))

txi.salmon1 <- tximport(files[1], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon1$counts, "sp1508quant.counts")


txi.salmon2 <- tximport(files[2], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon2$counts, "sp1510quant.counts")

txi.salmon3 <- tximport(files[3], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon3$counts, "sp1514quant.counts")


txi.salmon4 <- tximport(files[4], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon4$counts, "sp1515quant.counts")

txi.salmon5 <- tximport(files[5], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon5$counts, "sp1516quant.counts")

txi.salmon6 <- tximport(files[6], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon6$counts, "sp1517quant.counts")

txi.salmon7 <- tximport(files[7], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon7$counts, "sp1519quant.counts")

txi.salmon8 <- tximport(files[8], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon8$counts, "sp1520quant.counts")

txi.salmon9 <- tximport(files[9], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon9$counts, "sp1521quant.counts")

txi.salmon10 <- tximport(files[10], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon10$counts, "sp1522quant.counts")

txi.salmon11 <- tximport(files[11], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon11$counts, "sp1470quant.counts")

txi.salmon12 <- tximport(files[12], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon12$counts, "sp1466quant.counts")

txi.salmon13 <- tximport(files[13], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon13$counts, "sp1452quant.counts")

txi.salmon14 <- tximport(files[14], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon14$counts, "sp1476quant.counts")

txi.salmon15 <- tximport(files[15], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon15$counts, "sp1477quant.counts")

txi.salmon16 <- tximport(files[16], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon16$counts, "sp1479quant.counts")

txi.salmon17 <- tximport(files[17], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon17$counts, "sp1481quant.counts")

txi.salmon18 <- tximport(files[18], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon18$counts, "sp1483quant.counts")

txi.salmon19 <- tximport(files[19], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon19$counts, "sp1485quant.counts")

txi.salmon20 <- tximport(files[20], type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.salmon20$counts, "sp1487quant.counts")
