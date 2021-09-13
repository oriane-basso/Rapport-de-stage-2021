#Analyse de structure
library(LEA)
output=vcf2geno("Chr_01.vcf","CHR01genoformat1708.geno")

obj.at = snmf("CHR01genoformat1708.geno", K = 1:14, entropy = T, repetitions=10,
              CPU = 12, project = "new")

pdf(file = "plotCHR1.pdf")
plot(obj.at, col = "blue4", cex = 1.4, pch = 19)
dev.off()
export.snmfProject("CHR01genoformat1708.snmfProject")