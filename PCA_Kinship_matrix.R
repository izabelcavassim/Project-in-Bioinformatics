#Thomas said that we couldn't use the kinship matrix as an entry in the data, I just tested that we can use
#http://www.bioinformaticstutorials.com/?p=124
#https://tgmstat.wordpress.com/2013/11/28/computing-and-visualizing-pca-in-r/
#https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

setwd('C:/Users/MariaIzabel/Desktop/MASTER/PROJECT IN BIOINFORMATICS/prank_alignments_snps/')
require(xlsx)
strain_names = read.table('strain_ids_hand.txt',sep = '', row.names = NULL, header= FALSE )
origin = read.xlsx('Rhizobium_soiltypes.xlsx', sheetIndex = 1)
strain_names[,4] = paste(strain_names[,1], strain_names[,2])

colnames(strain_names) = c('Strain', 'Origin', 'Country', 'Mix')

install.packages('RcppCNPy')
library(RcppCNPy)

#Transforming the Kinship matrix in eigen vectors manually and then plot:
#Eigen vectors
kinship_matrix = read.table('kinship_matrix.csv', head = F, sep = ',')
row.names(kinship_matrix) = strain_names$Mix
eigenVecs <- eigen(kinship_matrix)$vectors
lines(x=c(0, eigenVecs[1,1]), y=c(0, eigenVecs[2,1]), col="red")
lines(x=c(0, eigenVecs[1,2]), y=c(0, eigenVecs[2,2]), col="red")

#Eigen values
eVal1 <- eigen(kinship_matrix)$values[1]
eVal2 <- eigen(kinship_matrix)$values[2]
lines(x=c(0, eVal1*eigenVecs[1,1]), y=c(0, eVal1*eigenVecs[2,1]), col="red")
lines(x=c(0, eVal2*eigenVecs[1,2]), y=c(0, eVal2*eigenVecs[2,2]), col="red")

rowFeatVec <- t(eigenVecs)
rowDataAdj <- t(kinship_matrix)
transFData <- rowFeatVec %*% rowDataAdj
finalData <- rowFeatVec %*% rowDataAdj
plot(t(finalData),  pch = c(16, 17, 18)[as.numeric(strain_names$Country)],  col = c("red", "green","blue")[as.numeric(strain_names$Country)])
barplot(eigen(kinship_matrix)$values)

require(ggfortify)
#Using the package prcomp and extract automatically the eigen values and plot:

pca_prco <- prcomp(kinship_matrix)
plot(pca_proc$x, pch=seq(9)) #the pca plot
plot(pca_proc) #the proportion of variance capture by each PC
autoplot(pca_proc, data = strain_names, colour = 'Country')