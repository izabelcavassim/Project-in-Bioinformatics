# PCA analysis of the Relationship Matrix between strains
install.packages('xlsx')
require(xlsx)
setwd('C:/Users/MariaIzabel/Desktop/MASTER/PROJECT IN BIOINFORMATICS/prank_alignments_snps/')

## For all genes not specifying anything
#Dowloading the kinship matrix
relationship_matrix = read.table(file="Kinship_matrix_all.csv",header=FALSE, sep =',' )
#Rename the matrix:
strain_names = read.table('strain_ids_hand.txt',sep = '', row.names = NULL, header= FALSE )
origin = read.xlsx('Rhizobium_soiltypes.xlsx', sheetIndex = 1)
strain_names[,4] = paste(strain_names[,1], strain_names[,2])

colnames(strain_names) = c('Strain', 'Origin', 'Country', 'Mix')
row.names(relationship_matrix) = strain_names$Mix

#Eigen vectors:
#Unlike princomp, variances are computed with the usual divisor N - 1
pca = prcomp(relationship_matrix, center =TRUE, scale=TRUE)
summary(pca)
ls(pca)
#PC are the sum of the total of the variance contained in the dataset
pca$sdev^2
screeplot(pca, type ='line')
biplot(pca, cex=c(1,0.3))

require(ggfortify)
autoplot(pca, data = strain_names, colour = 'Country', label =TRUE, label.size = 4, shape =TRUE)
autoplot(pca, data = strain_names, colour = 'Country', main = 'Population structure of strains', frame = T)


########################
# Plasmid specificity 
########################
#Used the normalized markers instead of kinship
relationship_matrix_plasmid = read.table(file="Gene_matrix_plasmid.csv",header=FALSE, sep =',' )
kinship_matrix_plasmid = read.table('Kinship_matrix_plasmid.csv', header = F, sep = ',')

pca = prcomp(relationship_matrix_plasmid)

row.names(relationship_matrix_plasmid) = strain_names$Mix
colnames(relationship_matrix_plasmid) = strain_names$Mix

#Eigen vectors:
#Unlike princomp, variances are computed with the usual divisor N - 1
pca_plasm = prcomp(relationship_matrix_plasmid, center =TRUE, scale=TRUE)
#PC are the sum of the total of the variance contained in the dataset
pca_plasm$sdev^2
screeplot(pca_plasm, type ='line')
biplot(pca_plasm, cex=c(1,0.3))

autoplot(pca_plasm, data = strain_names, colour = 'Country', label =TRUE, label.size = 4, shape =FALSE)
autoplot(pca_plasm, data = strain_names, colour = 'Country')

## Nod genes:
kinship_matrix_nod = read.table('Kinship_matrix_nod.csv', header = F, sep = ',')
row.names(kinship_matrix_nod) = strain_names$Mix
pca_nod = prcomp(kinship_matrix_nod, center =TRUE, scale=TRUE)
#PC are the sum of the total of the variance contained in the dataset
pca_nod$sdev^2
screeplot(pca_nod, type ='line')
biplot(pca_plasm, cex=c(1,0.3))

autoplot(pca_nod, data = strain_names, colour = 'Country', label =TRUE, label.size = 4, shape =FALSE)
autoplot(pca_nod, data = strain_names, colour = 'Country')

