# Producing matrix of similarity between strains for posterior PCA
import numpy as np
import pylab as pl
import glob
from matplotlib.mlab import PCA 

#This is the gene SNPs matrix
genes = [np.load(f) for f in glob.glob("*.snps.npy")] 

def normalize_columns(matrix):
    """Normalizes a matrix' columns to mean 0 and std 1."""
    mean = np.mean(matrix, axis=0)
    std = np.std(matrix, axis=0)
    matrix = (matrix - mean) / std
    return np.nan_to_num(matrix)

def kinship(maf, genes, ind, name):
	# Filtering by MAF
	genes_maf = list()
	for gene in genes:
		num_minor_per_column = np.sum(gene == 0, axis = 0)/float(gene.shape[0]) # zero is minor allele divided by the number of rows
		index_MAF = num_minor_per_column >= maf # MAF = 0.3
		gene = gene[:,index_MAF]
		genes_maf.append(gene)

	K = np.zeros((ind,ind)) # Kinship matrix initialized 
	CM = np.zeros((ind,ind)) # Count Markers initialized

	for G in genes_maf:
		NG = normalize_columns(G)
		K += np.dot(NG, NG.T)
		CM += NG.shape[1]

	K = K/CM
	# Saving the matrix
	np.savetxt(name, K, delimiter=',')
	
	# Matrix of individuals X Genes (their respective snps)
	#print genes_maf
	Gene_matrix_big = np.concatenate((genes_maf[1:50]), axis=1)
	print type(Gene_matrix_big)
	print Gene_matrix_big.shape
	#print Gene_matrix_big[:,Gene_matrix_big.shape[1]-1]
	#Gene_matrix = normalize_columns(Gene_matrix)

	#Duplicating the rows to be able to use in 
	print Gene_matrix_big.shape
	seq = range(0,44)
	repeats = np.repeat(seq, 2) 
	Gene_matrix_big = Gene_matrix_big[repeats,:]

	print 'the modified matrix has shape', Gene_matrix_big.shape
	np.savetxt('Gene_matrix_big.txt', Gene_matrix_big, delimiter=',', fmt='%.0f') # saving as txt
	np.save('Gene_matrix_snps.npy', Gene_matrix_big, allow_pickle=True, fix_imports=True) #saving as npy object
	return K

K = kinship(0.3, genes, 44, 'Kinship_matrix_all.csv')

# The mean of the diagonal should be close to 1
print 'The mean of the diagonal is:'
print np.mean(np.diag(K))

# The elements above the diagonal should be close to zero if they are not related  
print '\nThe sum of the elements above diagonal is:'
print np.sum(np.tril(K)/(43*42/2))

# Searching for extreme individuals:
print np.where(np.diag(K) > 2)

# A heatmap of the matrix
pl.pcolor(K)
pl.colorbar()
pl.title('All core genes of 44 strains')
pl.savefig('heat_map_allcoregenes.png')
pl.show()

# Using a package, it may scale
pca = PCA(K.T)

# Making the figure
pl.plot(pca.Y[0:9,0],pca.Y[0:9,1], 'o', markersize=7, color='blue', alpha=0.5, label='UK')
pl.plot(pca.Y[9:12,0], pca.Y[9:12,1], '^', markersize=7, color='red', alpha=0.5, label='FR')
pl.plot(pca.Y[12:43,0],pca.Y[12:43,1], 'o', markersize=7, color='green', alpha=0.5, label='DK')

pl.xlabel('PC1')
pl.ylabel('PC2')
#pl.xlim([-4,4])
#pl.ylim([-4,4])
pl.legend()
pl.title('Population structure based on coregenes')
pl.show()