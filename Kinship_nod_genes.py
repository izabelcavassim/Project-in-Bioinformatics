# Making a matrix of the interesting genes (nodA, nodD (both plasmid genes),  recA and rpoB are chromosomal genes)

# First take just the genes that all of the strains have
# Taking the file with the Nod genes and intersecting with the snp data

# Producing matrix of similarity between strains for posterior PCA
import numpy as np
import pylab as pl
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D

def normalize_columns(matrix):
    """Normalizes a matrix' columns to mean 0 and std 1."""
    mean = np.mean(matrix, axis=0)
    std = np.std(matrix, axis=0)
    matrix = (matrix - mean) / std
    return np.nan_to_num(matrix)

def coordinates(subset, sett):
	"""Indexes of the strains that contains that specific gene."""
	coord = list()
	boolean = list()
	for i in subset: #the list of strain names for the specific gene
		coord.append(int(i[:4]) - int(sett[0][:4])) #sett the entire list of strain_names
	return coord

def kinship(maf, genes, ind, name, strain_names_all, strain_names_sub, snp_name):
	# Filtering by MAF
	genes_maf = list()
	for gene in genes:
		num_minor_per_column = np.sum(gene == 0, axis = 0)/float(gene.shape[0]) # zero is minor allele
		num_major_per_column = np.sum(gene == 1, axis = 0) # one is the major allele
		index_MAF = num_minor_per_column >= maf # MAF = 0.3
		gene = gene[:,index_MAF]
		genes_maf.append(gene)

	K = np.zeros((ind, ind)) # Kinship matrix initialized 
	CM = np.zeros((ind, ind)) # Count Markers initialized

	for G in genes_maf:
		NG = normalize_columns(G)
		# Find the indexes of individuals that match with big matrix 
		ind = coordinates(strain_names_sub, strain_names_all)
		rows = [[x] for x in ind]
		#K[np.ix_(ind, ind)]
		K += np.dot(NG, NG.T)
		CM += NG.shape[1]
		#K[rows,ind] += np.dot(NG, NG.T)
		#CM[rows, ind] += NG.shape[1] # number of snps
	print CM
	K = K/CM	#element wise division
	# Saving the matrix
	np.savetxt(name, K, delimiter=',')
	
	# Matrix of strains X Genes (their respective snps)
	Gene_matrix = np.concatenate((genes_maf), axis=1)
	Gene_matrix = normalize_columns(Gene_matrix)
	np.savetxt(snp_name, Gene_matrix, fmt='%.18e', delimiter=',')
	# Save as numpy object 
	#np.save('Gene_matrix_snps.npy', Gene_matrix, allow_pickle=True, fix_imports=True)
	return K

def heat_map(K, show, name):
	pl.pcolor(K)
	pl.colorbar()
	pl.title(name)
	pl.savefig(name)
	if show == True:
		pl.show()

# Importing the strain names that will be important for locate the positions in the kinship matrix 
strain_names_all = open('all_strains.txt', 'r').read().split('\n')

# NodA genes
snps_NodA = [np.load('nodA.msa.prank.fna.snps.npy')]
strains_NodA = open('nodA.strains.txt', 'r').read().split('\n')

K_NodA = kinship(0.01, snps_NodA, len(strains_NodA), 'Kinship_matrix_NodA.csv', strain_names_all, strains_NodA, 'Relationship_matrix_nodA.csv')
heat_map(K_NodA, True, 'NodA_heatmap.png')


# NodD genes
snps_NodD = [np.load('nodD.msa.prank.fna.snps.npy')]
strains_NodD = open('nodD.strains.txt', 'r').read().split('\n')
K_NodD = kinship(0.01, snps_NodD, len(strains_NodD), 'Kinship_matrix_NodD.csv', strain_names_all, strains_NodD, 'Relationship_matrix_nodD.csv')
heat_map(K_NodD, True, 'NodD_heatmap.png')

# RecA genes
snps_recA = [np.load('recA.msa.prank.fna.snps.npy')]
strains_recA = open('recA.strains.txt', 'r').read().split('\n')
K_recA = kinship(0.01, snps_recA, len(strains_recA), 'Kinship_matrix_recA.csv', strain_names_all, strains_recA, 'Relationship_matrix_recA.csv')
heat_map(K_recA, True, 'recA_heatmap.png')


# RpoB genes
snps_rpoB = [np.load('rpoB.msa.prank.fna.snps.npy')]
strains_rpoB = open('rpoB.strains.txt', 'r').read().split('\n')
K_rpoB = kinship(0.01, snps_rpoB, len(strains_rpoB), 'Kinship_matrix_rpoB.csv', strain_names_all, strains_rpoB, 'Relationship_matrix_rpoB.csv')
heat_map(K_rpoB, True, 'rpoB_heatmap.png')


#########################
# Just using the strains that have both NodA NodD 
#########################

in_common = list(set(strains_NodA).intersection(strains_NodD))  
print len(in_common)