import numpy as np
import pylab as pl
import glob
from matplotlib.mlab import PCA 
from functools import reduce
from matplotlib_venn import venn3

# This is the gene SNPs matrix
genes = [np.load(f) for f in glob.glob("*.snps.npy")] 

print len(genes)

# Simulating missing data: in random rows input nan: it will be always gene.shape[1]/10 random strains that will not have the gene 
for gene in genes:
	gene[np.random.randint(44, size=gene.shape[1]/5),:] = np.nan

Matrix_counts = np.zeros((44,len(genes)+1))
print Matrix_counts.shape

col = 0 
for gene in genes:
	for i in range(len(gene)): 
		if np.isnan(gene[i,1]) == False: #looking just at the first collumn/snp
			Matrix_counts[i,col] = 1
	col += 1

pl.pcolor(Matrix_counts)
pl.colorbar()
pl.show()

# Using a package, it may scale
pca = PCA(Matrix_counts.T)

# Making the figure
pl.plot(pca.Y[0:8,0],pca.Y[0:8,1], 'o', markersize=7, color='blue', alpha=0.5, label='UK')
pl.plot(pca.Y[9:12,0], pca.Y[9:12,1], '^', markersize=7, color='red', alpha=0.5, label='FR')
pl.plot(pca.Y[12:43,0],pca.Y[12:43,1], 'o', markersize=7, color='green', alpha=0.5, label='DK')

pl.xlabel('x_values')
pl.ylabel('y_values')
#pl.xlim([-4,4])
#pl.ylim([-4,4])
pl.legend()
pl.title('Population structure based on presence-abscence of genes')
pl.show()

#####################################################
# Counting genes that are exclusive for certain groups
#####################################################

index_UK = np.where(Matrix_counts[0:8,:] >= 1)
index_FR = np.where(Matrix_counts[9:11,:] >= 1)
index_DK = np.where(Matrix_counts[12:43,:] >= 1)

print 'intersection of all the groups'
set4 = reduce(np.intersect1d, (index_FR,index_UK, index_DK)) # the intersection of all
set1 = np.intersect1d(index_FR,index_UK)
set1 = np.setdiff1d(set1, set4) # exclusive of set1 minus the intersection of the 3 together
set2 = np.intersect1d(index_FR,index_DK)
set2 = np.setdiff1d(set2, set4) # exclusive of set2 minus the intersection of the 3 together
set3 = np.intersect1d(index_UK, index_DK)
set3 = np.setdiff1d(set3, set4) # exclusive of set3 minus the intersection of the 3 together
set5 = np.setdiff1d(index_UK, np.append(set1, set3)) # exclusive elements of UK
set5 = np.setdiff1d(set5, set4)
set6 = np.setdiff1d(index_DK, np.append(set2, set3)) # exclusive elements of DK
set6 = np.setdiff1d(set6, set4)
set7 = np.setdiff1d(index_FR, np.append(set1, set2)) # exclusive elements of FR
set7 = np.setdiff1d(set7, set4)

print len(set1)
print len(set2)
print len(set3)
print 'intersection of all:' 
print len(set4)
print 'exclusive elements of UK:'
print len(set5)
print 'exclusive elements of DK:'
print len(set6)
print 'exclusive elements of FR:'
print len(set7)
print 'total sum:'
print sum([len(set1),len(set2),len(set3),len(set4),len(set5),len(set6),len(set7)])
# easy_install matplotlib-venn (install the module)

####################
# Making the plot Venn Diagram
######################
from matplotlib_venn import venn3, venn3_circles


# Subset sizes
s = (
    len(set5),    # UK
    len(set6),    # DK
    len(set3),    # UK/DK
    len(set7),    # FR
    len(set2),    # DK/FR
    len(set4),  # DK/FR/UK
    len(set1),    # UK/FR
)
print s