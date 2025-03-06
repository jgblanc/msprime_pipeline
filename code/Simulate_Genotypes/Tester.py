# Import modules 
#import pandas as pd
import numpy as np
from sklearn.decomposition import PCA


import numpy as np
from sklearn.decomposition import IncrementalPCA

import numpy as np
from sklearn.decomposition import PCA

# Define dataset size
n = 500000  # Number of individuals
L = 100000  # Number of SNPs

# Generate genotype matrix (simulated as binomial)
genotype_matrix = np.random.binomial(2, 0.5, size=(n, L)).astype(np.float32)

# Standardize genotype matrix (mean-center & scale)
means = np.mean(genotype_matrix, axis=0)
stds = np.std(genotype_matrix, axis=0, ddof=1)
stds[stds == 0] = 1  # Avoid division by zero
genotype_matrix_std = (genotype_matrix - means) / stds  # Standardization

# Perform PCA with randomized SVD for efficiency
pca = PCA(n_components=2, svd_solver='randomized', random_state=42)
pcs = pca.fit_transform(genotype_matrix_std)

# Save only the first 2 PCs
#np.savetxt("top2_pcs.txt", pcs[:, :2], delimiter="\t")
#print("Saved top 2 PCs to 'top2_pcs.txt'")


# Output first 5 individuals' PCs
print(pcs.shape)
