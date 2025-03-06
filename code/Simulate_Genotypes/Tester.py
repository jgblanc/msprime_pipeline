import numpy as np
from sklearn.decomposition import IncrementalPCA

# Define dataset size
n = 500000  # Number of individuals
L = 100000  # Number of SNPs
chunk_size = 10000  # Number of SNPs to process at a time

# Create Incremental PCA object to extract 2 principal components
ipca = IncrementalPCA(n_components=2)

# Process genotype matrix in chunks to avoid memory overload
for _ in range(L // chunk_size):
    # Simulate a chunk of genotype data (n x chunk_size)
    chunk = np.random.binomial(2, 0.5, size=(n, chunk_size)).astype(np.float32)
    
    # Standardize the chunk (mean-center & scale)
    means = np.mean(chunk, axis=0)
    stds = np.std(chunk, axis=0, ddof=1)
    stds[stds == 0] = 1  # Avoid division by zero
    chunk = (chunk - means) / stds  # Standardization

    # Incremental PCA update with the chunk
    ipca.partial_fit(chunk)

# After processing all chunks, transform the data and compute the first 2 PCs
pcs = ipca.transform(chunk)  # Transform the final chunk

# Save the first 2 PCs
#np.savetxt("top2_pcs.txt", pcs[:, :2], delimiter="\t")
#print("Saved top 2 PCs to 'top2_pcs.txt'")


# Output first 5 individuals' PCs
print(pcs.shape)
