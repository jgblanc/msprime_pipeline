import numpy as np
from sklearn.decomposition import IncrementalPCA
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
import h5py
from scipy.sparse import random
import os

# Define dataset size
n = 5000  # Number of individuals
L = 1000  # Number of SNPs
chunk_size = 500  # Number of individuals to process at a time
chunk_dir = 'chunks_data/'
os.makedirs(chunk_dir, exist_ok=True)

# Create Incremental PCA object to extract 2 principal components
ipca = IncrementalPCA(n_components=2)

# Create a scaler for standardization
scaler = StandardScaler()

# Process chunks in parallel
def process_chunk(start, chunk_size, n, L, chunk_dir, scaler, ipca):
    end = min(start + chunk_size, n)
    chunk = random(end - start, L, density=0.1, format='csr')
    chunk.data = np.random.binomial(2, 0.5, size=chunk.data.shape).astype(np.float32)

    # Standardize
    chunk_dense = scaler.fit_transform(chunk.toarray())

    # Save as HDF5
    chunk_filename = f'{chunk_dir}/chunk_{start}_{end}.h5'
    with h5py.File(chunk_filename, 'w') as hf:
        hf.create_dataset('chunk', data=chunk_dense)

    # Fit PCA incrementally
    ipca.partial_fit(chunk_dense)


# Parallelize chunk processing
Parallel(n_jobs=-1)(delayed(process_chunk)(start, chunk_size, n, L, chunk_dir, scaler, ipca) for start in range(0, n, chunk_size))
print("Done with fitting")

# Transform and concatenate results
pcs_list = []
for start in range(0, n, chunk_size):
    end = min(start + chunk_size, n)
    chunk_filename = f'{chunk_dir}/chunk_{start}_{end}.h5'

    with h5py.File(chunk_filename, 'r') as hf:
        chunk_dense = hf['chunk'][:]
        pcs_chunk = ipca.transform(chunk_dense)
        pcs_list.append(pcs_chunk)

pcs_all = np.vstack(pcs_list)
print(pcs_all.shape)








# Process genotype matrix in chunks to avoid memory overload
#for start in range(0, n, chunk_size):
#
#    print(start)
#    
#    end = min(start + chunk_size, n)  # Ensure the last chunk is not out of bounds
#    chunk = np.random.binomial(2, 0.5, size=(end - start, L)).astype(np.float32)
#
#    # Standardize the chunk
#    chunk = scaler.fit_transform(chunk)
#
#    # Save the chunk as an NPZ file
#    chunk_filename = f'{chunk_dir}/chunk_{start}_{end}.npz'
#    np.savez_compressed(chunk_filename, chunk=chunk)
#
#    # Fit Incremental PCA on the chunk incrementally
#    ipca.partial_fit(chunk)

#print("Done with fitting")
    
# Load in chunks and transform    
#pcs_list = []    
#for start in range(0, n, chunk_size):

#    print(start)
    
#    end = min(start + chunk_size, n)  # Ensure the last chunk is not out of bounds
#    chunk_filename = f'{chunk_dir}/chunk_{start}_{end}.npz'

    # Load chunk
#    loaded_chunk = np.load(chunk_filename)['chunk']

    # Transform and save PCs
#    pcs_chunk = ipca.transform(loaded_chunk)
#    pcs_list.append(pcs_chunk)

#print("done with PCA")
    
# Concatenate results for all individuals
#pcs_all = np.vstack(pcs_list)
#print(pcs_all.shape)
