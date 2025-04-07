from tqdm import tqdm
from joblib import Parallel, delayed
import numpy as np
from sklearn.decomposition import IncrementalPCA
from sklearn.preprocessing import StandardScaler
import h5py
from scipy.sparse import random
import os
import argparse

# Parse Inputs 
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

parser.add_argument("--M","-M",dest="M",help="GWAS sample size",type=int,default=5000)
parser.add_argument("--L","-L",dest="L",help="Number of independent SNPs",type=int,default=1000)
parser.add_argument("--theta","-t",dest="theta",help="theta parameter",type=float,default=0.05)
parser.add_argument("--fst1","-f1",dest="fst1",help="fst of first split",type=float,default=0.5)
parser.add_argument("--fst2","-f2",dest="fst2",help="fst of second split",type=float,default=0.05)
parser.add_argument("--outfile","-o",dest="outfile",help="outfile",type=str)
parser.add_argument("--genofile","-g",dest="genofile",help="path to chunck directory",type=str)
args=parser.parse_args()

print(args)

# Setup
M = args.M  # Number of individuals
L = args.L  # Number of SNPs
chunk_size_M = 5000  # Number of individuals to process at a time
num_chunks_L = 100
genofile = args.genofile
theta=args.theta
fst1=args.fst1
fst2=args.fst2
outfile = args.outfile
os.makedirs(os.path.dirname(genofile), exist_ok=True)

# Function to simulate genotypes 
def simulate_and_save_h5(M, L, theta, fst1, fst2, num_chunks, output_file):
  
    """Simulates genotype data in chunks and saves to an HDF5 file."""
    
    chunk_size = L // num_chunks  # Number of columns per chunk
    print(chunk_size)
    
    # Open an HDF5 file for writing
    with h5py.File(output_file, "w") as h5f:
        # Create an extendable dataset
        dset = h5f.create_dataset(
            "genotype_matrix", 
            shape=(M, 0),  # Start with zero columns
            maxshape=(M, L),  # Max possible size
            dtype="float32",
            chunks=(M, chunk_size)  # Store in column chunks
        )
        
        for i in range(num_chunks):
            print(f"Processing chunk {i+1}/{num_chunks}...")

            # Define the slice of L to process
            L_start = i * chunk_size
            L_end = (i + 1) * chunk_size if i < num_chunks - 1 else L  # Ensure last chunk covers remainder
            
            # Simulate allele frequencies for this chunk
            #panc = np.random.uniform(0.01, 0.99, L_end - L_start)
            panc = np.full(L_end - L_start, 0.5)
            pA = np.random.normal(panc, np.sqrt(fst1 * panc * (1 - panc)))
            pA = np.clip(pA, 0.01, 0.99)
            pint = np.random.normal(panc, np.sqrt((fst1 - fst2) * panc * (1 - panc)))
            pint = np.clip(pint, 0.01, 0.99)
            pB = np.random.normal(pint, np.sqrt(fst2 * pint * (1 - pint)))
            pC = np.random.normal(pint, np.sqrt(fst2 * pint * (1 - pint)))
            pB = np.clip(pB, 0.01, 0.99)
            pC = np.clip(pC, 0.01, 0.99)

            # Define population sizes
            size_B = round(theta * M * 0.5)
            size_C = round(theta * M * 0.5)
            size_A = M - size_B - size_C

            # Variance components
            sigma2A = (1 / np.log(2 * M)) * (((4 * fst1) - fst2) + 1)
            sigma2BC = (1 / np.log(2 * M)) * ((theta * fst2) + 1)

            # Generate genotype matrices
            Ga = np.random.normal(np.tile(pA, (size_A, 1)), np.sqrt(sigma2A), (size_A, L_end - L_start))
            Gb = np.random.normal(np.tile(pB, (size_B, 1)), np.sqrt(sigma2BC), (size_B, L_end - L_start))
            Gc = np.random.normal(np.tile(pC, (size_C, 1)), np.sqrt(sigma2BC), (size_C, L_end - L_start))

            # Stack and standardize
            G = np.vstack((Ga, Gb, Gc))
            G = (G - np.mean(G, axis=0)) / np.std(G, axis=0)
            
            # Resize dataset to append new chunk
            dset.resize((M, L_end))
            dset[:, L_start:L_end] = G

            print(f"Saved chunk {i+1}/{num_chunks} to {output_file}")

    print("Processing complete!")


# Function to run incremental PCA 
def process_chunk(start, chunk_size, M, L, outfile, ipca, pbar):
  
    end = min(start + chunk_size, M)
    
    # Read from HDF5 (each worker should open the file separately)
    with h5py.File(outfile, "r") as h5f:
        chunk = h5f["genotype_matrix"][start:end, :L]  # Read specific rows

    # Convert sparse matrix if needed
    if isinstance(chunk, np.ndarray):
        chunk_dense = chunk
    else:
        chunk_dense = chunk.toarray() 

    # Fit PCA incrementally
    ipca.partial_fit(chunk_dense)
    
    # Update progress bar
    pbar.update(1)


####### Main ########

# Make Genotype Matrix 
simulate_and_save_h5(M=M, L=L, theta=theta, fst1=fst1, fst2=fst2, num_chunks=num_chunks_L, output_file=genofile)

# Create Incremental PCA object to extract 2 principal components
ipca = IncrementalPCA(n_components=2)

# Do incremental PCA
with tqdm(total=M // chunk_size_M, desc="Processing Chunks") as pbar:
    Parallel(n_jobs=-1)(delayed(process_chunk)(start, chunk_size_M, M, L, genofile, ipca, pbar) 
                         for start in range(0, M, chunk_size_M))
                         
print("Done with PCA fitting")

# Open HDF5 file once
pcs_list = []
with h5py.File(genofile, "r") as h5f, tqdm(total=M // chunk_size_M, desc="Transforming Chunks") as pbar:
    for start in range(0, M, chunk_size_M):
        end = min(start + chunk_size_M, M)

        # Read chunk from HDF5
        chunk = h5f["genotype_matrix"][start:end, :L]  # Read specific rows

        # Ensure it's dense for PCA
        chunk_dense = chunk if isinstance(chunk, np.ndarray) else chunk.toarray()

        # Apply PCA transformation
        pcs_chunk = ipca.transform(chunk_dense)
        pcs_list.append(pcs_chunk)  # Store transformed chunk

        pbar.update(1)  # Update progress bar

# Concatenate all transformed chunks
pcs_all = np.vstack(pcs_list)
print("Final PCA shape:", pcs_all.shape)  

# Save the PCs 
np.savetxt(outfile, pcs_all, delimiter="\t")
