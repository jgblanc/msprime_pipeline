import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.decomposition import IncrementalPCA
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
import h5py
from scipy.sparse import random
import os


# Define the parameter values
theta_vals = np.exp(np.linspace(np.log(0.01), np.log(0.8), num=8))
fst1_vals = [0.6]
fst2_vals = [0.25, 0.5, 0.125]
reps = range(1, 21)  # 20 replicates
M = 500  # Population size
L = 1000  # Number of SNPs

# Generate all combinations of parameters and replicates using pandas
param_grid = pd.DataFrame(
    [(rep, theta, fst1, fst2) for rep in reps for theta in theta_vals for fst1 in fst1_vals for fst2 in fst2_vals],
    columns=["rep", "theta", "fst1", "fst2"]
)

# Function to apply simulate_normals from before
def simulate_normals(L, M, theta, fst1, fst2):
    sigma2A = (1 / np.log(2 * M)) * ((4 * fst1 - fst2) + 1)
    sigma2BC = (1 / np.log(2 * M)) * ((theta * fst2) + 1)
    panc = np.random.uniform(0.01, 0.99, L)
    pA = np.random.normal(panc, np.sqrt(fst1 * panc * (1 - panc)))
    pA = np.clip(pA, 0.01, 0.99)
    pint = np.random.normal(panc, np.sqrt((fst1 - fst2) * panc * (1 - panc)))
    pint = np.clip(pint, 0.01, 0.99)
    pB = np.random.normal(pint, np.sqrt(fst2 * pint * (1 - pint)))
    pC = np.random.normal(pint, np.sqrt(fst2 * pint * (1 - pint)))
    pB = np.clip(pB, 0.01, 0.99)
    pC = np.clip(pC, 0.01, 0.99)

    size_A = round((1 - theta) * M)
    size_B = round(theta * M * 0.5)
    size_C = round(theta * M * 0.5)

    Ga = np.random.normal(np.tile(pA, (size_A, 1)), np.sqrt(sigma2A), (size_A, L))
    Gb = np.random.normal(np.tile(pB, (size_B, 1)), np.sqrt(sigma2BC), (size_B, L))
    Gc = np.random.normal(np.tile(pC, (size_C, 1)), np.sqrt(sigma2BC), (size_C, L))

    G = np.vstack((Ga, Gb, Gc))
    G = (G - np.mean(G, axis=0)) / np.std(G, axis=0)

    ipca = IncrementalPCA(n_components=2)
    chunk_size = 10_000
    for i in range(L // chunk_size + 1):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, L)
        G_chunk = G[:, start_idx:end_idx]
        ipca.partial_fit(G_chunk)

    pcs = ipca.transform(G)
    sample_e1 = pcs[:, 0]
    sample_e2 = pcs[:, 1]

    pop_e1 = np.concatenate([np.repeat(-theta, size_A), np.repeat(1 - theta, size_B + size_C)]) / \
             (theta * np.sqrt(((theta ** 2 + 1) / theta) - 1))
    pop_e2 = np.concatenate([np.repeat(0, size_A), np.repeat(1, size_B), np.repeat(-1, size_C)]) / np.sqrt(theta)

    b2_e1 = np.corrcoef(sample_e1, pop_e1)[0, 1] ** 2
    b2_e2 = np.corrcoef(sample_e2, pop_e2)[0, 1] ** 2

    return {"b2_e1": b2_e1, "b2_e2": b2_e2}

# Run the simulations and store the results in the DataFrame
results = []
for _, row in tqdm(param_grid.iterrows(), total=param_grid.shape[0]):
    result = simulate_normals(L, M, row['theta'], row['fst1'], row['fst2'])
    results.append(result)

# Convert results to DataFrame
results_df = pd.DataFrame(results)
param_grid = param_grid.reset_index(drop=True)
final_results = pd.concat([param_grid, results_df], axis=1)

# Save the results to a text file
final_results.to_csv('simulation_results.txt', sep='\t', index=False)

# Print the results
print(final_results)



