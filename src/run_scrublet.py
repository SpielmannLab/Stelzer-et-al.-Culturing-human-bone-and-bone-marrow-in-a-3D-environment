# --- Load all the libraries 
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from numpy import savetxt

# --- Read all the arguments
print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))

# --- Parameters
sc_matrix_file = sys.argv[1]    # Name of count matrix (*.mtx) file generated after filtering
npcs = int(sys.argv[2])         # Number of principal components to use
exp_db_rate = float(sys.argv[3])# Expected rate of doublet (as a fraction)

print("sc_matrix_file: ", sc_matrix_file)
print("npcs: ", npcs)
print("exp_db_rate: ", exp_db_rate)

# --- Run
sc_matrix = scipy.io.mmread(sc_matrix_file).T.tocsc()
print('Counts matrix has the dimensions: {} rows, {} columns'.format(sc_matrix.shape[0], sc_matrix.shape[1]))

print("Now running scrublet:")
scrub = scr.Scrublet(sc_matrix, expected_doublet_rate=exp_db_rate)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=npcs)

#Run UMAP to show the scrublet scores in the UMAP
print('Running UMAP...')
umap = scr.get_umap(scrub.manifold_obs_, n_neighbors=10, min_dist=0.1)

# --- Save the output csv files
# Save duplet score for cells
filename=sc_matrix_file.replace(".mtx","_doublet_score.csv")
savetxt(filename, scrub.doublet_scores_obs_, delimiter=',')

# Save simulated duplets
filename=sc_matrix_file.replace(".mtx","_doublet_score_sim.csv")
savetxt(filename, scrub.doublet_scores_sim_, delimiter=',')

# run and save UMAP
filename=sc_matrix_file.replace(".mtx","_doublet_umap.csv")
savetxt(filename, umap, delimiter=',')
    
print("Finished running scrublet. Saved CSV files.")
