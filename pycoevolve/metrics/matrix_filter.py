#/usr/bin/env python
import numpy as np

def bivar_gauss(coev_matrix,epsilon=1e-8):
    resi_num = coev_matrix.shape[0]
    resi_mean = np.tile(coev_matrix.mean(axis=0),[resi_num,1])
    resi_var = np.tile(coev_matrix.var(axis=0),[resi_num,1])
    background_mean = ((resi_mean * resi_var.T) \
            + (resi_mean.T * resi_var)) / (resi_var + resi_var.T)
    background_var = (resi_var*resi_var.T)/(resi_var + resi_var.T)
    return (coev_matrix - background_mean)/np.sqrt(background_var + epsilon)
