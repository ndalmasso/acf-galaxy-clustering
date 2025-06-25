"""
Two-Point Correlation Function (TPCF) Analysis Library

Functions for calculating spatial correlation functions using bootstrap 
and jackknife resampling methods.

Parameter specifications:
- data: tuple of (x_coords, y_coords) in pixels
- rand: tuple of (x_coords, y_coords) in pixels  
- x: bin edges in arcsec
- px_scale: conversion factor from pixels to arcsec
- size: jackknife circular region radius in arcsec
- perc_data: fraction of data points to sample (0-1)
- perc_rand: fraction of random points to sample (0-1)
- hist: if True use np.histogram, if False use custom pair counting
"""

import warnings
from typing import Tuple
import numpy as np
import scipy.spatial.distance as distance
from tqdm import tqdm

warnings.simplefilter("ignore")


def bootstrap(data, rand, x, nboot, perc_data, perc_rand, px_scale, hist=False):
    """Bootstrap resampling for TPCF estimation."""
    progress_bar = tqdm(total=nboot, desc="Bootstrapping", unit="iteration")
    
    max_perc = len(rand[0]) / len(data[0])
    n_sample = int(perc_data * len(data[0]))
    r_sample = int(min(max_perc, perc_rand) * n_sample)
    
    results = {'est': [], 'centre': [], 'rr': [], 'dr': [], 'dd': []}
    valid_samples = 0

    while valid_samples < nboot:
        try:
            # Resample
            data_idx = np.random.choice(len(data[0]), n_sample, replace=False)
            rand_idx = np.random.choice(len(rand[0]), r_sample, replace=False)
            
            re_data = ([data[0][i] for i in data_idx], [data[1][i] for i in data_idx])
            re_rand = ([rand[0][i] for i in rand_idx], [rand[1][i] for i in rand_idx])
            
            # Calculate TPCF
            centre, est, dd, dr, rr = TPCF(re_data, re_rand, x, px_scale, hist)
            
            if not (np.any(np.isnan(est)) or np.any(np.isinf(est))):
                for key, val in zip(['centre', 'est', 'dd', 'dr', 'rr'], 
                                   [centre, est, dd, dr, rr]):
                    results[key].append(val)
                valid_samples += 1
                progress_bar.update(1)
        except:
            continue
    
    progress_bar.close()
    
    cov_mat = np.cov(np.array(results['est']).T)
    
    return (results['est'], 
            cov_mat,
            np.mean(results['est'], axis=0),
            np.mean(results['centre'], axis=0),
            np.mean(results['dd'], axis=0),
            np.mean(results['dr'], axis=0),
            np.mean(results['rr'], axis=0))


def jackknife(data, rand, x, njack, size, px_scale, hist=False):
    """Jackknife resampling using circular regions."""
    progress_bar = tqdm(total=njack, desc="Jackknife", unit="iteration")
    
    data, rand = np.asarray(data), np.asarray(rand)
    results = {'est': [], 'centre': [], 'rr': [], 'dr': [], 'dd': []}
    valid_samples = 0

    while valid_samples < njack:
        try:
            # Random circular region
            center_idx = np.random.randint(len(rand[0]))
            center_x, center_y = rand[0][center_idx], rand[1][center_idx]
            
            data_dist = np.sqrt((data[0] - center_x)**2 + (data[1] - center_y)**2)
            rand_dist = np.sqrt((rand[0] - center_x)**2 + (rand[1] - center_y)**2)
            
            sample_data = tuple(d[data_dist <= 2*(size/px_scale)] for d in data)
            sample_rand = tuple(r[rand_dist <= 2*(size/px_scale)] for r in rand)
            
            if len(sample_data[0]) < 2 or len(sample_rand[0]) < 2:
                continue
                
            centre, est, dd, dr, rr = TPCF(sample_data, sample_rand, x, px_scale, hist)
            
            if not (np.all(np.isnan(est)) or np.all(np.isinf(est))):
                for key, val in zip(['centre', 'est', 'dd', 'dr', 'rr'], 
                                   [centre, est, dd, dr, rr]):
                    results[key].append(val)
                valid_samples += 1
                progress_bar.update(1)
        except:
            continue
    
    progress_bar.close()
    
    cov_mat = np.cov(np.array(results['est']).T)
    
    return (results['est'],
            cov_mat,
            np.mean(results['est'], axis=0),
            np.mean(results['centre'], axis=0),
            np.mean(results['dd'], axis=0),
            np.mean(results['dr'], axis=0),
            np.mean(results['rr'], axis=0))


def TPCF(data, randp, x, px_scale, hist=False):
    """Calculate Two-Point Correlation Function using Landy-Szalay estimator."""
    r_points, d_points = np.array(randp).T, np.array(data).T
    Nr, Nd = r_points.shape[0], d_points.shape[0]
    factor = Nr / Nd

    # Calculate pairwise distances
    d_rr = distance.cdist(r_points, r_points, 'euclidean') * px_scale
    d_dd = distance.cdist(d_points, d_points, 'euclidean') * px_scale
    d_dr = distance.cdist(d_points, r_points, 'euclidean') * px_scale
    
    # Count pairs in bins
    if hist:
        dd_p, _ = np.histogram(d_dd, bins=x)
        dr_p, _ = np.histogram(d_dr, bins=x)
        rr_p, _ = np.histogram(d_rr, bins=x)
    else:
        dd_p, _ = pairs_count(d_dd.flatten(), x)
        dr_p, _ = pairs_count(d_dr.flatten(), x)
        rr_p, _ = pairs_count(d_rr.flatten(), x)

    # Landy-Szalay estimator
    rr_p = np.where(rr_p == 0, 1, rr_p)  # Avoid division by zero
    est = (dd_p * factor**2 - 2 * factor * dr_p + rr_p) / rr_p
    centre = (x[:-1] + x[1:]) / 2
    
    return centre, est, dd_p, dr_p, rr_p


def pairs_count(distances, x):
    """Count pairs within distance bins."""
    nbin = len(x) - 1
    pairs = np.zeros(nbin, dtype=np.int64)
    dist = np.zeros(nbin)
    
    for i in range(nbin):
        mask = (distances > x[i]) & (distances <= x[i + 1])
        pairs[i] = np.count_nonzero(mask)
        dist[i] = np.sum(distances[mask])
        
    return pairs, dist