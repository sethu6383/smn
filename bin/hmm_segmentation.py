#!/usr/bin/env python3

"""
hmm_segmentation.py - Hidden Markov Model segmentation for CNV detection
Usage: python hmm_segmentation.py <normalized_coverage_file> <output_dir> [--states 5]
"""

import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.mixture import GaussianMixture
import warnings
warnings.filterwarnings('ignore')

class CNV_HMM:
    """
    Hidden Markov Model for CNV detection with missing data handling
    """
    
    def __init__(self, n_states=5, max_iter=100, tol=1e-4, random_state=42):
        """
        Initialize HMM parameters
        
        Args:
            n_states: Number of copy number states (typically 5: 0,1,2,3,4+)
            max_iter: Maximum number of EM iterations
            tol: Convergence tolerance
            random_state: Random seed for reproducibility
        """
        self.n_states = n_states
        self.max_iter = max_iter
        self.tol = tol
        self.random_state = random_state
        
        # Initialize parameters
        self.transition_matrix = None
        self.emission_means = None
        self.emission_variances = None
        self.initial_probs = None
        
        # Copy number mapping
        self.state_to_cn = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}  # State -> Copy Number
        
    def _initialize_parameters(self, data):
        """Initialize HMM parameters from data"""
        np.random.seed(self.random_state)
        
        # Remove missing values for parameter initialization
        valid_data = data[~np.isnan(data)]
        
        if len(valid_data) < self.n_states:
            # Not enough data, use default parameters
            self.emission_means = np.array([0.0, 0.5, 1.0, 1.5, 2.0])[:self.n_states]
            self.emission_variances = np.ones(self.n_states) * 0.25
        else:
            # Use Gaussian Mixture Model to initialize emission parameters
            gmm = GaussianMixture(n_components=self.n_states, random_state=self.random_state)
            gmm.fit(valid_data.reshape(-1, 1))
            
            # Sort by means to ensure proper ordering
            sorted_indices = np.argsort(gmm.means_.flatten())
            self.emission_means = gmm.means_.flatten()[sorted_indices]
            self.emission_variances = gmm.covariances_.flatten()[sorted_indices]
        
        # Initialize transition matrix (distance-dependent)
        self.transition_matrix = self._create_transition_matrix()
        
        # Initialize state probabilities (favor normal copy number)
        self.initial_probs = np.ones(self.n_states) / self.n_states
        if self.n_states >= 3:
            self.initial_probs[2] = 0.5  # Favor CN=2 (normal)
            self.initial_probs /= self.initial_probs.sum()
    
    def _create_transition_matrix(self, persistence_prob=0.9):
        """
        Create distance-dependent transition matrix
        Closer states have higher transition probabilities
        """
        trans_matrix = np.zeros((self.n_states, self.n_states))
        
        for i in range(self.n_states):
            for j in range(self.n_states):
                if i == j:
                    # Stay in same state (persistence)
                    trans_matrix[i, j] = persistence_prob
                else:
                    # Transition probability decreases with state distance
                    distance = abs(i - j)
                    trans_matrix[i, j] = (1 - persistence_prob) * np.exp(-distance)
            
            # Normalize rows
            trans_matrix[i, :] /= trans_matrix[i, :].sum()
        
        return trans_matrix
    
    def _emission_probability(self, observation, state):
        """Calculate emission probability for observation given state"""
        if np.isnan(observation):
            return 1.0  # Equal probability for missing observations
        
        mean = self.emission_means[state]
        var = self.emission_variances[state]
        
        if var <= 0:
            var = 1e-6  # Avoid division by zero
        
        # Gaussian emission probability
        prob = stats.norm.pdf(observation, loc=mean, scale=np.sqrt(var))
        return max(prob, 1e-10)  # Avoid zero probabilities
    
    def _forward_algorithm(self, observations):
        """Forward algorithm for HMM"""
        T = len(observations)
        alpha = np.zeros((T, self.n_states))
        
        # Initialize
        for s in range(self.n_states):
            alpha[0, s] = (self.initial_probs[s] * 
                          self._emission_probability(observations[0], s))
        
        # Forward pass
        for t in range(1, T):
            for s in range(self.n_states):
                alpha[t, s] = (np.sum(alpha[t-1, :] * self.transition_matrix[:, s]) *
                              self._emission_probability(observations[t], s))
            
            # Normalize to prevent underflow
            alpha[t, :] /= (alpha[t, :].sum() + 1e-10)
        
        return alpha
    
    def _backward_algorithm(self, observations):
        """Backward algorithm for HMM"""
        T = len(observations)
        beta = np.zeros((T, self.n_states))
        
        # Initialize
        beta[T-1, :] = 1.0
        
        # Backward pass
        for t in range(T-2, -1, -1):
            for s in range(self.n_states):
                beta[t, s] = np.sum(
                    self.transition_matrix[s, :] * beta[t+1, :] *
                    [self._emission_probability(observations[t+1], s_next) 
                     for s_next in range(self.n_states)]
                )
            
            # Normalize to prevent underflow
            beta[t, :] /= (beta[
