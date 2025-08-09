#!/usr/bin/env python3

"""
missing_data_handler.py - Advanced missing data handling and interpolation
Usage: python missing_data_handler.py <coverage_file> <output_dir> [--method hmm]
"""

import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats, interpolate
from sklearn.impute import KNNImputer
from sklearn.ensemble import RandomForestRegressor
import warnings
warnings.filterwarnings('ignore')

class MissingDataHandler:
    """
    Advanced missing data handler for CNV detection with multiple interpolation methods
    """
    
    def __init__(self, method='hmm', window_size=3, min_coverage=5):
        """
        Initialize missing data handler
        
        Args:
            method: Interpolation method ('hmm', 'linear', 'knn', 'rf', 'median')
            window_size: Window size for local interpolation methods
            min_coverage: Minimum coverage to consider as valid
        """
        self.method = method
        self.window_size = window_size
        self.min_coverage = min_coverage
        
    def analyze_missing_patterns(self, coverage_data):
        """Analyze missing data patterns in coverage data"""
        print("Analyzing missing data patterns...")
        
        analysis_results = {}
        
        # Overall missing data statistics
        total_positions = coverage_data.size
        missing_positions = coverage_data.isna().sum().sum()
        missing_fraction = missing_positions / total_positions
        
        analysis_results['total_positions'] = total_positions
        analysis_results['missing_positions'] = missing_positions
        analysis_results['missing_fraction'] = missing_fraction
        
        # Per-sample missing data
        sample_missing = coverage_data.isna().sum(axis=1)
        analysis_results['samples_with_missing'] = (sample_missing > 0).sum()
        analysis_results['mean_missing_per_sample'] = sample_missing.mean()
        analysis_results['max_missing_per_sample'] = sample_missing.max()
        
        # Per-exon missing data
        exon_missing = coverage_data.isna().sum(axis=0)
        analysis_results['exons_with_missing'] = (exon_missing > 0).sum()
        analysis_results['mean_missing_per_exon'] = exon_missing.mean()
        analysis_results['max_missing_per_exon'] = exon_missing.max()
        
        # Identify problematic exons (>50% missing)
        problematic_exons = exon_missing[exon_missing > len(coverage_data) * 0.5].index.tolist()
        analysis_results['problematic_exons'] = problematic_exons
        
        # Missing data patterns
        analysis_results['complete_samples'] = (sample_missing == 0).sum()
        analysis_results['partially_missing_samples'] = ((sample_missing > 0) & (sample_missing < len(coverage_data.columns))).sum()
        analysis_results['completely_missing_samples'] = (sample_missing == len(coverage_data.columns)).sum()
        
        # Consecutive missing regions
        consecutive_missing = self._find_consecutive_missing_regions(coverage_data)
        analysis_results['consecutive_missing_regions'] = consecutive_missing
        
        print(f"  Total positions: {total_positions}")
        print(f"  Missing positions: {missing_positions} ({missing_fraction:.1%})")
        print(f"  Samples with missing data: {analysis_results['samples_with_missing']}")
        print(f"  Exons with missing data: {analysis_results['exons_with_missing']}")
        print(f"  Problematic exons (>50% missing): {len(problematic_exons)}")
        
        return analysis_results
    
    def _find_consecutive_missing_regions(self, coverage_data):
        """Find consecutive missing regions per sample"""
        consecutive_regions = []
        
        for sample_id in coverage_data.index:
            sample_data = coverage_data.loc[sample_id]
            missing_mask = sample_data.isna()
            
            # Find consecutive missing stretches
            in_missing_region = False
            region_start = None
            
            for i, is_missing in enumerate(missing_mask):
                if is_missing and not in_missing_region:
                    # Start of missing region
                    region_start = i
                    in_missing_region = True
                elif not is_missing and in_missing_region:
                    # End of missing region
                    region_length = i - region_start
                    consecutive_regions.append({
                        'sample_id': sample_id,
                        'start_exon': region_start,
                        'end_exon': i - 1,
                        'length': region_length
                    })
                    in_missing_region = False
            
            # Handle case where missing region extends to end
            if in_missing_region:
                region_length = len(missing_mask) - region_start
                consecutive_regions.append({
                    'sample_id': sample_id,
                    'start_exon': region_start,
                    'end_exon': len(missing_mask) - 1,
                    'length': region_length
                })
        
        return consecutive_regions
    
    def _linear_interpolation(self, coverage_data):
        """Linear interpolation for missing values"""
        interpolated_data = coverage_data.copy()
        
        for sample_id in coverage_data.index:
            sample_data = coverage_data.loc[sample_id].values
            
            # Create position array
            positions = np.arange(len(sample_data))
            valid_mask = ~np.isnan(sample_data)
            
            if valid_mask.sum() >= 2:  # Need at least 2 points for interpolation
                # Interpolate missing values
                f = interpolate.interp1d(positions[valid_mask], sample_data[valid_mask], 
                                       kind='linear', bounds_error=False, 
                                       fill_value='extrapolate')
                interpolated_values = f(positions)
                
                # Only replace missing values
                missing_mask = np.isnan(sample_data)
                sample_data[missing_mask] = interpolated_values[missing_mask]
                
                interpolated_data.loc[sample_id] = sample_data
        
        return interpolated_data
    
    def _knn_imputation(self, coverage_data, n_neighbors=5):
        """K-Nearest Neighbors imputation"""
        imputer = KNNImputer(n_neighbors=min(n_neighbors, len(coverage_data) - 1))
        
        try:
            imputed_array = imputer.fit_transform(coverage_data)
            interpolated_data = pd.DataFrame(imputed_array, 
                                           index=coverage_data.index, 
                                           columns=coverage_data.columns)
            return interpolated_data
        except Exception as e:
            print(f"Warning: KNN imputation failed ({e}), falling back to linear interpolation")
            return self._linear_interpolation(coverage_data)
    
    def _random_forest_imputation(self, coverage_data):
        """Random Forest-based imputation"""
        interpolated_data = coverage_data.copy()
        
        for col in coverage_data.columns:
            col_data = coverage_data[col]
            missing_mask = col_data.isna()
