#!/usr/bin/env python3
"""
missing_data_handler.py - Handle missing data and interpolation for SMN CNV Pipeline V3
"""

import pandas as pd
import numpy as np
import sys
import os
import argparse
import logging
from pathlib import Path
from scipy import interpolate
from sklearn.impute import KNNImputer
from sklearn.ensemble import RandomForestRegressor

def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def load_coverage_data(coverage_file):
    """Load coverage data with robust error handling"""
    logging.info(f"Loading coverage data from: {coverage_file}")
    
    if not os.path.exists(coverage_file):
        raise FileNotFoundError(f"Coverage file not found: {coverage_file}")
    
    try:
        df = pd.read_csv(coverage_file, sep='\t')
        logging.info(f"Loaded coverage data shape: {df.shape}")
        
        # Ensure required columns exist
        required_cols = ['sample_id', 'exon', 'coverage']
        if 'normalized_coverage' in df.columns:
            required_cols.append('normalized_coverage')
            
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Use normalized_coverage if available, otherwise use coverage
        if 'normalized_coverage' in df.columns:
            df['working_coverage'] = pd.to_numeric(df['normalized_coverage'], errors='coerce')
        else:
            df['working_coverage'] = pd.to_numeric(df['coverage'], errors='coerce')
            # Create normalized_coverage column for consistency
            df['normalized_coverage'] = df['working_coverage']
        
        # Handle any NaN values
        nan_count = df['working_coverage'].isna().sum()
        if nan_count > 0:
            logging.warning(f"Found {nan_count} NaN values in coverage data")
            # Fill with median of non-NaN values
            median_val = df['working_coverage'].median()
            if pd.isna(median_val):
                median_val = 1.0  # Fallback
            df['working_coverage'] = df['working_coverage'].fillna(median_val)
            df['normalized_coverage'] = df['normalized_coverage'].fillna(median_val)
        
        return df
        
    except Exception as e:
        logging.error(f"Error loading coverage data: {e}")
        raise

def detect_missing_data_patterns(df):
    """Detect patterns of missing data in the coverage matrix"""
    logging.info("Analyzing missing data patterns...")
    
    # Create pivot table for analysis
    pivot_df = df.pivot(index='sample_id', columns='exon', values='working_coverage')
    
    # Calculate missing data statistics
    total_cells = pivot_df.size
    missing_cells = pivot_df.isna().sum().sum()
    missing_percentage = (missing_cells / total_cells) * 100
    
    logging.info(f"Total data points: {total_cells}")
    logging.info(f"Missing data points: {missing_cells}")
    logging.info(f"Missing data percentage: {missing_percentage:.2f}%")
    
    # Identify samples and exons with high missing data
    sample_missing = pivot_df.isna().sum(axis=1)
    exon_missing = pivot_df.isna().sum(axis=0)
    
    high_missing_samples = sample_missing[sample_missing > len(pivot_df.columns) * 0.5]
    high_missing_exons = exon_missing[exon_missing > len(pivot_df.index) * 0.5]
    
    if len(high_missing_samples) > 0:
        logging.warning(f"Samples with >50% missing data: {list(high_missing_samples.index)}")
    
    if len(high_missing_exons) > 0:
        logging.warning(f"Exons with >50% missing data: {list(high_missing_exons.index)}")
    
    return {
        'total_cells': total_cells,
        'missing_cells': missing_cells,
        'missing_percentage': missing_percentage,
        'pivot_df': pivot_df,
        'sample_missing': sample_missing,
        'exon_missing': exon_missing
    }

def interpolate_linear(df, analysis_results):
    """Linear interpolation for missing data"""
    logging.info("Performing linear interpolation...")
    
    pivot_df = analysis_results['pivot_df'].copy()
    
    # Interpolate along exons (columns) for each sample
    for sample_id in pivot_df.index:
        sample_data = pivot_df.loc[sample_id, :]
        if sample_data.isna().any():
            # Linear interpolation
            interpolated = sample_data.interpolate(method='linear')
            # Fill any remaining NaN values at the edges
            interpolated = interpolated.fillna(method='bfill').fillna(method='ffill')
            pivot_df.loc[sample_id, :] = interpolated
    
    return pivot_df

def interpolate_knn(df, analysis_results, n_neighbors=3):
    """KNN imputation for missing data"""
    logging.info(f"Performing KNN interpolation with {n_neighbors} neighbors...")
    
    pivot_df = analysis_results['pivot_df'].copy()
    
    # Use KNN imputer
    imputer = KNNImputer(n_neighbors=min(n_neighbors, len(pivot_df.index)-1))
    imputed_values = imputer.fit_transform(pivot_df)
    
    # Create new dataframe with imputed values
    imputed_df = pd.DataFrame(
        imputed_values,
        index=pivot_df.index,
        columns=pivot_df.columns
    )
    
    return imputed_df

def interpolate_random_forest(df, analysis_results):
    """Random Forest imputation for missing data"""
    logging.info("Performing Random Forest interpolation...")
    
    pivot_df = analysis_results['pivot_df'].copy()
    
    # For each column with missing data, use other columns to predict
    for col in pivot_df.columns:
        if pivot_df[col].isna().any():
            # Get complete cases for training
            train_data = pivot_df.dropna()
            
            if len(train_data) == 0:
                # Fallback to median imputation
                pivot_df[col] = pivot_df[col].fillna(pivot_df[col].median())
                continue
            
            # Train Random Forest
            X_train = train_data.drop(columns=[col])
            y_train = train_data[col]
            
            if len(X_train.columns) == 0:
                # No features to train on, use median
                pivot_df[col] = pivot_df[col].fillna(pivot_df[col].median())
                continue
            
            rf = RandomForestRegressor(n_estimators=10, random_state=42, max_depth=3)
            rf.fit(X_train, y_train)
            
            # Predict missing values
            missing_idx = pivot_df[col].isna()
            if missing_idx.any():
                X_missing = pivot_df.loc[missing_idx, :].drop(columns=[col])
                # Handle any remaining NaN in features
                X_missing = X_missing.fillna(X_missing.median())
                predictions = rf.predict(X_missing)
                pivot_df.loc[missing_idx, col] = predictions
    
    return pivot_df

def interpolate_median(df, analysis_results):
    """Simple median imputation"""
    logging.info("Performing median interpolation...")
    
    pivot_df = analysis_results['pivot_df'].copy()
    
    # Fill missing values with column median
    for col in pivot_df.columns:
        median_val = pivot_df[col].median()
        if pd.isna(median_val):
            median_val = pivot_df[col].mean()  # Fallback to mean
        if pd.isna(median_val):
            median_val = 1.0  # Final fallback
        pivot_df[col] = pivot_df[col].fillna(median_val)
    
    return pivot_df

def interpolate_hmm_simple(df, analysis_results):
    """Simple HMM-like interpolation using local context"""
    logging.info("Performing HMM-like interpolation...")
    
    pivot_df = analysis_results['pivot_df'].copy()
    
    # For each sample, use neighboring values and sample context
    for sample_id in pivot_df.index:
        sample_data = pivot_df.loc[sample_id, :]
        
        if sample_data.isna().any():
            # First pass: linear interpolation
            interpolated = sample_data.interpolate(method='linear')
            
            # Second pass: use median of similar samples for remaining NaN
            if interpolated.isna().any():
                # Find most similar samples (based on available data)
                other_samples = pivot_df.drop(sample_id)
                
                for missing_col in interpolated[interpolated.isna()].index:
                    # Use median of this column from other samples
                    col_values = other_samples[missing_col].dropna()
                    if len(col_values) > 0:
                        interpolated[missing_col] = col_values.median()
                    else:
                        interpolated[missing_col] = 1.0  # Fallback
            
            pivot_df.loc[sample_id, :] = interpolated
    
    return pivot_df

def perform_interpolation(df, method, analysis_results):
    """Perform the specified interpolation method"""
    logging.info(f"Performing interpolation using method: {method}")
    
    method = method.lower()
    
    if method == 'linear':
        return interpolate_linear(df, analysis_results)
    elif method == 'knn':
        return interpolate_knn(df, analysis_results)
    elif method == 'rf' or method == 'random_forest':
        return interpolate_random_forest(df, analysis_results)
    elif method == 'median':
        return interpolate_median(df, analysis_results)
    elif method == 'hmm':
        return interpolate_hmm_simple(df, analysis_results)
    else:
        logging.warning(f"Unknown interpolation method: {method}. Using median.")
        return interpolate_median(df, analysis_results)

def convert_back_to_long_format(pivot_df, original_df):
    """Convert pivot table back to long format"""
    logging.info("Converting back to long format...")
    
    # Melt the pivot table
    melted = pivot_df.reset_index().melt(
        id_vars=['sample_id'],
        var_name='exon',
        value_name='interpolated_coverage'
    )
    
    # Merge with original data to preserve other columns
    result_df = original_df.merge(
        melted,
        on=['sample_id', 'exon'],
        how='left'
    )
    
    # Update the coverage columns
    result_df['coverage'] = result_df['interpolated_coverage']
    result_df['normalized_coverage'] = result_df['interpolated_coverage']
    
    # Add flag for interpolated values
    original_pivot = original_df.pivot(index='sample_id', columns='exon', values='working_coverage')
    was_missing = original_pivot.isna()
    was_missing_melted = was_missing.reset_index().melt(
        id_vars=['sample_id'],
        var_name='exon',
        value_name='was_interpolated'
    )
    
    result_df = result_df.merge(was_missing_melted, on=['sample_id', 'exon'], how='left')
    
    return result_df

def main():
    parser = argparse.ArgumentParser(description='Missing Data Handler for SMN CNV Pipeline V3')
    parser.add_argument('coverage_file', help='Input coverage file')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--method', default='hmm', 
                       choices=['linear', 'knn', 'rf', 'random_forest', 'median', 'hmm'],
                       help='Interpolation method')
    parser.add_argument('--validate', action='store_true', help='Perform validation')
    parser.add_argument('--no-plots', action='store_true', help='Skip plots')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Load coverage data
        df = load_coverage_data(args.coverage_file)
        
        # Analyze missing data patterns
        analysis_results = detect_missing_data_patterns(df)
        
        # Perform interpolation
        interpolated_pivot = perform_interpolation(df, args.method, analysis_results)
        
        # Convert back to long format
        result_df = convert_back_to_long_format(interpolated_pivot, df)
        
        # Save results
        output_file = os.path.join(args.output_dir, 'interpolated_coverage.txt')
        result_df.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Interpolated coverage saved to: {output_file}")
        
        # Save missing data analysis
        analysis_file = os.path.join(args.output_dir, 'missing_data_analysis.txt')
        with open(analysis_file, 'w') as f:
            f.write(f"Missing Data Analysis Report\n")
            f.write(f"============================\n\n")
            f.write(f"Total data points: {analysis_results['total_cells']}\n")
            f.write(f"Missing data points: {analysis_results['missing_cells']}\n")
            f.write(f"Missing percentage: {analysis_results['missing_percentage']:.2f}%\n")
            f.write(f"Interpolation method: {args.method}\n")
        
        logging.info("Missing data handling completed successfully")
        
    except Exception as e:
        logging.error(f"Missing data handling failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
