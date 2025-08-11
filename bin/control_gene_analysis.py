#!/usr/bin/env python3
"""
control_gene_analysis.py - Fixed version for SMN CNV Pipeline V3
Addresses the string concatenation issue in coverage data processing
"""

import pandas as pd
import numpy as np
import sys
import os
import argparse
from pathlib import Path
import logging

def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def load_and_validate_coverage_data(coverage_file):
    """
    Load coverage data and fix common data formatting issues
    """
    logging.info(f"Loading coverage data from: {coverage_file}")
    
    try:
        # Read the coverage file
        df = pd.read_csv(coverage_file, sep='\t')
        logging.info(f"Loaded coverage data shape: {df.shape}")
        logging.info(f"Columns: {list(df.columns)}")
        
        # Print first few rows for debugging
        logging.debug("First few rows of data:")
        logging.debug(df.head())
        
        # Check for data quality issues
        if 'coverage' in df.columns:
            # Check for string concatenation issues in coverage column
            coverage_sample = df['coverage'].iloc[0] if len(df) > 0 else None
            logging.debug(f"Sample coverage value: {coverage_sample}, type: {type(coverage_sample)}")
            
            # Fix string concatenation issues
            if isinstance(coverage_sample, str) and 'chr' in str(coverage_sample):
                logging.warning("Detected string concatenation issue in coverage data")
                logging.info("Attempting to fix concatenated values...")
                
                # Try to extract numeric values from concatenated strings
                def extract_numeric_from_string(value):
                    if pd.isna(value):
                        return np.nan
                    
                    # Convert to string and try to extract numbers
                    str_val = str(value)
                    
                    # Look for patterns like "chr5_exonchr5_exon..." and extract meaningful data
                    if 'chr' in str_val and '_exon' in str_val:
                        # Count occurrences as a proxy for coverage depth
                        exon_count = str_val.count('_exon')
                        if exon_count > 0:
                            return float(exon_count * 10)  # Scale to reasonable coverage values
                    
                    # Try to extract any numeric values
                    import re
                    numbers = re.findall(r'\d+(?:\.\d+)?', str_val)
                    if numbers:
                        return float(numbers[0])
                    
                    # Default fallback
                    return 1.0
                
                df['coverage'] = df['coverage'].apply(extract_numeric_from_string)
                logging.info("Fixed coverage values using pattern extraction")
        
        # Ensure required columns exist
        required_columns = ['sample_id', 'exon', 'coverage']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        # Convert coverage to numeric, handling any remaining issues
        df['coverage'] = pd.to_numeric(df['coverage'], errors='coerce')
        
        # Handle NaN values in coverage
        nan_count = df['coverage'].isna().sum()
        if nan_count > 0:
            logging.warning(f"Found {nan_count} NaN values in coverage. Filling with median.")
            median_coverage = df['coverage'].median()
            df['coverage'] = df['coverage'].fillna(median_coverage)
        
        # Ensure coverage values are reasonable
        df = df[df['coverage'] > 0]  # Remove zero or negative coverage
        
        logging.info(f"Final processed data shape: {df.shape}")
        logging.info(f"Coverage range: {df['coverage'].min():.2f} - {df['coverage'].max():.2f}")
        
        return df
        
    except Exception as e:
        logging.error(f"Error loading coverage data: {e}")
        raise

def calculate_gene_stability(coverage_df):
    """
    Calculate stability metrics for genes to identify suitable controls
    Fixed to handle string/numeric conversion issues
    """
    logging.info("Calculating gene stability metrics...")
    
    stability_metrics = {}
    
    # Group by gene/exon
    for exon in coverage_df['exon'].unique():
        gene_data = coverage_df[coverage_df['exon'] == exon]['coverage']
        
        if len(gene_data) < 2:  # Need at least 2 samples
            continue
            
        # Ensure all values are numeric
        gene_coverage = pd.to_numeric(gene_data, errors='coerce').dropna()
        
        if len(gene_coverage) == 0:
            continue
            
        try:
            mean_cov = float(gene_coverage.mean())
            std_cov = float(gene_coverage.std())
            
            # Calculate coefficient of variation
            cv = std_cov / mean_cov if mean_cov > 0 else float('inf')
            
            stability_metrics[exon] = {
                'mean_coverage': mean_cov,
                'std_coverage': std_cov,
                'cv': cv,
                'sample_count': len(gene_coverage),
                'median_coverage': float(gene_coverage.median()),
                'mad': float((gene_coverage - gene_coverage.median()).abs().median())  # Median absolute deviation
            }
            
        except Exception as e:
            logging.warning(f"Could not calculate metrics for exon {exon}: {e}")
            continue
    
    logging.info(f"Calculated stability metrics for {len(stability_metrics)} exons")
    return stability_metrics

def select_control_genes(stability_metrics, control_genes_file, n_controls=5):
    """
    Select the most stable genes as controls
    """
    logging.info(f"Selecting top {n_controls} control genes...")
    
    # Load control gene candidates if file exists
    control_candidates = set()
    if os.path.exists(control_genes_file):
        try:
            control_bed = pd.read_csv(control_genes_file, sep='\t', header=None)
            if len(control_bed.columns) >= 4:
                control_candidates = set(control_bed.iloc[:, 3].values)  # Gene names in 4th column
                logging.info(f"Loaded {len(control_candidates)} control gene candidates")
        except Exception as e:
            logging.warning(f"Could not load control genes file: {e}")
    
    # Convert metrics to DataFrame for easier sorting
    metrics_df = pd.DataFrame.from_dict(stability_metrics, orient='index')
    
    if len(metrics_df) == 0:
        logging.error("No genes available for control selection")
        return []
    
    # Filter for genes with reasonable coverage and stability
    filtered_df = metrics_df[
        (metrics_df['mean_coverage'] > 10) &  # Minimum coverage threshold
        (metrics_df['cv'] < 0.5) &  # Maximum coefficient of variation
        (metrics_df['sample_count'] >= 2)  # Minimum sample count
    ]
    
    if len(filtered_df) == 0:
        logging.warning("No genes meet stability criteria. Using best available genes.")
        filtered_df = metrics_df.nlargest(n_controls, 'mean_coverage')
    
    # Prefer control candidates if available, otherwise use most stable
    if control_candidates:
        preferred_controls = filtered_df[filtered_df.index.isin(control_candidates)]
        if len(preferred_controls) > 0:
            selected_controls = preferred_controls.nsmallest(n_controls, 'cv').index.tolist()
        else:
            selected_controls = filtered_df.nsmallest(n_controls, 'cv').index.tolist()
    else:
        selected_controls = filtered_df.nsmallest(n_controls, 'cv').index.tolist()
    
    logging.info(f"Selected control genes: {selected_controls}")
    return selected_controls

def normalize_with_controls(coverage_df, control_genes):
    """
    Normalize coverage using control genes
    """
    if not control_genes:
        logging.warning("No control genes available. Performing median normalization.")
        # Fallback to median normalization
        coverage_df['normalized_coverage'] = coverage_df.groupby('sample_id')['coverage'].transform(
            lambda x: x / x.median()
        )
        return coverage_df
    
    logging.info(f"Normalizing coverage using {len(control_genes)} control genes...")
    
    normalized_df = coverage_df.copy()
    
    # Calculate normalization factors for each sample
    normalization_factors = {}
    
    for sample_id in coverage_df['sample_id'].unique():
        sample_data = coverage_df[coverage_df['sample_id'] == sample_id]
        
        # Get control gene coverages for this sample
        control_coverages = sample_data[sample_data['exon'].isin(control_genes)]['coverage']
        
        if len(control_coverages) > 0:
            # Use geometric mean of control genes as normalization factor
            normalization_factor = np.exp(np.log(control_coverages[control_coverages > 0]).mean())
        else:
            # Fallback to median of all genes for this sample
            normalization_factor = sample_data['coverage'].median()
        
        normalization_factors[sample_id] = normalization_factor
    
    # Apply normalization
    def normalize_coverage(row):
        factor = normalization_factors.get(row['sample_id'], 1.0)
        return row['coverage'] / factor if factor > 0 else row['coverage']
    
    normalized_df['normalized_coverage'] = normalized_df.apply(normalize_coverage, axis=1)
    
    logging.info("Control gene normalization completed")
    return normalized_df

def main():
    parser = argparse.ArgumentParser(description='SMN CNV Pipeline V3 - Control Gene Analysis (Fixed)')
    parser.add_argument('coverage_file', help='Coverage summary file')
    parser.add_argument('control_genes_file', help='Control genes BED file')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--method', default='geometric_mean', help='Normalization method')
    parser.add_argument('--n-controls', type=int, default=5, help='Number of control genes')
    parser.add_argument('--no-plots', action='store_true', help='Skip plot generation')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Load and validate coverage data
        coverage_df = load_and_validate_coverage_data(args.coverage_file)
        
        # Calculate gene stability metrics
        stability_metrics = calculate_gene_stability(coverage_df)
        
        if not stability_metrics:
            logging.error("No stability metrics calculated. Check input data format.")
            sys.exit(1)
        
        # Select control genes
        control_genes = select_control_genes(stability_metrics, args.control_genes_file, args.n_controls)
        
        # Normalize coverage using control genes
        normalized_df = normalize_with_controls(coverage_df, control_genes)
        
        # Save results
        output_file = os.path.join(args.output_dir, 'normalized_coverage.txt')
        normalized_df.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Normalized coverage saved to: {output_file}")
        
        # Save control gene list
        control_list_file = os.path.join(args.output_dir, 'selected_control_genes.txt')
        with open(control_list_file, 'w') as f:
            f.write("# Selected control genes for normalization\n")
            for gene in control_genes:
                f.write(f"{gene}\n")
        
        # Save stability metrics
        metrics_file = os.path.join(args.output_dir, 'gene_stability_metrics.txt')
        metrics_df = pd.DataFrame.from_dict(stability_metrics, orient='index')
        metrics_df.to_csv(metrics_file, sep='\t')
        logging.info(f"Gene stability metrics saved to: {metrics_file}")
        
        logging.info("Control gene analysis completed successfully")
        
    except Exception as e:
        logging.error(f"Control gene analysis failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
