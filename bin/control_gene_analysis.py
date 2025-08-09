#!/usr/bin/env python3

"""
control_gene_analysis.py - Control gene selection and normalization for CNV detection
Usage: python control_gene_analysis.py <coverage_file> <control_genes_bed> <output_dir> [--method auto]
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
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

# Default stable control genes for normalization
DEFAULT_CONTROL_GENES = [
    'GAPDH', 'ACTB', 'HPRT1', 'GUSB', 'TBP', 'YWHAZ', 
    'PPIA', 'RPL13A', 'SDHA', 'UBC'
]

def read_control_genes_bed(bed_file):
    """Read control genes BED file and return gene coordinates."""
    control_genes = {}
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 4:
                    chrom, start, end, gene_name = parts[:4]
                    
                    if gene_name not in control_genes:
                        control_genes[gene_name] = []
                    
                    control_genes[gene_name].append({
                        'chrom': chrom,
                        'start': int(start),
                        'end': int(end),
                        'length': int(end) - int(start)
                    })
    except FileNotFoundError:
        print(f"Warning: Control genes BED file not found: {bed_file}")
        return {}
    
    return control_genes

def calculate_gene_stability(coverage_data, min_samples=5):
    """Calculate stability metrics for potential control genes."""
    stability_metrics = {}
    
    # Group coverage by gene
    for gene in coverage_data.columns:
        if gene == 'sample_id':
            continue
            
        gene_coverage = coverage_data[gene].dropna()
        
        if len(gene_coverage) < min_samples:
            continue
        
        # Calculate stability metrics
        cv = gene_coverage.std() / gene_coverage.mean() if gene_coverage.mean() > 0 else float('inf')
        mad = np.median(np.abs(gene_coverage - gene_coverage.median()))
        q75, q25 = np.percentile(gene_coverage, [75, 25])
        iqr_cv = (q75 - q25) / gene_coverage.median() if gene_coverage.median() > 0 else float('inf')
        
        # Outlier detection
        z_scores = np.abs(stats.zscore(gene_coverage))
        outlier_fraction = np.sum(z_scores > 2.5) / len(gene_coverage)
        
        stability_metrics[gene] = {
            'mean_coverage': gene_coverage.mean(),
            'median_coverage': gene_coverage.median(),
            'cv': cv,  # Coefficient of variation
            'mad': mad,  # Median absolute deviation
            'iqr_cv': iqr_cv,  # IQR coefficient of variation
            'outlier_fraction': outlier_fraction,
            'n_samples': len(gene_coverage),
            'stability_score': 1.0 / (1.0 + cv + iqr_cv + outlier_fraction)  # Higher is more stable
        }
    
    return stability_metrics

def select_control_genes(stability_metrics, target_genes=None, n_controls=5, 
                        min_coverage=10, max_cv=0.3, max_outliers=0.1):
    """Select optimal control genes based on stability metrics."""
    
    # Filter genes based on quality criteria
    candidate_genes = {}
    for gene, metrics in stability_metrics.items():
        if (metrics['mean_coverage'] >= min_coverage and 
            metrics['cv'] <= max_cv and 
            metrics['outlier_fraction'] <= max_outliers and
            metrics['n_samples'] >= 3):
            candidate_genes[gene] = metrics
    
    # Exclude target genes from control selection
    if target_genes:
        for target in target_genes:
            candidate_genes.pop(target, None)
    
    if len(candidate_genes) == 0:
        print("Warning: No genes meet control gene criteria")
        return []
    
    # Sort by stability score and select top N
    sorted_genes = sorted(candidate_genes.items(), 
                         key=lambda x: x[1]['stability_score'], 
                         reverse=True)
    
    selected_genes = [gene for gene, _ in sorted_genes[:n_controls]]
    
    return selected_genes

def calculate_normalization_factors(coverage_data, control_genes, method='geometric_mean'):
    """Calculate sample-specific normalization factors using control genes."""
    
    if not control_genes:
        print("Warning: No control genes available for normalization")
        return pd.Series(1.0, index=coverage_data.index)
    
    # Extract control gene coverage
    available_controls = [gene for gene in control_genes if gene in coverage_data.columns]
    
    if not available_controls:
        print("Warning: No control genes found in coverage data")
        return pd.Series(1.0, index=coverage_data.index)
    
    control_coverage = coverage_data[available_controls]
    
    # Calculate normalization factors
    normalization_factors = pd.Series(index=coverage_data.index, dtype=float)
    
    if method == 'geometric_mean':
        # Geometric mean of control genes (similar to DESeq2 size factors)
        log_coverage = np.log(control_coverage + 1)  # Add pseudocount
        geo_means = np.exp(log_coverage.mean(axis=1))
        
        # Calculate size factors relative to overall geometric mean
        overall_geo_mean = np.exp(log_coverage.values.mean())
        normalization_factors = geo_means / overall_geo_mean
        
    elif method == 'median_ratio':
        # Median ratio method (similar to DESeq2)
        # Calculate gene-wise geometric means across samples
        gene_geo_means = np.exp(np.log(control_coverage + 1).mean(axis=0))
        
        # Calculate ratios for each sample
        ratios = control_coverage.div(gene_geo_means, axis=1)
        
        # Use median ratio as normalization factor
        normalization_factors = ratios.median(axis=1)
        
    elif method == 'tmm':
        # Trimmed mean of M-values (TMM) - simplified version
        ref_sample_idx = control_coverage.sum(axis=1).idxmax()  # Sample with highest total coverage
        ref_coverage = control_coverage.loc[ref_sample_idx]
        
        for sample_idx in control_coverage.index:
            sample_coverage = control_coverage.loc[sample_idx]
            
            # Calculate M-values (log-fold changes)
            m_values = np.log2((sample_coverage + 1) / (ref_coverage + 1))
            
            # Calculate A-values (average log intensities)  
            a_values = 0.5 * np.log2((sample_coverage + 1) * (ref_coverage + 1))
            
            # Trim extreme values
            valid_idx = ~(np.isnan(m_values) | np.isnan(a_values) | np.isinf(m_values) | np.isinf(a_values))
            if valid_idx.sum() > 0:
                m_trimmed = np.trim_zeros(np.sort(m_values[valid_idx]))[len(m_values[valid_idx])//10:-len(m_values[valid_idx])//10]
                if len(m_trimmed) > 0:
                    tmm_factor = np.mean(m_trimmed)
                    normalization_factors.loc[sample_idx] = 2 ** tmm_factor
                else:
                    normalization_factors.loc[sample_idx] = 1.0
            else:
                normalization_factors.loc[sample_idx] = 1.0
    
    else:
        # Simple mean normalization
        control_means = control_coverage.mean(axis=1)
        overall_mean = control_means.mean()
        normalization_factors = control_means / overall_mean
    
    # Handle edge cases
    normalization_factors = normalization_factors.fillna(1.0)
    normalization_factors[normalization_factors <= 0] = 1.0
    normalization_factors[np.isinf(normalization_factors)] = 1.0
    
    return normalization_factors

def detect_batch_effects(coverage_data, sample_info, control_genes):
    """Detect potential batch effects using control genes."""
    
    if not control_genes or 'batch' not in sample_info.columns:
        return {'batch_detected': False, 'batch_pvalue': 1.0, 'batch_effect_size': 0.0}
    
    # Extract control gene coverage
    available_controls = [gene for gene in control_genes if gene in coverage_data.columns]
    if not available_controls:
        return {'batch_detected': False, 'batch_pvalue': 1.0, 'batch_effect_size': 0.0}
    
    control_coverage = coverage_data[available_controls]
    
    # Merge with sample info
    merged_data = control_coverage.merge(sample_info[['sample_id', 'batch']], 
                                        left_index=True, right_on='sample_id', how='inner')
    
    if merged_data['batch'].nunique() < 2:
        return {'batch_detected': False, 'batch_pvalue': 1.0, 'batch_effect_size': 0.0}
    
    # PCA analysis to detect batch effects
    scaler = StandardScaler()
    scaled_coverage = scaler.fit_transform(control_coverage.fillna(0))
    
    pca = PCA(n_components=min(2, scaled_coverage.shape[1]))
    pca_result = pca.fit_transform(scaled_coverage)
    
    # Test for batch effect on first PC
    batches = sample_info.loc[sample_info['sample_id'].isin(coverage_data.index), 'batch'].values
    
    if len(np.unique(batches)) >= 2:
        # ANOVA test for batch effect
        batch_groups = [pca_result[batches == batch, 0] for batch in np.unique(batches)]
        f_stat, p_value = stats.f_oneway(*batch_groups)
        
        # Effect size calculation (eta squared)
        ss_between = sum(len(group) * (np.mean(group) - np.mean(pca_result[:, 0]))**2 for group in batch_groups)
        ss_total = np.sum((pca_result[:, 0] - np.mean(pca_result[:, 0]))**2)
        eta_squared = ss_between / ss_total if ss_total > 0 else 0
        
        batch_detected = p_value < 0.05 and eta_squared > 0.1
        
        return {
            'batch_detected': batch_detected,
            'batch_pvalue': p_value,
            'batch_effect_size': eta_squared,
            'pc1_variance_explained': pca.explained_variance_ratio_[0]
        }
    
    return {'batch_detected': False, 'batch_pvalue': 1.0, 'batch_effect_size': 0.0}

def apply_control_normalization(coverage_data, normalization_factors, target_genes):
    """Apply control gene normalization to target genes."""
    
    normalized_data = coverage_data.copy()
    
    # Normalize target genes using control-derived factors
    for gene in target_genes:
        if gene in normalized_data.columns:
            normalized_data[gene] = normalized_data[gene] / normalization_factors
    
    return normalized_data

def create_control_analysis_plots(stability_metrics, selected_controls, normalization_factors, 
                                 coverage_data, output_dir):
    """Create visualization plots for control gene analysis."""
    
    plot_dir = Path(output_dir) / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    # Plot 1: Control gene stability metrics
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    genes = list(stability_metrics.keys())
    cvs = [stability_metrics[gene]['cv'] for gene in genes]
    stability_scores = [stability_metrics[gene]['stability_score'] for gene in genes]
    outlier_fractions = [stability_metrics[gene]['outlier_fraction'] for gene in genes]
    mean_coverages = [stability_metrics[gene]['mean_coverage'] for gene in genes]
    
    # CV distribution
    axes[0, 0].hist(cvs, bins=20, alpha=0.7, edgecolor='black')
    axes[0, 0].axvline(x=0.3, color='red', linestyle='--', label='CV threshold')
    axes[0, 0].set_xlabel('Coefficient of Variation')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Control Gene CV Distribution')
    axes[0, 0].legend()
    
    # Stability score vs CV
    colors = ['red' if gene in selected_controls else 'blue' for gene in genes]
    axes[0, 1].scatter(cvs, stability_scores, c=colors, alpha=0.7)
    axes[0, 1].set_xlabel('Coefficient of Variation')
    axes[0, 1].set_ylabel('Stability Score')
    axes[0, 1].set_title('Gene Stability Analysis')
    
    # Mean coverage vs outlier fraction
    axes[1, 0].scatter(mean_coverages, outlier_fractions, c=colors, alpha=0.7)
    axes[1, 0].set_xlabel('Mean Coverage')
    axes[1, 0].set_ylabel('Outlier Fraction')
    axes[1, 0].set_title('Coverage vs Outlier Analysis')
    axes[1, 0].axhline(y=0.1, color='red', linestyle='--', alpha=0.7)
    
    # Normalization factors distribution
    axes[1, 1].hist(normalization_factors, bins=20, alpha=0.7, edgecolor='black')
    axes[1, 1].axvline(x=1.0, color='red', linestyle='-', alpha=0.7, label='No normalization')
    axes[1, 1].set_xlabel('Normalization Factor')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Sample Normalization Factors')
    axes[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'control_gene_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Before/after normalization comparison
    if selected_controls:
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # Select a representative control gene for visualization
        repr_control = selected_controls[0] if selected_controls[0] in coverage_data.columns else None
        
        if repr_control:
            original_coverage = coverage_data[repr_control]
            normalized_coverage = coverage_data[repr_control] / normalization_factors
            
            # Before normalization
            axes[0].hist(original_coverage, bins=20, alpha=0.7, edgecolor='black', color='blue')
            axes[0].set_xlabel('Coverage')
            axes[0].set_ylabel('Frequency')
            axes[0].set_title(f'Before Normalization - {repr_control}')
            axes[0].axvline(x=original_coverage.median(), color='red', linestyle='--', 
                           label=f'Median: {original_coverage.median():.1f}')
            axes[0].legend()
            
            # After normalization
            axes[1].hist(normalized_coverage, bins=20, alpha=0.7, edgecolor='black', color='green')
            axes[1].set_xlabel('Normalized Coverage')
            axes[1].set_ylabel('Frequency')
            axes[1].set_title(f'After Normalization - {repr_control}')
            axes[1].axvline(x=normalized_coverage.median(), color='red', linestyle='--',
                           label=f'Median: {normalized_coverage.median():.1f}')
            axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(plot_dir / 'normalization_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Control analysis plots saved to: {plot_dir}")

def main():
    parser = argparse.ArgumentParser(description='Control gene analysis and normalization')
    parser.add_argument('coverage_file', help='Coverage data file')
    parser.add_argument('control_genes_bed', help='Control genes BED file')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--method', choices=['geometric_mean', 'median_ratio', 'tmm', 'mean'], 
                       default='geometric_mean', help='Normalization method')
    parser.add_argument('--target-genes', nargs='+', default=['SMN1', 'SMN2'], 
                       help='Target genes to normalize')
    parser.add_argument('--n-controls', type=int, default=5, 
                       help='Number of control genes to select')
    parser.add_argument('--min-coverage', type=float, default=10, 
                       help='Minimum mean coverage for control genes')
    parser.add_argument('--max-cv', type=float, default=0.3, 
                       help='Maximum CV for control genes')
    parser.add_argument('--no-plots', action='store_true', help='Skip creating plots')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading coverage data...")
    
    # Read coverage data
    try:
        if args.coverage_file.endswith('.txt') or args.coverage_file.endswith('.tsv'):
            coverage_df = pd.read_csv(args.coverage_file, sep='\t', index_col=0)
        else:
            coverage_df = pd.read_csv(args.coverage_file, index_col=0)
        print(f"Loaded coverage data: {coverage_df.shape}")
    except Exception as e:
        print(f"Error loading coverage file: {e}")
        sys.exit(1)
    
    # Read control genes
    control_genes_coords = read_control_genes_bed(args.control_genes_bed)
    potential_controls = list(control_genes_coords.keys()) + DEFAULT_CONTROL_GENES
    potential_controls = [gene for gene in potential_controls if gene in coverage_df.columns]
    
    print(f"Found {len(potential_controls)} potential control genes in data")
    
    # Calculate stability metrics
    print("Calculating gene stability metrics...")
    stability_metrics = calculate_gene_stability(coverage_df)
    
    # Select optimal control genes
    print("Selecting optimal control genes...")
    selected_controls = select_control_genes(
        stability_metrics, 
        target_genes=args.target_genes,
        n_controls=args.n_controls,
        min_coverage=args.min_coverage,
        max_cv=args.max_cv
    )
    
    print(f"Selected control genes: {selected_controls}")
    
    # Calculate normalization factors
    print(f"Calculating normalization factors using {args.method} method...")
    normalization_factors = calculate_normalization_factors(
        coverage_df, selected_controls, method=args.method
    )
    
    # Apply normalization to target genes
    print("Applying control gene normalization...")
    normalized_data = apply_control_normalization(
        coverage_df, normalization_factors, args.target_genes
    )
    
    # Save results
    print("Saving results...")
    
    # Control gene selection results
    control_selection_file = output_dir / 'control_selection.txt'
    with open(control_selection_file, 'w') as f:
        f.write("# Control Gene Selection Results\n")
        f.write("gene\tstability_score\tmean_coverage\tcv\toutlier_fraction\tn_samples\n")
        for gene in selected_controls:
            if gene in stability_metrics:
                metrics = stability_metrics[gene]
                f.write(f"{gene}\t{metrics['stability_score']:.4f}\t{metrics['mean_coverage']:.2f}\t"
                       f"{metrics['cv']:.4f}\t{metrics['outlier_fraction']:.4f}\t{metrics['n_samples']}\n")
    
    # Normalization factors
    normalization_file = output_dir / 'normalization_factors.txt'
    norm_df = pd.DataFrame({
        'sample_id': normalization_factors.index,
        'normalization_factor': normalization_factors.values,
        'method': args.method
    })
    norm_df.to_csv(normalization_file, sep='\t', index=False)
    
    # Normalized coverage data
    normalized_file = output_dir / 'normalized_coverage.txt'
    normalized_data.to_csv(normalized_file, sep='\t')
    
    # Stability metrics for all genes
    stability_file = output_dir / 'gene_stability_metrics.txt'
    stability_df = pd.DataFrame.from_dict(stability_metrics, orient='index')
    stability_df.index.name = 'gene'
    stability_df.to_csv(stability_file, sep='\t')
    
    # Create plots
    if not args.no_plots:
        try:
            create_control_analysis_plots(
                stability_metrics, selected_controls, normalization_factors, 
                coverage_df, output_dir
            )
        except Exception as e:
            print(f"Warning: Could not create plots: {e}")
    
    # Summary statistics
    print("\nControl Gene Analysis Summary:")
    print(f"  Total genes analyzed: {len(stability_metrics)}")
    print(f"  Control genes selected: {len(selected_controls)}")
    print(f"  Normalization method: {args.method}")
    print(f"  Normalization factor range: {normalization_factors.min():.3f} - {normalization_factors.max():.3f}")
    print(f"  Normalization factor median: {normalization_factors.median():.3f}")
    
    print(f"\nResults saved to: {output_dir}")
    print(f"  Control selection: {control_selection_file}")
    print(f"  Normalization factors: {normalization_file}")
    print(f"  Normalized coverage: {normalized_file}")
    print(f"  Stability metrics: {stability_file}")

if __name__ == "__main__":
    main()
