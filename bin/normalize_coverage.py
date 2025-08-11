#!/usr/bin/env python3

"""
normalize_coverage_v3.py - Enhanced coverage normalization with V3 features
Usage: python normalize_coverage_v3.py <coverage_file> <sample_info_file> <output_file> [--control-factors control_factors.txt]
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
import yaml
import warnings
warnings.filterwarnings('ignore')

def load_control_normalization_factors(control_factors_file):
    """Load pre-computed control gene normalization factors"""
    try:
        factors_df = pd.read_csv(control_factors_file, sep='\t', index_col='sample_id')
        factors = factors_df['normalization_factor'].to_dict()
        print(f"Loaded control normalization factors for {len(factors)} samples")
        return factors
    except Exception as e:
        print(f"Warning: Could not load control factors ({e}). Using standard normalization.")
        return {}

def load_v3_parameters(params_file):
    """Load V3 segmentation parameters from YAML config"""
    default_params = {
        'quality_control': {
            'coverage': {
                'min_depth': 10,
                'min_breadth': 0.8,
                'max_cv': 2.0
            }
        },
        'control_normalization': {
            'method': 'geometric_mean'
        }
    }
    
    if params_file and os.path.exists(params_file):
        try:
            with open(params_file, 'r') as f:
                params = yaml.safe_load(f)
            return params
        except Exception as e:
            print(f"Warning: Could not load parameters file ({e}). Using defaults.")
    
    return default_params

def apply_control_normalization(coverage_df, control_factors):
    """Apply pre-computed control gene normalization factors"""
    if not control_factors:
        return coverage_df
    
    normalized_df = coverage_df.copy()
    
    for sample_id in coverage_df.index:
        if sample_id in control_factors:
            factor = control_factors[sample_id]
            if factor > 0 and not np.isnan(factor):
                normalized_df.loc[sample_id] = normalized_df.loc[sample_id] / factor
            else:
                print(f"Warning: Invalid normalization factor for {sample_id}: {factor}")
    
    return normalized_df

def enhanced_reference_stats(coverage_df, reference_samples, method='robust'):
    """Calculate enhanced reference statistics with multiple methods"""
    ref_stats = {}
    
    # Filter for reference samples
    ref_df = coverage_df[coverage_df.index.isin(reference_samples)]
    
    if ref_df.empty:
        print("Warning: No reference samples found in coverage data!")
        return ref_stats
    
    print(f"Calculating reference statistics using {len(reference_samples)} reference samples")
    
    # Calculate statistics for each exon
    exons = coverage_df.columns
    
    for exon in exons:
        exon_data = ref_df[exon].dropna()
        
        if len(exon_data) > 0:
            if method == 'robust':
                # Robust statistics using median and MAD
                median_val = np.median(exon_data)
                mad = np.median(np.abs(exon_data - median_val))
                
                # Convert MAD to std equivalent
                mad_std = mad * 1.4826  # MAD to std conversion factor
                
                # Use MAD-based outlier detection
                outlier_threshold = 2.5 * mad_std
                outliers = np.abs(exon_data - median_val) > outlier_threshold
                
                # Filter outliers
                filtered_data = exon_data[~outliers]
                
                if len(filtered_data) > 1:
                    ref_stats[exon] = {
                        'mean': np.mean(filtered_data),
                        'std': np.std(filtered_data, ddof=1),
                        'median': np.median(filtered_data),
                        'mad': np.median(np.abs(filtered_data - np.median(filtered_data))),
                        'n_samples': len(filtered_data),
                        'n_outliers': len(exon_data) - len(filtered_data),
                        'min': np.min(filtered_data),
                        'max': np.max(filtered_data),
                        'method': 'robust'
                    }
                else:
                    ref_stats[exon] = {
                        'mean': median_val,
                        'std': mad_std if mad_std > 0 else 1.0,
                        'median': median_val,
                        'mad': mad,
                        'n_samples': len(exon_data),
                        'n_outliers': 0,
                        'min': median_val,
                        'max': median_val,
                        'method': 'robust_fallback'
                    }
            
            elif method == 'trimmed':
                # Trimmed mean approach
                trim_fraction = 0.1  # Trim 10% from each end
                sorted_data = np.sort(exon_data)
                n_trim = int(len(sorted_data) * trim_fraction)
                
                if n_trim > 0 and len(sorted_data) > 2 * n_trim:
                    trimmed_data = sorted_data[n_trim:-n_trim]
                else:
                    trimmed_data = sorted_data
                
                ref_stats[exon] = {
                    'mean': np.mean(trimmed_data),
                    'std': np.std(trimmed_data, ddof=1) if len(trimmed_data) > 1 else 1.0,
                    'median': np.median(trimmed_data),
                    'mad': np.median(np.abs(trimmed_data - np.median(trimmed_data))),
                    'n_samples': len(trimmed_data),
                    'n_outliers': len(exon_data) - len(trimmed_data),
                    'min': np.min(trimmed_data),
                    'max': np.max(trimmed_data),
                    'method': 'trimmed'
                }
            
            else:
                # Standard approach (IQR-based outlier removal)
                Q1 = np.percentile(exon_data, 25)
                Q3 = np.percentile(exon_data, 75)
                IQR = Q3 - Q1
                lower_bound = Q1 - 1.5 * IQR
                upper_bound = Q3 + 1.5 * IQR
                
                filtered_data = exon_data[(exon_data >= lower_bound) & (exon_data <= upper_bound)]
                
                if len(filtered_data) > 1:
                    ref_stats[exon] = {
                        'mean': np.mean(filtered_data),
                        'std': np.std(filtered_data, ddof=1),
                        'median': np.median(filtered_data),
                        'mad': np.median(np.abs(filtered_data - np.median(filtered_data))),
                        'n_samples': len(filtered_data),
                        'n_outliers': len(exon_data) - len(filtered_data),
                        'min': np.min(filtered_data),
                        'max': np.max(filtered_data),
                        'method': 'standard'
                    }
                else:
                    ref_stats[exon] = {
                        'mean': np.mean(exon_data),
                        'std': np.std(exon_data, ddof=1) if len(exon_data) > 1 else 1.0,
                        'median': np.median(exon_data),
                        'mad': np.median(np.abs(exon_data - np.median(exon_data))),
                        'n_samples': len(exon_data),
                        'n_outliers': 0,
                        'min': np.min(exon_data),
                        'max': np.max(exon_data),
                        'method': 'standard_fallback'
                    }
        else:
            print(f"Warning: No coverage data for exon {exon} in reference samples")
    
    return ref_stats

def calculate_enhanced_z_scores(coverage_df, ref_stats, method='standard'):
    """Calculate Z-scores with enhanced methods"""
    z_score_data = []
    
    for sample_id in coverage_df.index:
        for exon in coverage_df.columns:
            coverage = coverage_df.loc[sample_id, exon]
            
            if exon in ref_stats:
                ref_mean = ref_stats[exon]['mean']
                ref_std = ref_stats[exon]['std']
                ref_median = ref_stats[exon]['median']
                ref_mad = ref_stats[exon]['mad']
                
                # Calculate different types of normalized scores
                if method == 'robust':
                    # Use median and MAD for robust z-score
                    mad_std = ref_mad * 1.4826
                    if mad_std > 0:
                        robust_z = (coverage - ref_median) / mad_std
                    else:
                        robust_z = 0.0
                    z_score = robust_z
                    
                elif method == 'modified':
                    # Modified z-score using MAD
                    if ref_mad > 0:
                        modified_z = 0.6745 * (coverage - ref_median) / ref_mad
                    else:
                        modified_z = 0.0
                    z_score = modified_z
                    
                else:
                    # Standard z-score
                    if ref_std > 0:
                        standard_z = (coverage - ref_mean) / ref_std
                    else:
                        standard_z = 0.0
                    z_score = standard_z
                
                # Calculate additional quality metrics
                percentile_rank = stats.percentileofscore(
                    [ref_stats[exon]['min'], ref_stats[exon]['max']], coverage
                ) / 100.0
                
                z_score_data.append({
                    'sample_id': sample_id,
                    'exon': exon,
                    'raw_coverage': coverage,
                    'ref_mean': ref_mean,
                    'ref_std': ref_std,
                    'ref_median': ref_median,
                    'ref_mad': ref_mad,
                    'z_score': z_score,
                    'percentile_rank': percentile_rank,
                    'ref_n_samples': ref_stats[exon]['n_samples'],
                    'normalization_method': method
                })
            else:
                print(f"Warning: No reference statistics for exon {exon}")
                z_score_data.append({
                    'sample_id': sample_id,
                    'exon': exon,
                    'raw_coverage': coverage,
                    'ref_mean': np.nan,
                    'ref_std': np.nan,
                    'ref_median': np.nan,
                    'ref_mad': np.nan,
                    'z_score': np.nan,
                    'percentile_rank': np.nan,
                    'ref_n_samples': 0,
                    'normalization_method': method
                })
    
    return pd.DataFrame(z_score_data)

def assess_coverage_quality(coverage_df, params):
    """Assess coverage quality using V3 parameters"""
    quality_assessment = {}
    qc_params = params.get('quality_control', {}).get('coverage', {})
    
    min_depth = qc_params.get('min_depth', 10)
    min_breadth = qc_params.get('min_breadth', 0.8)
    max_cv = qc_params.get('max_cv', 2.0)
    
    for sample_id in coverage_df.index:
        sample_data = coverage_df.loc[sample_id]
        valid_data = sample_data.dropna()
        
        if len(valid_data) > 0:
            mean_coverage = valid_data.mean()
            cv = valid_data.std() / mean_coverage if mean_coverage > 0 else float('inf')
            breadth = len(valid_data) / len(sample_data)
            low_coverage_fraction = (valid_data < min_depth).sum() / len(valid_data)
            
            # Quality flags
            flags = []
            if mean_coverage < min_depth:
                flags.append('low_depth')
            if breadth < min_breadth:
                flags.append('low_breadth')
            if cv > max_cv:
                flags.append('high_variability')
            if low_coverage_fraction > 0.3:
                flags.append('many_low_coverage_positions')
            
            # Overall quality score (0-1)
            depth_score = min(1.0, mean_coverage / (2 * min_depth))
            breadth_score = breadth
            cv_score = max(0.0, 1.0 - cv / max_cv)
            
            quality_score = (depth_score + breadth_score + cv_score) / 3.0
            
            quality_assessment[sample_id] = {
                'mean_coverage': mean_coverage,
                'cv': cv,
                'breadth': breadth,
                'low_coverage_fraction': low_coverage_fraction,
                'quality_score': quality_score,
                'quality_flags': flags,
                'pass_qc': len(flags) == 0
            }
        else:
            quality_assessment[sample_id] = {
                'mean_coverage': 0.0,
                'cv': float('inf'),
                'breadth': 0.0,
                'low_coverage_fraction': 1.0,
                'quality_score': 0.0,
                'quality_flags': ['no_valid_data'],
                'pass_qc': False
            }
    
    return quality_assessment

def create_v3_normalization_plots(z_scores_df, ref_stats, quality_assessment, output_dir):
    """Create enhanced visualization plots for V3 normalization"""
    
    plot_dir = Path(output_dir) / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    # Plot 1: Enhanced Z-score distributions with quality overlay
    plt.figure(figsize=(15, 10))
    exons = z_scores_df['exon'].unique()
    
    n_exons = len(exons)
    n_cols = min(2, n_exons)
    n_rows = (n_exons + n_cols - 1) // n_cols
    
    for i, exon in enumerate(sorted(exons)):
        plt.subplot(n_rows, n_cols, i+1)
        exon_data = z_scores_df[z_scores_df['exon'] == exon]
        
        # Color points by quality
        colors = []
        for sample_id in exon_data['sample_id']:
            if sample_id in quality_assessment:
                if quality_assessment[sample_id]['pass_qc']:
                    colors.append('green')
                elif quality_assessment[sample_id]['quality_score'] > 0.5:
                    colors.append('orange')
                else:
                    colors.append('red')
            else:
                colors.append('gray')
        
        plt.scatter(range(len(exon_data)), exon_data['z_score'], c=colors, alpha=0.7)
        
        # Add reference lines
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        plt.axhline(y=-2, color='red', linestyle='--', alpha=0.7, label='CN boundary')
        plt.axhline(y=2, color='red', linestyle='--', alpha=0.7)
        
        plt.title(f'{exon} Z-scores (Quality Colored)')
        plt.xlabel('Sample Index')
        plt.ylabel('Z-score')
        
        if i == 0:
            # Add legend only to first subplot
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='green', label='Pass QC'),
                Patch(facecolor='orange', label='Warning'),
                Patch(facecolor='red', label='Fail QC')
            ]
            plt.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'v3_z_score_quality_overview.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Quality assessment summary
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    quality_df = pd.DataFrame.from_dict(quality_assessment, orient='index')
    
    # Quality score distribution
    axes[0, 0].hist(quality_df['quality_score'], bins=20, alpha=0.7, edgecolor='black')
    axes[0, 0].axvline(x=0.7, color='red', linestyle='--', label='QC threshold')
    axes[0, 0].set_xlabel('Quality Score')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Sample Quality Score Distribution')
    axes[0, 0].legend()
    
    # Coverage depth vs CV
    axes[0, 1].scatter(quality_df['mean_coverage'], quality_df['cv'], alpha=0.7)
    axes[0, 1].set_xlabel('Mean Coverage')
    axes[0, 1].set_ylabel('Coefficient of Variation')
    axes[0, 1].set_title('Coverage Depth vs Variability')
    axes[0, 1].set_yscale('log')
    
    # Breadth distribution
    axes[1, 0].hist(quality_df['breadth'], bins=20, alpha=0.7, edgecolor='black')
    axes[1, 0].axvline(x=0.8, color='red', linestyle='--', label='Min breadth')
    axes[1, 0].set_xlabel('Coverage Breadth')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Coverage Breadth Distribution')
    axes[1, 0].legend()
    
    # QC pass rates
    pass_counts = quality_df['pass_qc'].value_counts()
    axes[1, 1].pie(pass_counts.values, labels=['Fail QC', 'Pass QC'], autopct='%1.1f%%',
                  colors=['red', 'green'], alpha=0.7)
    axes[1, 1].set_title('Quality Control Pass/Fail Rates')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'v3_quality_assessment.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Reference statistics comparison
    if ref_stats:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        exons = list(ref_stats.keys())
        ref_means = [ref_stats[exon]['mean'] for exon in exons]
        ref_stds = [ref_stats[exon]['std'] for exon in exons]
        ref_cvs = [ref_stats[exon]['std'] / ref_stats[exon]['mean'] 
                   if ref_stats[exon]['mean'] > 0 else 0 for exon in exons]
        n_samples = [ref_stats[exon]['n_samples'] for exon in exons]
        
        # Reference mean coverage by exon
        bars = axes[0, 0].bar(range(len(exons)), ref_means, alpha=0.7, color='blue')
        axes[0, 0].set_xlabel('Exons')
        axes[0, 0].set_ylabel('Mean Coverage')
        axes[0, 0].set_title('Reference Mean Coverage by Exon')
        axes[0, 0].set_xticks(range(len(exons)))
        axes[0, 0].set_xticklabels([str(e) for e in exons], rotation=45)
        
        # Error bars for standard deviation
        axes[0, 0].errorbar(range(len(exons)), ref_means, yerr=ref_stds, 
                           fmt='none', color='red', alpha=0.7)
        
        # Reference CV by exon
        axes[0, 1].bar(range(len(exons)), ref_cvs, alpha=0.7, color='orange')
        axes[0, 1].axhline(y=0.3, color='red', linestyle='--', label='High CV threshold')
        axes[0, 1].set_xlabel('Exons')
        axes[0, 1].set_ylabel('Coefficient of Variation')
        axes[0, 1].set_title('Reference CV by Exon')
        axes[0, 1].set_xticks(range(len(exons)))
        axes[0, 1].set_xticklabels([str(e) for e in exons], rotation=45)
        axes[0, 1].legend()
        
        # Number of reference samples per exon
        axes[1, 0].bar(range(len(exons)), n_samples, alpha=0.7, color='green')
        axes[1, 0].set_xlabel('Exons')
        axes[1, 0].set_ylabel('Number of Reference Samples')
        axes[1, 0].set_title('Reference Sample Count by Exon')
        axes[1, 0].set_xticks(range(len(exons)))
        axes[1, 0].set_xticklabels([str(e) for e in exons], rotation=45)
        
        # Mean vs STD relationship
        axes[1, 1].scatter(ref_means, ref_stds, alpha=0.7)
        axes[1, 1].set_xlabel('Reference Mean')
        axes[1, 1].set_ylabel('Reference Standard Deviation')
        axes[1, 1].set_title('Mean-Variance Relationship')
        
        # Add exon labels
        for i, exon in enumerate(exons):
            axes[1, 1].annotate(str(exon), (ref_means[i], ref_stds[i]), 
                               xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        plt.tight_layout()
        plt.savefig(plot_dir / 'v3_reference_statistics.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"V3 normalization plots saved to: {plot_dir}")

def main():
    parser = argparse.ArgumentParser(description='Enhanced V3 coverage normalization')
    parser.add_argument('coverage_file', help='Coverage data file')
    parser.add_argument('sample_info_file', help='Sample information file')
    parser.add_argument('output_file', help='Output file for Z-scores')
    parser.add_argument('--control-factors', help='Control gene normalization factors file')
    parser.add_argument('--params-file', help='V3 parameters YAML file')
    parser.add_argument('--method', choices=['standard', 'robust', 'modified'], 
                       default='robust', help='Z-score calculation method')
    parser.add_argument('--ref-method', choices=['standard', 'robust', 'trimmed'],
                       default='robust', help='Reference statistics method')
    parser.add_argument('--no-plots', action='store_true', help='Skip creating plots')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load V3 parameters
    v3_params = load_v3_parameters(args.params_file)
    
    print("Loading input files...")
    
    # Read coverage data
    try:
        coverage_df = pd.read_csv(args.coverage_file, sep='\t', index_col=0)
        print(f"Loaded coverage data for {len(coverage_df)} samples, {len(coverage_df.columns)} exons")
    except Exception as e:
        print(f"Error reading coverage file: {e}")
        sys.exit(1)
    
    # Read sample information
    try:
        sample_df = pd.read_csv(args.sample_info_file, sep='\t')
        samples_info = {}
        for _, row in sample_df.iterrows():
            samples_info[row['sample_id']] = {
                'bam_path': row.get('bam_path', ''),
                'sample_type': row.get('sample_type', 'unknown')
            }
    except Exception as e:
        print(f"Warning: Could not read sample info file ({e}). Using auto-detection.")
        samples_info = {}
        for sample_id in coverage_df.index:
            if any(keyword in sample_id.lower() for keyword in ['ref', 'control', 'normal']):
                sample_type = 'reference'
            else:
                sample_type = 'test'
            samples_info[sample_id] = {'sample_type': sample_type, 'bam_path': ''}
    
    # Apply control gene normalization if available
    control_factors = load_control_normalization_factors(args.control_factors) if args.control_factors else {}
    if control_factors:
        print("Applying control gene normalization...")
        coverage_df = apply_control_normalization(coverage_df, control_factors)
    
    # Identify reference samples
    reference_samples = [sid for sid, info in samples_info.items() 
                        if info['sample_type'] == 'reference']
    
    print(f"Found {len(reference_samples)} reference samples")
    
    if len(reference_samples) < 2:
        print("Warning: Very few reference samples. Results may be unreliable.")
    
    # Calculate enhanced reference statistics
    print(f"Calculating reference statistics using {args.ref_method} method...")
    ref_stats = enhanced_reference_stats(coverage_df, reference_samples, method=args.ref_method)
    
    # Calculate enhanced Z-scores
    print(f"Calculating Z-scores using {args.method} method...")
    z_scores_df = calculate_enhanced_z_scores(coverage_df, ref_stats, method=args.method)
    
    # Add sample type information
    z_scores_df['sample_type'] = z_scores_df['sample_id'].map(
        lambda x: samples_info.get(x, {}).get('sample_type', 'unknown')
    )
    
    # Assess coverage quality
    print("Assessing coverage quality...")
    quality_assessment = assess_coverage_quality(coverage_df, v3_params)
    
    # Add quality information to Z-scores
    z_scores_df['quality_score'] = z_scores_df['sample_id'].map(
        lambda x: quality_assessment.get(x, {}).get('quality_score', 0.0)
    )
    z_scores_df['pass_qc'] = z_scores_df['sample_id'].map(
        lambda x: quality_assessment.get(x, {}).get('pass_qc', False)
    )
    
    # Save results
    z_scores_df.to_csv(args.output_file, index=False, sep='\t')
    
    # Save enhanced reference statistics
    ref_stats_file = args.output_file.replace('.txt', '_ref_stats.txt')
    ref_stats_df = pd.DataFrame.from_dict(ref_stats, orient='index').reset_index()
    ref_stats_df.rename(columns={'index': 'exon'}, inplace=True)
    ref_stats_df.to_csv(ref_stats_file, index=False, sep='\t')
    
    # Save quality assessment
    quality_file = args.output_file.replace('.txt', '_quality_assessment.txt')
    quality_df = pd.DataFrame.from_dict(quality_assessment, orient='index').reset_index()
    quality_df.rename(columns={'index': 'sample_id'}, inplace=True)
    quality_df.to_csv(quality_file, index=False, sep='\t')
    
    # Create enhanced plots
    if not args.no_plots:
        try:
            create_v3_normalization_plots(z_scores_df, ref_stats, quality_assessment, output_dir)
        except Exception as e:
            print(f"Warning: Could not create plots: {e}")
    
    # Print V3 summary
    print(f"\nV3 Enhanced Normalization Summary:")
    print(f"  Method: {args.method} Z-scores with {args.ref_method} reference statistics")
    print(f"  Control normalization: {'Applied' if control_factors else 'Not applied'}")
    print(f"  Total samples: {len(z_scores_df['sample_id'].unique())}")
    print(f"  Reference samples: {len(reference_samples)}")
    
    # Quality summary
    qc_pass_count = sum(1 for qa in quality_assessment.values() if qa['pass_qc'])
    print(f"  Quality control pass rate: {qc_pass_count}/{len(quality_assessment)} ({100*qc_pass_count/len(quality_assessment):.1f}%)")
    
    # Coverage summary
    mean_quality_score = np.mean([qa['quality_score'] for qa in quality_assessment.values()])
    print(f"  Mean quality score: {mean_quality_score:.3f}")
    
    print(f"\nResults saved to:")
    print(f"  Z-scores: {args.output_file}")
    print(f"  Reference statistics: {ref_stats_file}")
    print(f"  Quality assessment: {quality_file}")

if __name__ == "__main__":
    main()
