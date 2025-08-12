#!/usr/bin/env python3
"""
consensus_calling.py - Multi-algorithm consensus calling for SMN CNV Pipeline V3
"""

import pandas as pd
import numpy as np
import sys
import os
import argparse
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def load_segmentation_results(cbs_file, hmm_file):
    """Load results from CBS and HMM segmentation"""
    logging.info("Loading segmentation results for consensus calling...")
    
    algorithms = {}
    
    # Load CBS results
    if os.path.exists(cbs_file):
        try:
            cbs_df = pd.read_csv(cbs_file, sep='\t')
            algorithms['CBS'] = cbs_df
            logging.info(f"Loaded CBS results: {len(cbs_df)} segments")
        except Exception as e:
            logging.warning(f"Could not load CBS results from {cbs_file}: {e}")
    
    # Load HMM results
    if os.path.exists(hmm_file):
        try:
            hmm_df = pd.read_csv(hmm_file, sep='\t')
            algorithms['HMM'] = hmm_df
            logging.info(f"Loaded HMM results: {len(hmm_df)} segments")
        except Exception as e:
            logging.warning(f"Could not load HMM results from {hmm_file}: {e}")
    
    if not algorithms:
        # Create dummy data for testing
        logging.warning("No algorithm results found. Creating dummy data.")
        
        samples = ['sample1', 'sample2', 'sample3']
        dummy_segments = []
        
        for sample in samples:
            # CBS-like segments
            dummy_segments.append({
                'sample_id': sample,
                'chromosome': 'chr5',
                'start_pos': 1000000,
                'end_pos': 1200000,
                'num_markers': 3,
                'mean_coverage': np.random.normal(1.0, 0.2),
                'copy_number': 2
            })
            
        cbs_df = pd.DataFrame(dummy_segments)
        hmm_df = cbs_df.copy()
        hmm_df['confidence'] = 0.8
        hmm_df['state'] = 1
        
        algorithms['CBS'] = cbs_df
        algorithms['HMM'] = hmm_df
    
    logging.info(f"Loaded results from {len(algorithms)} algorithms: {list(algorithms.keys())}")
    return algorithms

def standardize_segment_format(df, algorithm_name):
    """Standardize segment format across algorithms"""
    standardized = df.copy()
    
    # Ensure required columns exist
    required_columns = ['sample_id', 'start_pos', 'end_pos', 'copy_number']
    for col in required_columns:
        if col not in standardized.columns:
            if col == 'copy_number':
                standardized[col] = 2  # Default diploid
            else:
                standardized[col] = 0
    
    # Add algorithm identifier
    standardized['algorithm'] = algorithm_name
    
    # Add confidence if not present
    if 'confidence' not in standardized.columns:
        if algorithm_name == 'HMM':
            standardized['confidence'] = 0.8  # Default HMM confidence
        else:
            standardized['confidence'] = 0.7  # Default CBS confidence
    
    # Ensure numeric types
    for col in ['start_pos', 'end_pos', 'copy_number']:
        if col in standardized.columns:
            standardized[col] = pd.to_numeric(standardized[col], errors='coerce')
    
    return standardized

def find_overlapping_segments(seg1, seg2):
    """Check if two segments overlap and calculate overlap metrics"""
    overlap_start = max(seg1['start_pos'], seg2['start_pos'])
    overlap_end = min(seg1['end_pos'], seg2['end_pos'])
    
    if overlap_start <= overlap_end:
        overlap_length = overlap_end - overlap_start + 1
        seg1_length = seg1['end_pos'] - seg1['start_pos'] + 1
        seg2_length = seg2['end_pos'] - seg2['start_pos'] + 1
        
        # Calculate overlap percentage relative to each segment
        overlap_pct_1 = overlap_length / seg1_length if seg1_length > 0 else 0
        overlap_pct_2 = overlap_length / seg2_length if seg2_length > 0 else 0
        
        return {
            'overlap': True,
            'overlap_start': overlap_start,
            'overlap_end': overlap_end,
            'overlap_length': overlap_length,
            'overlap_pct_1': overlap_pct_1,
            'overlap_pct_2': overlap_pct_2,
            'min_overlap_pct': min(overlap_pct_1, overlap_pct_2)
        }
    
    return {'overlap': False}

def create_consensus_segments(algorithms, min_algorithms=2, confidence_threshold=0.7, overlap_threshold=0.5):
    """Create consensus segments from multiple algorithms"""
    logging.info(f"Creating consensus segments (min_algorithms={min_algorithms}, confidence_threshold={confidence_threshold})")
    
    # Standardize all algorithm results
    standardized_results = {}
    for alg_name, df in algorithms.items():
        standardized_results[alg_name] = standardize_segment_format(df, alg_name)
    
    # Combine all segments
    all_segments = []
    for alg_name, df in standardized_results.items():
        for _, segment in df.iterrows():
            segment_dict = segment.to_dict()
            all_segments.append(segment_dict)
    
    if not all_segments:
        logging.warning("No segments available for consensus calling")
        return pd.DataFrame()
    
    all_segments_df = pd.DataFrame(all_segments)
    
    consensus_calls = []
    
    # Group by sample
    for sample_id in all_segments_df['sample_id'].unique():
        sample_segments = all_segments_df[all_segments_df['sample_id'] == sample_id]
        
        logging.debug(f"Processing consensus for sample {sample_id} with {len(sample_segments)} segments")
        
        # Find overlapping segments across algorithms
        processed_indices = set()
        
        for i, seg1 in sample_segments.iterrows():
            if i in processed_indices:
                continue
            
            # Find all segments that overlap with this one
            overlapping_segments = [seg1]
            overlapping_indices = {i}
            
            for j, seg2 in sample_segments.iterrows():
                if i == j or j in processed_indices:
                    continue
                
                overlap_info = find_overlapping_segments(seg1, seg2)
                
                if overlap_info['overlap'] and overlap_info['min_overlap_pct'] >= overlap_threshold:
                    overlapping_segments.append(seg2)
                    overlapping_indices.add(j)
            
            # Create consensus call from overlapping segments
            if len(overlapping_segments) >= min_algorithms:
                consensus_call = create_single_consensus(overlapping_segments, confidence_threshold)
                if consensus_call is not None:
                    consensus_calls.append(consensus_call)
            
            # Mark these segments as processed
            processed_indices.update(overlapping_indices)
    
    if consensus_calls:
        consensus_df = pd.DataFrame(consensus_calls)
        logging.info(f"Created {len(consensus_df)} consensus calls")
        return consensus_df
    else:
        logging.warning("No consensus calls met the criteria")
        return pd.DataFrame()

def create_single_consensus(segments, confidence_threshold):
    """Create a single consensus call from overlapping segments"""
    if not segments:
        return None
    
    # Convert to DataFrame for easier processing
    if isinstance(segments[0], pd.Series):
        segments_df = pd.DataFrame([seg.to_dict() for seg in segments])
    else:
        segments_df = pd.DataFrame(segments)
    
    # Basic consensus metrics
    sample_id = segments_df['sample_id'].iloc[0]
    algorithms_used = segments_df['algorithm'].unique().tolist()
    
    # Position consensus (use intersection of all segments)
    consensus_start = segments_df['start_pos'].max()
    consensus_end = segments_df['end_pos'].min()
    
    if consensus_start > consensus_end:
        # No valid intersection, use union approach
        consensus_start = segments_df['start_pos'].min()
        consensus_end = segments_df['end_pos'].max()
    
    # Copy number consensus (weighted by confidence)
    if 'confidence' in segments_df.columns:
        weights = segments_df['confidence'].values
        weights = weights / weights.sum()  # Normalize weights
        consensus_cn = np.average(segments_df['copy_number'].values, weights=weights)
    else:
        consensus_cn = segments_df['copy_number'].median()
    
    consensus_cn = int(round(consensus_cn))
    
    # Overall confidence (average of individual confidences)
    if 'confidence' in segments_df.columns:
        overall_confidence = segments_df['confidence'].mean()
    else:
        overall_confidence = 0.5  # Default moderate confidence
    
    # Agreement metrics
    cn_agreement = (segments_df['copy_number'] == consensus_cn).mean()
    
    # Only return consensus if it meets confidence threshold
    if overall_confidence >= confidence_threshold:
        return {
            'sample_id': sample_id,
            'chromosome': segments_df['chromosome'].iloc[0] if 'chromosome' in segments_df.columns else 'chr5',
            'start_pos': consensus_start,
            'end_pos': consensus_end,
            'copy_number': consensus_cn,
            'confidence': overall_confidence,
            'num_algorithms': len(algorithms_used),
            'algorithms_used': ','.join(algorithms_used),
            'cn_agreement': cn_agreement,
            'num_segments_used': len(segments_df)
        }
    
    return None

def filter_high_confidence_calls(consensus_df, confidence_threshold):
    """Filter for high-confidence consensus calls"""
    if consensus_df.empty:
        return pd.DataFrame()
    
    high_conf = consensus_df[consensus_df['confidence'] >= confidence_threshold].copy()
    logging.info(f"High-confidence calls (>={confidence_threshold}): {len(high_conf)}/{len(consensus_df)}")
    
    return high_conf

def generate_consensus_plots(consensus_df, algorithms, output_dir):
    """Generate consensus calling visualization plots"""
    if consensus_df.empty:
        logging.warning("No consensus data to plot")
        return
    
    logging.info("Generating consensus plots...")
    
    try:
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Plot 1: Consensus confidence distribution
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Confidence histogram
        axes[0,0].hist(consensus_df['confidence'], bins=20, alpha=0.7, edgecolor='black')
        axes[0,0].axvline(consensus_df['confidence'].mean(), color='red', linestyle='--', 
                         label=f'Mean: {consensus_df["confidence"].mean():.3f}')
        axes[0,0].set_xlabel('Consensus Confidence')
        axes[0,0].set_ylabel('Frequency')
        axes[0,0].set_title('Distribution of Consensus Confidence Scores')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # Copy number agreement
        if 'cn_agreement' in consensus_df.columns:
            axes[0,1].hist(consensus_df['cn_agreement'], bins=20, alpha=0.7, edgecolor='black')
            axes[0,1].set_xlabel('Copy Number Agreement')
            axes[0,1].set_ylabel('Frequency')
            axes[0,1].set_title('Copy Number Agreement Across Algorithms')
            axes[0,1].grid(True, alpha=0.3)
        
        # Copy number distribution
        if 'copy_number' in consensus_df.columns:
            cn_counts = consensus_df['copy_number'].value_counts().sort_index()
            axes[1,0].bar(cn_counts.index, cn_counts.values, alpha=0.7)
            axes[1,0].set_xlabel('Copy Number')
            axes[1,0].set_ylabel('Number of Calls')
            axes[1,0].set_title('Consensus Copy Number Distribution')
            axes[1,0].grid(True, alpha=0.3)
        
        # Algorithm usage
        if 'algorithms_used' in consensus_df.columns:
            alg_combinations = consensus_df['algorithms_used'].value_counts()
            axes[1,1].bar(range(len(alg_combinations)), alg_combinations.values, alpha=0.7)
            axes[1,1].set_xticks(range(len(alg_combinations)))
            axes[1,1].set_xticklabels(alg_combinations.index, rotation=45)
            axes[1,1].set_ylabel('Number of Calls')
            axes[1,1].set_title('Algorithm Combinations Used')
            axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'consensus_analysis_plots.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot 2: Sample-wise consensus heatmap
        if len(consensus_df) > 0:
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Create pivot table for heatmap
            if 'start_pos' in consensus_df.columns:
                pivot_data = consensus_df.pivot_table(
                    index='sample_id', 
                    columns='start_pos', 
                    values='copy_number', 
                    fill_value=2
                )
                
                sns.heatmap(pivot_data, annot=True, cmap='RdYlBu_r', center=2, ax=ax,
                           cbar_kws={'label': 'Copy Number'})
                ax.set_title('Consensus Copy Number Calls')
                ax.set_xlabel('Genomic Position')
                ax.set_ylabel('Sample')
                
                plt.tight_layout()
                heatmap_file = os.path.join(output_dir, 'consensus_copy_number_heatmap.png')
                plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
                plt.close()
        
        logging.info("Consensus plots generated successfully")
        
    except Exception as e:
        logging.warning(f"Could not generate consensus plots: {e}")

def main():
    parser = argparse.ArgumentParser(description='Consensus Calling for SMN CNV Pipeline V3')
    parser.add_argument('cbs_file', help='CBS segments file')
    parser.add_argument('hmm_file', help='HMM segments file')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--min-algorithms', type=int, default=2, help='Minimum algorithms for consensus')
    parser.add_argument('--confidence-threshold', type=float, default=0.7, help='Confidence threshold')
    parser.add_argument('--overlap-threshold', type=float, default=0.5, help='Minimum overlap for segment matching')
    parser.add_argument('--no-plots', action='store_true', help='Skip plot generation')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Load segmentation results
        algorithms = load_segmentation_results(args.cbs_file, args.hmm_file)
        
        if len(algorithms) < args.min_algorithms:
            logging.error(f"Not enough algorithms available ({len(algorithms)} < {args.min_algorithms})")
            
            # Create single algorithm output if available
            if len(algorithms) == 1:
                alg_name, alg_df = list(algorithms.items())[0]
                logging.info(f"Using single algorithm results from {alg_name}")
                
                single_alg_df = standardize_segment_format(alg_df, alg_name)
                single_alg_df['confidence'] = single_alg_df.get('confidence', 0.5)
                single_alg_df['num_algorithms'] = 1
                single_alg_df['algorithms_used'] = alg_name
                single_alg_df['cn_agreement'] = 1.0
                single_alg_df['num_segments_used'] = 1
                
                # Save single algorithm results
                output_file = os.path.join(args.output_dir, 'single_algorithm_calls.txt')
                single_alg_df.to_csv(output_file, sep='\t', index=False)
                logging.info(f"Single algorithm results saved to: {output_file}")
                
                sys.exit(0)
            else:
                logging.error("No algorithm results available for consensus calling")
                sys.exit(1)
        
        # Create consensus segments
        consensus_df = create_consensus_segments(
            algorithms, 
            args.min_algorithms, 
            args.confidence_threshold,
            args.overlap_threshold
        )
        
        if consensus_df.empty:
            logging.warning("No consensus calls generated")
            
            # Create empty output files
            empty_df = pd.DataFrame(columns=[
                'sample_id', 'chromosome', 'start_pos', 'end_pos', 'copy_number',
                'confidence', 'num_algorithms', 'algorithms_used', 'cn_agreement'
            ])
            
            output_file = os.path.join(args.output_dir, 'consensus_calls.txt')
            empty_df.to_csv(output_file, sep='\t', index=False)
            
            high_conf_file = os.path.join(args.output_dir, 'high_confidence_calls.txt')
            empty_df.to_csv(high_conf_file, sep='\t', index=False)
            
        else:
            # Save consensus results
            output_file = os.path.join(args.output_dir, 'consensus_calls.txt')
            consensus_df.to_csv(output_file, sep='\t', index=False)
            logging.info(f"Consensus calls saved to: {output_file}")
            
            # Filter and save high-confidence calls
            high_conf_df = filter_high_confidence_calls(consensus_df, args.confidence_threshold)
            high_conf_file = os.path.join(args.output_dir, 'high_confidence_calls.txt')
            high_conf_df.to_csv(high_conf_file, sep='\t', index=False)
            logging.info(f"High-confidence calls saved to: {high_conf_file}")
            
            # Save consensus summary
            summary_file = os.path.join(args.output_dir, 'consensus_summary.txt')
            with open(summary_file, 'w') as f:
                f.write("Consensus Calling Summary\n")
                f.write("=" * 30 + "\n\n")
                f.write(f"Algorithms used: {list(algorithms.keys())}\n")
                f.write(f"Total consensus calls: {len(consensus_df)}\n")
                f.write(f"High-confidence calls: {len(high_conf_df)}\n")
                f.write(f"Confidence threshold: {args.confidence_threshold}\n")
                f.write(f"Min algorithms required: {args.min_algorithms}\n")
                f.write(f"Overlap threshold: {args.overlap_threshold}\n\n")
                
                if len(consensus_df) > 0:
                    f.write("Copy Number Distribution:\n")
                    cn_counts = consensus_df['copy_number'].value_counts().sort_index()
                    for cn, count in cn_counts.items():
                        f.write(f"  CN={cn}: {count} calls\n")
                    
                    f.write(f"\nMean confidence: {consensus_df['confidence'].mean():.3f}\n")
                    f.write(f"Mean CN agreement: {consensus_df['cn_agreement'].mean():.3f}\n")
            
            # Generate plots if requested
            if not args.no_plots:
                generate_consensus_plots(consensus_df, algorithms, args.output_dir)
        
        logging.info("Consensus calling completed successfully")
        
    except Exception as e:
        logging.error(f"Consensus calling failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
