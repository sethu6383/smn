#!/usr/bin/env python3

"""
consensus_calling.py - Multi-algorithm consensus CNV calling
Usage: python consensus_calling.py <cbs_segments> <hmm_segments> <output_dir> [--min-algorithms 2]
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
import warnings
warnings.filterwarnings('ignore')

class ConsensusCallerV3:
    """
    Advanced consensus CNV caller combining CBS and HMM results
    """
    
    def __init__(self, min_algorithms=2, confidence_threshold=0.7, 
                 overlap_threshold=0.5, max_cn_difference=1):
        """
        Initialize consensus caller parameters
        
        Args:
            min_algorithms: Minimum number of algorithms that must agree
            confidence_threshold: Minimum confidence for high-quality calls
            overlap_threshold: Minimum overlap for segment matching
            max_cn_difference: Maximum copy number difference for consensus
        """
        self.min_algorithms = min_algorithms
        self.confidence_threshold = confidence_threshold
        self.overlap_threshold = overlap_threshold
        self.max_cn_difference = max_cn_difference
    
    def _calculate_overlap(self, seg1, seg2):
        """Calculate overlap fraction between two segments"""
        start_max = max(seg1['start_exon'], seg2['start_exon'])
        end_min = min(seg1['end_exon'], seg2['end_exon'])
        
        if start_max > end_min:
            return 0.0
        
        overlap_length = end_min - start_max + 1
        seg1_length = seg1['end_exon'] - seg1['start_exon'] + 1
        seg2_length = seg2['end_exon'] - seg2['start_exon'] + 1
        
        # Return overlap as fraction of smaller segment
        min_length = min(seg1_length, seg2_length)
        return overlap_length / min_length if min_length > 0 else 0.0
    
    def _match_segments(self, cbs_segments, hmm_segments):
        """Match segments between CBS and HMM results"""
        matched_pairs = []
        
        for _, cbs_seg in cbs_segments.iterrows():
            best_match = None
            best_overlap = 0
            
            for _, hmm_seg in hmm_segments.iterrows():
                if cbs_seg['sample_id'] != hmm_seg['sample_id']:
                    continue
                
                overlap = self._calculate_overlap(cbs_seg, hmm_seg)
                
                if overlap >= self.overlap_threshold and overlap > best_overlap:
                    # Check if copy numbers are compatible
                    cn_diff = abs(cbs_seg['copy_number'] - hmm_seg['copy_number'])
                    if cn_diff <= self.max_cn_difference:
                        best_overlap = overlap
                        best_match = hmm_seg
            
            if best_match is not None:
                matched_pairs.append({
                    'cbs_segment': cbs_seg,
                    'hmm_segment': best_match,
                    'overlap': best_overlap
                })
        
        return matched_pairs
    
    def _calculate_consensus_confidence(self, cbs_seg, hmm_seg, overlap):
        """Calculate consensus confidence score"""
        
        # Base confidence from individual algorithms
        cbs_conf = 1.0 - cbs_seg.get('p_value', 0.5)  # Convert p-value to confidence
        hmm_conf = hmm_seg.get('confidence', 0.5)
        
        # Overlap bonus
        overlap_bonus = overlap * 0.2
        
        # Copy number agreement bonus
        cn_diff = abs(cbs_seg['copy_number'] - hmm_seg['copy_number'])
        agreement_bonus = (1.0 - cn_diff / 4.0) * 0.2  # Max bonus for exact match
        
        # Coverage consistency bonus
        cbs_coverage = cbs_seg.get('mean_value', 1.0)
        hmm_coverage = hmm_seg.get('mean_coverage', 1.0)
        
        if not (np.isnan(cbs_coverage) or np.isnan(hmm_coverage)) and hmm_coverage != 0:
            coverage_ratio = min(cbs_coverage, hmm_coverage) / max(cbs_coverage, hmm_coverage)
            coverage_bonus = coverage_ratio * 0.1
        else:
            coverage_bonus = 0.0
        
        # Combine confidences
        consensus_confidence = (
            0.4 * cbs_conf +           # CBS contribution
            0.4 * hmm_conf +           # HMM contribution
            overlap_bonus +            # Overlap bonus
            agreement_bonus +          # Agreement bonus
            coverage_bonus             # Coverage consistency
        )
        
        return min(1.0, consensus_confidence)
    
    def _resolve_copy_number(self, cbs_seg, hmm_seg, consensus_confidence):
        """Resolve consensus copy number from multiple algorithms"""
        
        cbs_cn = cbs_seg['copy_number']
        hmm_cn = hmm_seg['copy_number']
        
        # If exact agreement, use that
        if cbs_cn == hmm_cn:
            return cbs_cn
        
        # Weight by algorithm-specific confidence
        cbs_weight = 1.0 - cbs_seg.get('p_value', 0.5)
        hmm_weight = hmm_seg.get('confidence', 0.5)
        
        # Additional weighting based on segment characteristics
        # Favor HMM for segments with missing data
        missing_data_penalty = hmm_seg.get('num_missing', 0) / hmm_seg.get('num_markers', 1)
        if missing_data_penalty > 0.3:  # > 30% missing
            hmm_weight *= 1.2  # Boost HMM (designed for missing data)
        
        # Favor CBS for segments with high statistical significance
        if cbs_seg.get('p_value', 1.0) < 0.001:
            cbs_weight *= 1.2
        
        # Weighted average (rounded to nearest integer)
        total_weight = cbs_weight + hmm_weight
        if total_weight > 0:
            weighted_cn = (cbs_cn * cbs_weight + hmm_cn * hmm_weight) / total_weight
            consensus_cn = int(round(weighted_cn))
        else:
            # Default to more conservative estimate (closer to normal)
            consensus_cn = min(cbs_cn, hmm_cn, key=lambda x: abs(x - 2))
        
        return max(0, min(4, consensus_cn))  # Clamp to valid range
    
    def _create_consensus_segment(self, cbs_seg, hmm_seg, overlap):
        """Create consensus segment from matched CBS and HMM segments"""
        
        consensus_confidence = self._calculate_consensus_confidence(cbs_seg, hmm_seg, overlap)
        consensus_cn = self._resolve_copy_number(cbs_seg, hmm_seg, consensus_confidence)
        
        # Use the broader segment boundaries
        start_exon = min(cbs_seg['start_exon'], hmm_seg['start_exon'])
        end_exon = max(cbs_seg['end_exon'], hmm_seg['end_exon'])
        
        # Calculate consensus coverage
        cbs_coverage = cbs_seg.get('mean_value', np.nan)
        hmm_coverage = hmm_seg.get('mean_coverage', np.nan)
        
        if not np.isnan(cbs_coverage) and not np.isnan(hmm_coverage):
            consensus_coverage = (cbs_coverage + hmm_coverage) / 2
        elif not np.isnan(cbs_coverage):
            consensus_coverage = cbs_coverage
        elif not np.isnan(hmm_coverage):
            consensus_coverage = hmm_coverage
        else:
            consensus_coverage = np.nan
        
        # Determine quality classification
        if consensus_confidence >= 0.9:
            quality = 'high'
        elif consensus_confidence >= 0.7:
            quality = 'medium'
        else:
            quality = 'low'
        
        # Create consensus segment
        consensus_segment = {
            'sample_id': cbs_seg['sample_id'],
            'start_exon': start_exon,
            'end_exon': end_exon,
            'copy_number': consensus_cn,
            'consensus_confidence': consensus_confidence,
            'consensus_coverage': consensus_coverage,
            'quality': quality,
            'num_markers': end_exon - start_exon + 1,
            'algorithms_agree': 2,  # CBS + HMM
            'cbs_copy_number': cbs_seg['copy_number'],
            'hmm_copy_number': hmm_seg['copy_number'],
            'cbs_pvalue': cbs_seg.get('p_value', np.nan),
            'hmm_confidence': hmm_seg.get('confidence', np.nan),
            'overlap_fraction': overlap,
            'has_missing_data': hmm_seg.get('num_missing', 0) > 0
        }
        
        return consensus_segment
    
    def _handle_unmatched_segments(self, segments, algorithm_name):
        """Handle segments that don't have matches from other algorithms"""
        unmatched_consensus = []
        
        for _, seg in segments.iterrows():
            # Lower confidence for single-algorithm calls
            if algorithm_name == 'CBS':
                confidence = (1.0 - seg.get('p_value', 0.5)) * 0.7  # Reduce confidence
                coverage = seg.get('mean_value', np.nan)
            else:  # HMM
                confidence = seg.get('confidence', 0.5) * 0.7
                coverage = seg.get('mean_coverage', np.nan)
            
            # Only include if meets minimum confidence
            if confidence >= self.confidence_threshold * 0.8:  # Slightly lower threshold
                quality = 'medium' if confidence >= 0.6 else 'low'
                
                consensus_segment = {
                    'sample_id': seg['sample_id'],
                    'start_exon': seg['start_exon'],
                    'end_exon': seg['end_exon'],
                    'copy_number': seg['copy_number'],
                    'consensus_confidence': confidence,
                    'consensus_coverage': coverage,
                    'quality': quality,
                    'num_markers': seg.get('num_markers', seg['end_exon'] - seg['start_exon'] + 1),
                    'algorithms_agree': 1,  # Single algorithm
                    'cbs_copy_number': seg['copy_number'] if algorithm_name == 'CBS' else np.nan,
                    'hmm_copy_number': seg['copy_number'] if algorithm_name == 'HMM' else np.nan,
                    'cbs_pvalue': seg.get('p_value', np.nan) if algorithm_name == 'CBS' else np.nan,
                    'hmm_confidence': seg.get('confidence', np.nan) if algorithm_name == 'HMM' else np.nan,
                    'overlap_fraction': np.nan,
                    'has_missing_data': seg.get('num_missing', 0) > 0 if algorithm_name == 'HMM' else False,
                    'single_algorithm': algorithm_name
                }
                
                unmatched_consensus.append(consensus_segment)
        
        return unmatched_consensus
    
    def call_consensus(self, cbs_segments, hmm_segments):
        """Generate consensus CNV calls from CBS and HMM results"""
        print("Generating consensus CNV calls...")
        
        consensus_calls = []
        
        # Match segments between algorithms
        matched_pairs = self._match_segments(cbs_segments, hmm_segments)
        print(f"Found {len(matched_pairs)} matched segment pairs")
        
        # Create consensus for matched segments
        for pair in matched_pairs:
            consensus_seg = self._create_consensus_segment(
                pair['cbs_segment'], pair['hmm_segment'], pair['overlap']
            )
            consensus_calls.append(consensus_seg)
        
        # Track which segments were matched
        matched_cbs_indices = set()
        matched_hmm_indices = set()
        
        for pair in matched_pairs:
            matched_cbs_indices.add(pair['cbs_segment'].name)
            matched_hmm_indices.add(pair['hmm_segment'].name)
        
        # Handle unmatched CBS segments
        unmatched_cbs = cbs_segments[~cbs_segments.index.isin(matched_cbs_indices)]
        unmatched_cbs_consensus = self._handle_unmatched_segments(unmatched_cbs, 'CBS')
        consensus_calls.extend(unmatched_cbs_consensus)
        
        # Handle unmatched HMM segments
        unmatched_hmm = hmm_segments[~hmm_segments.index.isin(matched_hmm_indices)]
        unmatched_hmm_consensus = self._handle_unmatched_segments(unmatched_hmm, 'HMM')
        consensus_calls.extend(unmatched_hmm_consensus)
        
        # Convert to DataFrame
        consensus_df = pd.DataFrame(consensus_calls)
        
        # Filter by minimum algorithms requirement
        if self.min_algorithms > 1:
            high_confidence_calls = consensus_df[consensus_df['algorithms_agree'] >= self.min_algorithms]
            print(f"High-confidence consensus calls (≥{self.min_algorithms} algorithms): {len(high_confidence_calls)}")
        else:
            high_confidence_calls = consensus_df
        
        return consensus_df, high_confidence_calls
    
    def generate_quality_metrics(self, consensus_df):
        """Generate quality metrics for consensus calls"""
        metrics = {}
        
        # Overall statistics
        metrics['total_calls'] = len(consensus_df)
        metrics['high_quality_calls'] = len(consensus_df[consensus_df['quality'] == 'high'])
        metrics['medium_quality_calls'] = len(consensus_df[consensus_df['quality'] == 'medium'])
        metrics['low_quality_calls'] = len(consensus_df[consensus_df['quality'] == 'low'])
        
        # Algorithm agreement
        metrics['multi_algorithm_calls'] = len(consensus_df[consensus_df['algorithms_agree'] >= 2])
        metrics['single_algorithm_calls'] = len(consensus_df[consensus_df['algorithms_agree'] == 1])
        
        # Confidence distribution
        if len(consensus_df) > 0:
            metrics['mean_confidence'] = consensus_df['consensus_confidence'].mean()
            metrics['median_confidence'] = consensus_df['consensus_confidence'].median()
            metrics['confidence_std'] = consensus_df['consensus_confidence'].std()
        
        # Copy number distribution
        cn_distribution = consensus_df['copy_number'].value_counts().sort_index().to_dict()
        metrics['copy_number_distribution'] = cn_distribution
        
        # Missing data impact
        if 'has_missing_data' in consensus_df.columns:
            missing_data_calls = consensus_df[consensus_df['has_missing_data']]
            metrics['calls_with_missing_data'] = len(missing_data_calls)
            if len(missing_data_calls) > 0:
                metrics['missing_data_confidence'] = missing_data_calls['consensus_confidence'].mean()
        
        return metrics

def create_consensus_plots(consensus_df, cbs_df, hmm_df, output_dir):
    """Create visualization plots for consensus calling results"""
    
    plot_dir = Path(output_dir) / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    # Plot 1: Consensus overview for selected samples
    samples = consensus_df['sample_id'].unique()
    n_samples = min(3, len(samples))
    
    if n_samples > 0:
        fig, axes = plt.subplots(n_samples, 1, figsize=(15, 4*n_samples))
        if n_samples == 1:
            axes = [axes]
        
        for i, sample_id in enumerate(samples[:n_samples]):
            ax = axes[i]
            
            # Get segments for this sample
            sample_consensus = consensus_df[consensus_df['sample_id'] == sample_id]
            sample_cbs = cbs_df[cbs_df['sample_id'] == sample_id] if 'sample_id' in cbs_df.columns else pd.DataFrame()
            sample_hmm = hmm_df[hmm_df['sample_id'] == sample_id] if 'sample_id' in hmm_df.columns else pd.DataFrame()
            
            y_positions = {'CBS': 0.1, 'HMM': 0.2, 'Consensus': 0.3}
            colors = ['red', 'orange', 'green', 'blue', 'purple']
            
            # Plot CBS segments
            for _, seg in sample_cbs.iterrows():
                color = colors[min(seg.get('copy_number', 2), 4)]
                ax.barh(y_positions['CBS'], seg['end_exon'] - seg['start_exon'] + 1, 
                       left=seg['start_exon'], height=0.05, color=color, alpha=0.6)
                
                # Add copy number label
                mid_pos = (seg['start_exon'] + seg['end_exon']) / 2
                ax.text(mid_pos, y_positions['CBS'], f"CN={seg.get('copy_number', '?')}", 
                       ha='center', va='center', fontsize=8)
            
            # Plot HMM segments
            for _, seg in sample_hmm.iterrows():
                color = colors[min(seg.get('copy_number', 2), 4)]
                ax.barh(y_positions['HMM'], seg['end_exon'] - seg['start_exon'] + 1, 
                       left=seg['start_exon'], height=0.05, color=color, alpha=0.6)
                
                mid_pos = (seg['start_exon'] + seg['end_exon']) / 2
                ax.text(mid_pos, y_positions['HMM'], f"CN={seg.get('copy_number', '?')}", 
                       ha='center', va='center', fontsize=8)
            
            # Plot consensus segments
            for _, seg in sample_consensus.iterrows():
                color = colors[min(seg['copy_number'], 4)]
                # Width proportional to confidence
                alpha = 0.4 + 0.6 * seg['consensus_confidence']
                
                ax.barh(y_positions['Consensus'], seg['end_exon'] - seg['start_exon'] + 1, 
                       left=seg['start_exon'], height=0.05, color=color, alpha=alpha)
                
                mid_pos = (seg['start_exon'] + seg['end_exon']) / 2
                label = f"CN={seg['copy_number']}\n({seg['consensus_confidence']:.2f})"
                ax.text(mid_pos, y_positions['Consensus'], label, 
                       ha='center', va='center', fontsize=8, weight='bold')
            
            ax.set_yticks(list(y_positions.values()))
            ax.set_yticklabels(list(y_positions.keys()))
            ax.set_xlabel('Exon Position')
            ax.set_title(f'Consensus Calling Overview - {sample_id}')
            ax.grid(True, axis='x', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plot_dir / 'consensus_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Plot 2: Quality and confidence metrics
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Confidence distribution by quality
    if 'quality' in consensus_df.columns:
        for quality in ['high', 'medium', 'low']:
            quality_data = consensus_df[consensus_df['quality'] == quality]['consensus_confidence']
            if len(quality_data) > 0:
                axes[0, 0].hist(quality_data, alpha=0.7, label=f'{quality} quality', bins=15)
        axes[0, 0].set_xlabel('Consensus Confidence')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Confidence Distribution by Quality')
        axes[0, 0].legend()
    
    # Algorithm agreement
    if 'algorithms_agree' in consensus_df.columns:
        agreement_counts = consensus_df['algorithms_agree'].value_counts().sort_index()
        axes[0, 1].bar(agreement_counts.index, agreement_counts.values, alpha=0.7)
        axes[0, 1].set_xlabel('Number of Algorithms Agreeing')
        axes[0, 1].set_ylabel('Number of Calls')
        axes[0, 1].set_title('Algorithm Agreement Distribution')
    
    # Copy number comparison
    if 'cbs_copy_number' in consensus_df.columns and 'hmm_copy_number' in consensus_df.columns:
        # Only plot where both algorithms have calls
        both_calls = consensus_df.dropna(subset=['cbs_copy_number', 'hmm_copy_number'])
        if len(both_calls) > 0:
            scatter = axes[1, 0].scatter(both_calls['cbs_copy_number'], both_calls['hmm_copy_number'], 
                                       c=both_calls['consensus_confidence'], cmap='viridis', alpha=0.7)
            axes[1, 0].plot([0, 4], [0, 4], 'r--', alpha=0.5)  # Perfect agreement line
            axes[1, 0].set_xlabel('CBS Copy Number')
            axes[1, 0].set_ylabel('HMM Copy Number')
            axes[1, 0].set_title('Algorithm Copy Number Comparison')
            plt.colorbar(scatter, ax=axes[1, 0], label='Consensus Confidence')
    
    # Segment length vs confidence
    axes[1, 1].scatter(consensus_df['num_markers'], consensus_df['consensus_confidence'], alpha=0.6)
    axes[1, 1].set_xlabel('Segment Length (markers)')
    axes[1, 1].set_ylabel('Consensus Confidence')
    axes[1, 1].set_title('Segment Length vs Confidence')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'consensus_quality_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Copy number distribution comparison
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # CBS copy number distribution
    if not cbs_df.empty and 'copy_number' in cbs_df.columns:
        cbs_cn_counts = cbs_df['copy_number'].value_counts().sort_index()
        axes[0].bar(cbs_cn_counts.index, cbs_cn_counts.values, alpha=0.7, color='blue')
        axes[0].set_title('CBS Copy Number Distribution')
        axes[0].set_xlabel('Copy Number')
        axes[0].set_ylabel('Number of Segments')
    
    # HMM copy number distribution
    if not hmm_df.empty and 'copy_number' in hmm_df.columns:
        hmm_cn_counts = hmm_df['copy_number'].value_counts().sort_index()
        axes[1].bar(hmm_cn_counts.index, hmm_cn_counts.values, alpha=0.7, color='orange')
        axes[1].set_title('HMM Copy Number Distribution')
        axes[1].set_xlabel('Copy Number')
        axes[1].set_ylabel('Number of Segments')
    
    # Consensus copy number distribution
    consensus_cn_counts = consensus_df['copy_number'].value_counts().sort_index()
    bars = axes[2].bar(consensus_cn_counts.index, consensus_cn_counts.values, alpha=0.7, color='green')
    axes[2].set_title('Consensus Copy Number Distribution')
    axes[2].set_xlabel('Copy Number')
    axes[2].set_ylabel('Number of Calls')
    
    # Add count labels
    for bar, count in zip(bars, consensus_cn_counts.values):
        axes[2].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                    f'{count}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'algorithm_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Consensus calling plots saved to: {plot_dir}")

def main():
    parser = argparse.ArgumentParser(description='Multi-algorithm consensus CNV calling')
    parser.add_argument('cbs_segments', help='CBS segmentation results file')
    parser.add_argument('hmm_segments', help='HMM segmentation results file')
    parser.add_argument('output_dir', help='Output directory for consensus results')
    parser.add_argument('--min-algorithms', type=int, default=2, 
                       help='Minimum algorithms required for consensus')
    parser.add_argument('--confidence-threshold', type=float, default=0.7, 
                       help='Minimum confidence threshold')
    parser.add_argument('--overlap-threshold', type=float, default=0.5, 
                       help='Minimum overlap for segment matching')
    parser.add_argument('--max-cn-difference', type=int, default=1, 
                       help='Maximum copy number difference for consensus')
    parser.add_argument('--no-plots', action='store_true', help='Skip creating plots')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load segmentation results
    print("Loading segmentation results...")
    try:
        cbs_df = pd.read_csv(args.cbs_segments, sep='\t')
        print(f"Loaded CBS results: {len(cbs_df)} segments")
    except Exception as e:
        print(f"Error loading CBS file: {e}")
        sys.exit(1)
    
    try:
        hmm_df = pd.read_csv(args.hmm_segments, sep='\t')
        print(f"Loaded HMM results: {len(hmm_df)} segments")
    except Exception as e:
        print(f"Error loading HMM file: {e}")
        sys.exit(1)
    
    # Initialize consensus caller
    consensus_caller = ConsensusCallerV3(
        min_algorithms=args.min_algorithms,
        confidence_threshold=args.confidence_threshold,
        overlap_threshold=args.overlap_threshold,
        max_cn_difference=args.max_cn_difference
    )
    
    # Generate consensus calls
    consensus_df, high_confidence_df = consensus_caller.call_consensus(cbs_df, hmm_df)
    
    # Generate quality metrics
    quality_metrics = consensus_caller.generate_quality_metrics(consensus_df)
    
    # Save results
    output_file = output_dir / 'consensus_calls.txt'
    consensus_df.to_csv(output_file, sep='\t', index=False)
    
    high_conf_file = output_dir / 'high_confidence_calls.txt'
    high_confidence_df.to_csv(high_conf_file, sep='\t', index=False)
    
    # Save quality metrics
    metrics_file = output_dir / 'quality_metrics.txt'
    with open(metrics_file, 'w') as f:
        f.write("# Consensus Calling Quality Metrics\n")
        for key, value in quality_metrics.items():
            f.write(f"{key}\t{value}\n")
    
    # Save parameters
    params_file = output_dir / 'consensus_parameters.txt'
    with open(params_file, 'w') as f:
        f.write("# Consensus Calling Parameters\n")
        f.write(f"min_algorithms\t{args.min_algorithms}\n")
        f.write(f"confidence_threshold\t{args.confidence_threshold}\n")
        f.write(f"overlap_threshold\t{args.overlap_threshold}\n")
        f.write(f"max_cn_difference\t{args.max_cn_difference}\n")
    
    # Create plots
    if not args.no_plots and len(consensus_df) > 0:
        try:
            create_consensus_plots(consensus_df, cbs_df, hmm_df, output_dir)
        except Exception as e:
            print(f"Warning: Could not create plots: {e}")
    
    # Print summary
    print(f"\nConsensus Calling Results:")
    print(f"  Total consensus calls: {len(consensus_df)}")
    print(f"  High-confidence calls: {len(high_confidence_df)}")
    print(f"  Multi-algorithm agreement: {quality_metrics.get('multi_algorithm_calls', 0)}")
    print(f"  Single-algorithm calls: {quality_metrics.get('single_algorithm_calls', 0)}")
    
    if len(consensus_df) > 0:
        print(f"\nQuality distribution:")
        print(f"  High quality: {quality_metrics.get('high_quality_calls', 0)}")
        print(f"  Medium quality: {quality_metrics.get('medium_quality_calls', 0)}")
        print(f"  Low quality: {quality_metrics.get('low_quality_calls', 0)}")
        
        print(f"\nConfidence statistics:")
        print(f"  Mean confidence: {quality_metrics.get('mean_confidence', 0):.3f}")
        print(f"  Median confidence: {quality_metrics.get('median_confidence', 0):.3f}")
        
        print(f"\nCopy number distribution:")
        cn_dist = quality_metrics.get('copy_number_distribution', {})
        for cn, count in sorted(cn_dist.items()):
            print(f"  CN={cn}: {count} calls")
    
    print(f"\nResults saved to:")
    print(f"  All consensus calls: {output_file}")
    print(f"  High-confidence calls: {high_conf_file}")
    print(f"  Quality metrics: {metrics_file}")
    print(f"  Parameters: {params_file}")

if __name__ == "__main__":
    main()
