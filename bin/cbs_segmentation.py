#!/usr/bin/env python3

"""
cbs_segmentation.py - Circular Binary Segmentation for CNV detection
Usage: python cbs_segmentation.py <normalized_coverage_file> <output_dir> [--alpha 0.01]
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
from collections import namedtuple
import warnings
warnings.filterwarnings('ignore')

# Define segment structure
Segment = namedtuple('Segment', ['sample_id', 'start_pos', 'end_pos', 'start_exon', 'end_exon', 
                                'mean_value', 'num_markers', 'p_value', 'copy_number'])

class CircularBinarySegmentation:
    """
    Circular Binary Segmentation algorithm for detecting copy number variations
    Based on Olshen et al. (2004) Biostatistics
    """
    
    def __init__(self, alpha=0.01, min_segment_markers=3, nperm=10000, trim=0.025, 
                 undo_SD=2.0, undo_splits=4):
        """
        Initialize CBS parameters
        
        Args:
            alpha: Significance level for change-point detection
            min_segment_markers: Minimum number of markers per segment
            nperm: Number of permutations for p-value calculation
            trim: Trimming fraction for outlier removal
            undo_SD: Standard deviation threshold for undoing splits
            undo_splits: Maximum number of splits to undo
        """
        self.alpha = alpha
        self.min_segment_markers = min_segment_markers
        self.nperm = nperm
        self.trim = trim
        self.undo_SD = undo_SD
        self.undo_splits = undo_splits
    
    def _calculate_test_statistic(self, x, i, j, k):
        """
        Calculate CBS test statistic for change-point at position k between i and j
        """
        if k <= i or k >= j or j - i < 2:
            return 0.0
        
        # Left segment: i to k-1
        left_data = x[i:k]
        # Right segment: k to j
        right_data = x[k:j+1]
        
        if len(left_data) == 0 or len(right_data) == 0:
            return 0.0
        
        n1, n2 = len(left_data), len(right_data)
        mean1, mean2 = np.mean(left_data), np.mean(right_data)
        
        # Calculate test statistic
        if n1 > 0 and n2 > 0:
            # Use t-test statistic weighted by segment sizes
            pooled_var = (np.sum((left_data - mean1)**2) + np.sum((right_data - mean2)**2)) / (n1 + n2 - 2)
            if pooled_var > 0:
                t_stat = (mean1 - mean2) / np.sqrt(pooled_var * (1/n1 + 1/n2))
                # Convert to CBS statistic
                cbs_stat = abs(t_stat) * np.sqrt(n1 * n2 / (n1 + n2))
                return cbs_stat
        
        return 0.0
    
    def _find_changepoint(self, x, i, j):
        """
        Find the optimal change-point between positions i and j
        """
        if j - i < 2 * self.min_segment_markers:
            return None, 0.0
        
        max_stat = 0.0
        best_k = None
        
        # Search for optimal change-point
        for k in range(i + self.min_segment_markers, j - self.min_segment_markers + 1):
            stat = self._calculate_test_statistic(x, i, j, k)
            if stat > max_stat:
                max_stat = stat
                best_k = k
        
        return best_k, max_stat
    
    def _permutation_test(self, x, i, j, observed_stat):
        """
        Calculate p-value using permutation test
        """
        if observed_stat == 0 or self.nperm == 0:
            return 1.0
        
        segment_data = x[i:j+1].copy()
        n = len(segment_data)
        
        if n < 2 * self.min_segment_markers:
            return 1.0
        
        # Perform permutations
        perm_stats = []
        for _ in range(self.nperm):
            # Shuffle the data
            perm_data = np.random.permutation(segment_data)
            
            # Find best change-point in permuted data
            max_perm_stat = 0.0
            for k in range(self.min_segment_markers, n - self.min_segment_markers):
                stat = self._calculate_test_statistic(perm_data, 0, n-1, k)
                max_perm_stat = max(max_perm_stat, stat)
            
            perm_stats.append(max_perm_stat)
        
        # Calculate p-value
        p_value = np.sum(np.array(perm_stats) >= observed_stat) / self.nperm
        return p_value
    
    def _segment_recursive(self, x, positions, sample_id, i=0, j=None):
        """
        Recursive segmentation algorithm
        """
        if j is None:
            j = len(x) - 1
        
        if j - i < 2 * self.min_segment_markers:
            # Create single segment
            mean_val = np.mean(x[i:j+1])
            return [Segment(
                sample_id=sample_id,
                start_pos=positions[i] if positions is not None else i,
                end_pos=positions[j] if positions is not None else j,
                start_exon=i,
                end_exon=j,
                mean_value=mean_val,
                num_markers=j - i + 1,
                p_value=1.0,
                copy_number=self._estimate_copy_number(mean_val)
            )]
        
        # Find change-point
        k, test_stat = self._find_changepoint(x, i, j)
        
        if k is None:
            # No significant change-point found
            mean_val = np.mean(x[i:j+1])
            return [Segment(
                sample_id=sample_id,
                start_pos=positions[i] if positions is not None else i,
                end_pos=positions[j] if positions is not None else j,
                start_exon=i,
                end_exon=j,
                mean_value=mean_val,
                num_markers=j - i + 1,
                p_value=1.0,
                copy_number=self._estimate_copy_number(mean_val)
            )]
        
        # Calculate p-value for the change-point
        p_value = self._permutation_test(x, i, j, test_stat)
        
        if p_value > self.alpha:
            # Change-point not significant
            mean_val = np.mean(x[i:j+1])
            return [Segment(
                sample_id=sample_id,
                start_pos=positions[i] if positions is not None else i,
                end_pos=positions[j] if positions is not None else j,
                start_exon=i,
                end_exon=j,
                mean_value=mean_val,
                num_markers=j - i + 1,
                p_value=p_value,
                copy_number=self._estimate_copy_number(mean_val)
            )]
        
        # Recursively segment left and right parts
        left_segments = self._segment_recursive(x, positions, sample_id, i, k-1)
        right_segments = self._segment_recursive(x, positions, sample_id, k, j)
        
        return left_segments + right_segments
    
    def _estimate_copy_number(self, normalized_value):
        """
        Convert normalized coverage value to copy number estimate
        Assumes normalized value of 1.0 corresponds to CN=2
        """
        if normalized_value <= 0.25:
            return 0
        elif normalized_value <= 0.75:
            return 1
        elif normalized_value <= 1.5:
            return 2
        elif normalized_value <= 2.25:
            return 3
        else:
            return 4
    
    def _undo_splits(self, segments):
        """
        Undo splits that don't meet the undo criteria
        """
        if len(segments) <= 1 or self.undo_splits == 0:
            return segments
        
        # Calculate standard deviation across all segments
        all_values = []
        for seg in segments:
            all_values.extend([seg.mean_value] * seg.num_markers)
        
        global_sd = np.std(all_values) if len(all_values) > 1 else 0
        
        if global_sd == 0:
            return segments
        
        # Merge adjacent segments if difference is small
        merged_segments = []
        current_seg = segments[0]
        
        for next_seg in segments[1:]:
            # Calculate difference in standard deviation units
            diff = abs(current_seg.mean_value - next_seg.mean_value) / global_sd
            
            if diff < self.undo_SD:
                # Merge segments
                total_markers = current_seg.num_markers + next_seg.num_markers
                weighted_mean = (current_seg.mean_value * current_seg.num_markers + 
                               next_seg.mean_value * next_seg.num_markers) / total_markers
                
                current_seg = Segment(
                    sample_id=current_seg.sample_id,
                    start_pos=current_seg.start_pos,
                    end_pos=next_seg.end_pos,
                    start_exon=current_seg.start_exon,
                    end_exon=next_seg.end_exon,
                    mean_value=weighted_mean,
                    num_markers=total_markers,
                    p_value=max(current_seg.p_value, next_seg.p_value),
                    copy_number=self._estimate_copy_number(weighted_mean)
                )
            else:
                merged_segments.append(current_seg)
                current_seg = next_seg
        
        merged_segments.append(current_seg)
        
        return merged_segments
    
    def segment(self, coverage_data, exon_positions=None):
        """
        Perform CBS segmentation on coverage data
        
        Args:
            coverage_data: DataFrame with samples as rows, exons as columns
            exon_positions: Optional list of genomic positions for exons
            
        Returns:
            List of Segment objects
        """
        all_segments = []
        
        for sample_id in coverage_data.index:
            sample_data = coverage_data.loc[sample_id].values
            
            # Remove NaN values and track positions
            valid_indices = ~np.isnan(sample_data)
            valid_data = sample_data[valid_indices]
            
            if exon_positions is not None:
                valid_positions = np.array(exon_positions)[valid_indices]
            else:
                valid_positions = np.arange(len(valid_data))
            
            if len(valid_data) < 2 * self.min_segment_markers:
                # Not enough data for segmentation
                if len(valid_data) > 0:
                    mean_val = np.mean(valid_data)
                    segment = Segment(
                        sample_id=sample_id,
                        start_pos=valid_positions[0] if len(valid_positions) > 0 else 0,
                        end_pos=valid_positions[-1] if len(valid_positions) > 0 else 0,
                        start_exon=0,
                        end_exon=len(valid_data) - 1,
                        mean_value=mean_val,
                        num_markers=len(valid_data),
                        p_value=1.0,
                        copy_number=self._estimate_copy_number(mean_val)
                    )
                    all_segments.append(segment)
                continue
            
            # Perform segmentation
            segments = self._segment_recursive(valid_data, valid_positions, sample_id)
            
            # Apply undo splits
            segments = self._undo_splits(segments)
            
            all_segments.extend(segments)
        
        return all_segments

def load_coverage_data(coverage_file):
    """Load normalized coverage data from file"""
    try:
        if coverage_file.endswith(('.txt', '.tsv')):
            df = pd.read_csv(coverage_file, sep='\t', index_col=0)
        else:
            df = pd.read_csv(coverage_file, index_col=0)
        
        print(f"Loaded coverage data: {df.shape[0]} samples, {df.shape[1]} exons")
        return df
    except Exception as e:
        print(f"Error loading coverage file: {e}")
        sys.exit(1)

def create_segmentation_plots(segments_df, coverage_data, output_dir):
    """Create visualization plots for CBS segmentation results"""
    
    plot_dir = Path(output_dir) / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    # Plot 1: Segmentation overview
    samples = segments_df['sample_id'].unique()
    n_samples = len(samples)
    
    fig, axes = plt.subplots(min(4, n_samples), 1, figsize=(15, 3*min(4, n_samples)))
    if n_samples == 1:
        axes = [axes]
    
    for i, sample_id in enumerate(samples[:4]):  # Show first 4 samples
        ax = axes[i] if i < len(axes) else axes[-1]
        
        # Get sample data
        sample_segments = segments_df[segments_df['sample_id'] == sample_id]
        sample_coverage = coverage_data.loc[sample_id] if sample_id in coverage_data.index else None
        
        if sample_coverage is not None:
            # Plot original coverage
            exon_positions = range(len(sample_coverage))
            ax.scatter(exon_positions, sample_coverage, alpha=0.6, s=30, color='gray', label='Coverage')
            
            # Plot segments
            for _, segment in sample_segments.iterrows():
                start_idx = segment['start_exon']
                end_idx = segment['end_exon']
                mean_val = segment['mean_value']
                
                # Draw segment line
                ax.plot([start_idx, end_idx], [mean_val, mean_val], 
                       color='red', linewidth=3, alpha=0.8)
                
                # Annotate copy number
                mid_pos = (start_idx + end_idx) / 2
                ax.annotate(f"CN={segment['copy_number']}", 
                           xy=(mid_pos, mean_val), xytext=(mid_pos, mean_val + 0.1),
                           ha='center', fontsize=8, color='red')
        
        ax.set_title(f'CBS Segmentation - {sample_id}')
        ax.set_xlabel('Exon Position')
        ax.set_ylabel('Normalized Coverage')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'cbs_segmentation_overview.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Copy number distribution
    plt.figure(figsize=(10, 6))
    
    cn_counts = segments_df['copy_number'].value_counts().sort_index()
    colors = ['red', 'orange', 'green', 'blue', 'purple']
    
    bars = plt.bar(cn_counts.index, cn_counts.values, 
                   color=[colors[min(cn, 4)] for cn in cn_counts.index], alpha=0.7)
    
    plt.title('Copy Number Distribution from CBS Segmentation')
    plt.xlabel('Copy Number')
    plt.ylabel('Number of Segments')
    
    # Add count labels on bars
    for bar, count in zip(bars, cn_counts.values):
        plt.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                f'{count}', ha='center', va='bottom')
    
    plt.savefig(plot_dir / 'cbs_copy_number_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Segment statistics
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Segment length distribution
    axes[0, 0].hist(segments_df['num_markers'], bins=20, alpha=0.7, edgecolor='black')
    axes[0, 0].set_xlabel('Segment Length (markers)')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Segment Length Distribution')
    
    # P-value distribution
    axes[0, 1].hist(segments_df['p_value'], bins=20, alpha=0.7, edgecolor='black')
    axes[0, 1].axvline(x=0.05, color='red', linestyle='--', label='α = 0.05')
    axes[0, 1].set_xlabel('P-value')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Changepoint P-value Distribution')
    axes[0, 1].legend()
    
    # Mean value vs copy number
    axes[1, 0].scatter(segments_df['copy_number'], segments_df['mean_value'], alpha=0.6)
    axes[1, 0].set_xlabel('Copy Number')
    axes[1, 0].set_ylabel('Mean Normalized Coverage')
    axes[1, 0].set_title('Copy Number vs Coverage')
    
    # Segments per sample
    segments_per_sample = segments_df['sample_id'].value_counts()
    axes[1, 1].hist(segments_per_sample, bins=20, alpha=0.7, edgecolor='black')
    axes[1, 1].set_xlabel('Number of Segments')
    axes[1, 1].set_ylabel('Number of Samples')
    axes[1, 1].set_title('Segments per Sample Distribution')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'cbs_segment_statistics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"CBS segmentation plots saved to: {plot_dir}")

def main():
    parser = argparse.ArgumentParser(description='Circular Binary Segmentation for CNV detection')
    parser.add_argument('coverage_file', help='Normalized coverage data file')
    parser.add_argument('output_dir', help='Output directory for CBS results')
    parser.add_argument('--alpha', type=float, default=0.01, 
                       help='Significance level for change-point detection')
    parser.add_argument('--min-markers', type=int, default=3, 
                       help='Minimum number of markers per segment')
    parser.add_argument('--nperm', type=int, default=10000, 
                       help='Number of permutations for p-value calculation')
    parser.add_argument('--trim', type=float, default=0.025, 
                       help='Trimming fraction for outlier removal')
    parser.add_argument('--undo-SD', type=float, default=2.0, 
                       help='Standard deviation threshold for undoing splits')
    parser.add_argument('--undo-splits', type=int, default=4, 
                       help='Maximum number of splits to undo')
    parser.add_argument('--no-plots', action='store_true', help='Skip creating plots')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load coverage data
    print("Loading normalized coverage data...")
    coverage_data = load_coverage_data(args.coverage_file)
    
    # Initialize CBS
    print("Initializing Circular Binary Segmentation...")
    cbs = CircularBinarySegmentation(
        alpha=args.alpha,
        min_segment_markers=args.min_markers,
        nperm=args.nperm,
        trim=args.trim,
        undo_SD=args.undo_SD,
        undo_splits=args.undo_splits
    )
    
    # Perform segmentation
    print("Performing CBS segmentation...")
    segments = cbs.segment(coverage_data)
    
    # Convert to DataFrame
    segments_data = []
    for segment in segments:
        segments_data.append({
            'sample_id': segment.sample_id,
            'start_pos': segment.start_pos,
            'end_pos': segment.end_pos,
            'start_exon': segment.start_exon,
            'end_exon': segment.end_exon,
            'mean_value': segment.mean_value,
            'num_markers': segment.num_markers,
            'p_value': segment.p_value,
            'copy_number': segment.copy_number
        })
    
    segments_df = pd.DataFrame(segments_data)
    
    # Save results
    output_file = output_dir / 'cbs_segments.txt'
    segments_df.to_csv(output_file, sep='\t', index=False)
    
    # Save parameters
    params_file = output_dir / 'cbs_parameters.txt'
    with open(params_file, 'w') as f:
        f.write("# CBS Segmentation Parameters\n")
        f.write(f"alpha\t{args.alpha}\n")
        f.write(f"min_markers\t{args.min_markers}\n")
        f.write(f"nperm\t{args.nperm}\n")
        f.write(f"trim\t{args.trim}\n")
        f.write(f"undo_SD\t{args.undo_SD}\n")
        f.write(f"undo_splits\t{args.undo_splits}\n")
    
    # Create plots
    if not args.no_plots and len(segments_df) > 0:
        try:
            create_segmentation_plots(segments_df, coverage_data, output_dir)
        except Exception as e:
            print(f"Warning: Could not create plots: {e}")
    
    # Print summary
    print(f"\nCBS Segmentation Results:")
    print(f"  Total segments detected: {len(segments_df)}")
    print(f"  Samples analyzed: {segments_df['sample_id'].nunique()}")
    print(f"  Average segments per sample: {len(segments_df) / segments_df['sample_id'].nunique():.1f}")
    
    if len(segments_df) > 0:
        print(f"\nCopy number distribution:")
        cn_counts = segments_df['copy_number'].value_counts().sort_index()
        for cn, count in cn_counts.items():
            print(f"  CN={cn}: {count} segments")
        
        print(f"\nSegment statistics:")
        print(f"  Mean segment length: {segments_df['num_markers'].mean():.1f} markers")
        print(f"  Median p-value: {segments_df['p_value'].median():.4f}")
    
    print(f"\nResults saved to:")
    print(f"  Segments: {output_file}")
    print(f"  Parameters: {params_file}")

if __name__ == "__main__":
    main()
