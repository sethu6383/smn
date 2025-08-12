#!/usr/bin/env python3
"""
SMN CNV Detection Pipeline V3 - Depth Calculation Module
========================================================

Advanced depth calculation with quality filtering, normalization preparation,
and missing data detection for the SMN CNV detection pipeline.

Features:
- Per-exon coverage calculation with quality filtering
- Missing data detection and flagging
- Coverage statistics and quality metrics
- Compatibility with V3 normalization framework
- Robust handling of low-coverage regions

Author: SMN CNV Pipeline V3
Version: 3.0
Date: 2025-08-09
"""

import argparse
import os
import sys
import logging
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import tempfile
import json
from typing import Dict, List, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

class DepthCalculator:
    """Advanced depth calculation for SMN CNV detection pipeline V3."""
    
    def __init__(self, input_bam_dir: str, bed_file: str, output_dir: str, 
                 min_mapq: int = 20, min_baseq: int = 20, verbose: bool = False):
        """
        Initialize the depth calculator.
        
        Args:
            input_bam_dir: Directory containing BAM files
            bed_file: BED file with target regions (SMN exons)
            output_dir: Output directory for results
            min_mapq: Minimum mapping quality (default: 20)
            min_baseq: Minimum base quality (default: 20)
            verbose: Enable verbose logging
        """
        self.input_bam_dir = Path(input_bam_dir)
        self.bed_file = Path(bed_file)
        self.output_dir = Path(output_dir)
        self.min_mapq = min_mapq
        self.min_baseq = min_baseq
        self.verbose = verbose
        
        # V3 specific settings
        self.missing_data_threshold = 0.1  # 10% minimum coverage for valid exon
        self.quality_thresholds = {
            'excellent': 50,
            'good': 20,
            'poor': 5,
            'insufficient': 1
        }
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if verbose:
            logger.setLevel(logging.DEBUG)
            
        logger.info(f"Initialized DepthCalculator V3")
        logger.info(f"Input BAM directory: {self.input_bam_dir}")
        logger.info(f"BED file: {self.bed_file}")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Quality filters: MAPQ≥{self.min_mapq}, BaseQ≥{self.min_baseq}")

    def validate_inputs(self) -> bool:
        """Validate input files and directories."""
        logger.info("Validating inputs...")
        
        # Check BAM directory
        if not self.input_bam_dir.exists():
            logger.error(f"BAM directory not found: {self.input_bam_dir}")
            return False
            
        # Check for BAM files
        bam_files = list(self.input_bam_dir.glob("*.bam"))
        if not bam_files:
            logger.error(f"No BAM files found in {self.input_bam_dir}")
            return False
            
        logger.info(f"Found {len(bam_files)} BAM files")
        
        # Check BED file
        if not self.bed_file.exists():
            logger.error(f"BED file not found: {self.bed_file}")
            return False
            
        # Validate BED file format
        try:
            bed_df = pd.read_csv(self.bed_file, sep='\t', header=None, 
                               names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
            if len(bed_df) == 0:
                logger.error("BED file is empty")
                return False
            logger.info(f"BED file contains {len(bed_df)} regions")
        except Exception as e:
            logger.error(f"Invalid BED file format: {e}")
            return False
            
        # Check samtools availability
        try:
            subprocess.run(['samtools', '--version'], capture_output=True, check=True)
            logger.info("samtools found and accessible")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error("samtools not found or not accessible")
            return False
            
        return True

    def load_bed_regions(self) -> pd.DataFrame:
        """Load and parse BED file regions."""
        logger.info("Loading BED regions...")
        
        try:
            # Read BED file with flexible column handling
            bed_df = pd.read_csv(self.bed_file, sep='\t', header=None)
            
            # Handle different BED formats
            if len(bed_df.columns) >= 4:
                bed_df.columns = ['chrom', 'start', 'end', 'name'] + [f'col_{i}' for i in range(4, len(bed_df.columns))]
            else:
                bed_df.columns = ['chrom', 'start', 'end'] + [f'col_{i}' for i in range(3, len(bed_df.columns))]
                bed_df['name'] = bed_df.apply(lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1)
            
            # Ensure coordinates are integers
            bed_df['start'] = bed_df['start'].astype(int)
            bed_df['end'] = bed_df['end'].astype(int)
            
            # Add region length
            bed_df['length'] = bed_df['end'] - bed_df['start']
            
            logger.info(f"Loaded {len(bed_df)} regions")
            logger.debug(f"Regions: {bed_df['name'].tolist()}")
            
            return bed_df
            
        except Exception as e:
            logger.error(f"Failed to load BED regions: {e}")
            raise

    def calculate_sample_depth(self, bam_file: Path, bed_df: pd.DataFrame) -> Dict:
        """Calculate depth for a single sample."""
        sample_name = bam_file.stem
        logger.info(f"Processing sample: {sample_name}")
        
        try:
            # Check if BAM index exists
            index_file = bam_file.with_suffix('.bam.bai')
            if not index_file.exists():
                logger.warning(f"BAM index not found for {sample_name}, attempting to create...")
                subprocess.run(['samtools', 'index', str(bam_file)], check=True)
                
            sample_results = {
                'sample_id': sample_name,
                'regions': [],
                'total_reads': 0,
                'mapped_reads': 0,
                'quality_summary': defaultdict(int)
            }
            
            # Get basic BAM statistics
            try:
                stats_output = subprocess.run(
                    ['samtools', 'flagstat', str(bam_file)],
                    capture_output=True, text=True, check=True
                )
                
                for line in stats_output.stdout.split('\n'):
                    if 'total' in line and 'read1' not in line:
                        sample_results['total_reads'] = int(line.split()[0])
                    elif 'mapped' in line and '%' in line:
                        sample_results['mapped_reads'] = int(line.split()[0])
                        
            except subprocess.CalledProcessError:
                logger.warning(f"Could not get BAM statistics for {sample_name}")
            
            # Calculate depth for each region
            for _, region in bed_df.iterrows():
                region_result = self.calculate_region_depth(
                    bam_file, region, sample_name
                )
                sample_results['regions'].append(region_result)
                
                # Update quality summary
                quality_category = self.categorize_coverage_quality(region_result['mean_coverage'])
                sample_results['quality_summary'][quality_category] += 1
                
            logger.info(f"Completed sample {sample_name}: "
                       f"{len(sample_results['regions'])} regions processed")
            
            return sample_results
            
        except Exception as e:
            logger.error(f"Failed to process sample {sample_name}: {e}")
            raise

    def calculate_region_depth(self, bam_file: Path, region: pd.Series, sample_name: str) -> Dict:
        """Calculate depth for a specific region."""
        region_name = region['name']
        chrom = region['chrom']
        start = region['start']
        end = region['end']
        length = region['length']
        
        logger.debug(f"Processing region {region_name} for {sample_name}")
        
        try:
            # Use samtools depth with quality filtering
            cmd = [
                'samtools', 'depth',
                '-r', f"{chrom}:{start+1}-{end}",  # samtools uses 1-based coordinates
                '-q', str(self.min_mapq),
                '-Q', str(self.min_baseq),
                str(bam_file)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse depth output
            depths = []
            coverage_positions = set()
            
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        pos = int(parts[1])
                        depth = int(parts[2])
                        depths.append(depth)
                        coverage_positions.add(pos)
            
            # Calculate statistics
            if depths:
                mean_coverage = np.mean(depths)
                median_coverage = np.median(depths)
                std_coverage = np.std(depths)
                min_coverage = np.min(depths)
                max_coverage = np.max(depths)
                
                # Calculate coverage uniformity
                cv = std_coverage / mean_coverage if mean_coverage > 0 else float('inf')
                
                # Calculate fraction of bases covered
                covered_bases = len(coverage_positions)
                coverage_fraction = covered_bases / length
                
                # Detect missing data
                is_missing = coverage_fraction < (1 - self.missing_data_threshold)
                
            else:
                # No coverage detected
                mean_coverage = 0
                median_coverage = 0
                std_coverage = 0
                min_coverage = 0
                max_coverage = 0
                cv = float('inf')
                coverage_fraction = 0
                is_missing = True
                
            # Prepare region result
            region_result = {
                'sample_id': sample_name,
                'region_name': region_name,
                'chromosome': chrom,
                'start': start,
                'end': end,
                'length': length,
                'mean_coverage': round(mean_coverage, 2),
                'median_coverage': round(median_coverage, 2),
                'std_coverage': round(std_coverage, 2),
                'min_coverage': min_coverage,
                'max_coverage': max_coverage,
                'coefficient_of_variation': round(cv, 3) if cv != float('inf') else 'inf',
                'coverage_fraction': round(coverage_fraction, 3),
                'is_missing_data': is_missing,
                'quality_category': self.categorize_coverage_quality(mean_coverage)
            }
            
            logger.debug(f"Region {region_name}: mean_cov={mean_coverage:.1f}, "
                        f"fraction={coverage_fraction:.2f}, missing={is_missing}")
            
            return region_result
            
        except subprocess.CalledProcessError as e:
            logger.warning(f"samtools depth failed for {region_name} in {sample_name}: {e}")
            # Return zero coverage result
            return {
                'sample_id': sample_name,
                'region_name': region_name,
                'chromosome': chrom,
                'start': start,
                'end': end,
                'length': length,
                'mean_coverage': 0,
                'median_coverage': 0,
                'std_coverage': 0,
                'min_coverage': 0,
                'max_coverage': 0,
                'coefficient_of_variation': 'inf',
                'coverage_fraction': 0,
                'is_missing_data': True,
                'quality_category': 'insufficient'
            }

    def categorize_coverage_quality(self, mean_coverage: float) -> str:
        """Categorize coverage quality based on mean coverage."""
        if mean_coverage >= self.quality_thresholds['excellent']:
            return 'excellent'
        elif mean_coverage >= self.quality_thresholds['good']:
            return 'good'
        elif mean_coverage >= self.quality_thresholds['poor']:
            return 'poor'
        else:
            return 'insufficient'

    def save_coverage_summary(self, all_results: List[Dict]) -> None:
        """Save coverage summary in V3 format."""
        logger.info("Saving coverage summary...")
        
        # Prepare data for coverage summary
        summary_data = []
        for sample_result in all_results:
            for region_result in sample_result['regions']:
                summary_data.append({
                    'sample_id': region_result['sample_id'],
                    'exon': region_result['region_name'],  # V3 expects 'exon' column
                    'coverage': region_result['mean_coverage'],
                    'normalized_coverage': region_result['mean_coverage']  # Will be updated by normalization step
                })
        
        # Save main coverage summary (required by V3 pipeline)
        summary_df = pd.DataFrame(summary_data)
        summary_file = self.output_dir / 'coverage_summary.txt'
        summary_df.to_csv(summary_file, sep='\t', index=False)
        logger.info(f"Saved coverage summary: {summary_file}")
        
        # Save detailed coverage metrics for V3 quality assessment
        detailed_data = []
        for sample_result in all_results:
            for region_result in sample_result['regions']:
                detailed_data.append(region_result)
        
        detailed_df = pd.DataFrame(detailed_data)
        detailed_file = self.output_dir / 'detailed_coverage_metrics.txt'
        detailed_df.to_csv(detailed_file, sep='\t', index=False)
        logger.info(f"Saved detailed metrics: {detailed_file}")

    def save_quality_report(self, all_results: List[Dict]) -> None:
        """Save quality assessment report."""
        logger.info("Generating quality assessment report...")
        
        quality_report = {
            'pipeline_version': '3.0',
            'analysis_date': pd.Timestamp.now().isoformat(),
            'total_samples': len(all_results),
            'quality_filters': {
                'min_mapq': self.min_mapq,
                'min_baseq': self.min_baseq,
                'missing_data_threshold': self.missing_data_threshold
            },
            'sample_summary': [],
            'overall_statistics': {}
        }
        
        # Collect sample-level statistics
        all_coverages = []
        missing_data_samples = 0
        
        for sample_result in all_results:
            sample_coverages = [r['mean_coverage'] for r in sample_result['regions']]
            missing_regions = sum(1 for r in sample_result['regions'] if r['is_missing_data'])
            
            sample_summary = {
                'sample_id': sample_result['sample_id'],
                'total_reads': sample_result['total_reads'],
                'mapped_reads': sample_result['mapped_reads'],
                'mean_coverage': round(np.mean(sample_coverages), 2),
                'median_coverage': round(np.median(sample_coverages), 2),
                'missing_regions': missing_regions,
                'quality_distribution': dict(sample_result['quality_summary'])
            }
            
            quality_report['sample_summary'].append(sample_summary)
            all_coverages.extend(sample_coverages)
            
            if missing_regions > 0:
                missing_data_samples += 1
        
        # Overall statistics
        quality_report['overall_statistics'] = {
            'mean_coverage_all_regions': round(np.mean(all_coverages), 2),
            'median_coverage_all_regions': round(np.median(all_coverages), 2),
            'std_coverage_all_regions': round(np.std(all_coverages), 2),
            'samples_with_missing_data': missing_data_samples,
            'missing_data_rate': round(missing_data_samples / len(all_results), 3)
        }
        
        # Save quality report
        quality_file = self.output_dir / 'depth_quality_report.json'
        with open(quality_file, 'w') as f:
            json.dump(quality_report, f, indent=2)
            
        logger.info(f"Saved quality report: {quality_file}")
        
        # Also save a human-readable summary
        summary_file = self.output_dir / 'depth_quality_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("SMN CNV Pipeline V3 - Depth Calculation Quality Summary\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Analysis Date: {quality_report['analysis_date']}\n")
            f.write(f"Total Samples: {quality_report['total_samples']}\n")
            f.write(f"Mean Coverage (all regions): {quality_report['overall_statistics']['mean_coverage_all_regions']}\n")
            f.write(f"Samples with Missing Data: {quality_report['overall_statistics']['samples_with_missing_data']}\n")
            f.write(f"Missing Data Rate: {quality_report['overall_statistics']['missing_data_rate']:.1%}\n\n")
            
            f.write("Quality Filters Applied:\n")
            f.write(f"  Minimum Mapping Quality: {self.min_mapq}\n")
            f.write(f"  Minimum Base Quality: {self.min_baseq}\n")
            f.write(f"  Missing Data Threshold: {self.missing_data_threshold:.1%}\n\n")
            
            f.write("Sample-wise Summary:\n")
            for sample in quality_report['sample_summary']:
                f.write(f"  {sample['sample_id']}: "
                       f"mean_cov={sample['mean_coverage']}, "
                       f"missing={sample['missing_regions']}\n")
        
        logger.info(f"Saved readable summary: {summary_file}")

    def run(self) -> bool:
        """Run the complete depth calculation pipeline."""
        logger.info("Starting SMN CNV Pipeline V3 depth calculation...")
        
        try:
            # Validate inputs
            if not self.validate_inputs():
                return False
            
            # Load BED regions
            bed_df = self.load_bed_regions()
            
            # Get list of BAM files
            bam_files = sorted(self.input_bam_dir.glob("*.bam"))
            logger.info(f"Processing {len(bam_files)} BAM files...")
            
            # Process each BAM file
            all_results = []
            for i, bam_file in enumerate(bam_files, 1):
                logger.info(f"Processing file {i}/{len(bam_files)}: {bam_file.name}")
                sample_result = self.calculate_sample_depth(bam_file, bed_df)
                all_results.append(sample_result)
            
            # Save results
            self.save_coverage_summary(all_results)
            self.save_quality_report(all_results)
            
            logger.info("Depth calculation completed successfully!")
            logger.info(f"Results saved to: {self.output_dir}")
            
            return True
            
        except Exception as e:
            logger.error(f"Depth calculation failed: {e}")
            if self.verbose:
                import traceback
                logger.error(traceback.format_exc())
            return False


def main():
    """Main entry point for the depth calculation script."""
    parser = argparse.ArgumentParser(
        description="SMN CNV Pipeline V3 - Advanced Depth Calculation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic depth calculation
  python calculate_depth.py /path/to/bam/files/ smn_exons.bed /output/dir/
  
  # With custom quality filters
  python calculate_depth.py /path/to/bam/files/ smn_exons.bed /output/dir/ \\
    --min-mapq 30 --min-baseq 25 --verbose
        """
    )
    
    parser.add_argument(
        'input_bam_dir',
        help='Directory containing BAM files to analyze'
    )
    
    parser.add_argument(
        'bed_file',
        help='BED file with target regions (SMN exons)'
    )
    
    parser.add_argument(
        'output_dir',
        help='Output directory for depth calculation results'
    )
    
    parser.add_argument(
        '--min-mapq',
        type=int,
        default=20,
        help='Minimum mapping quality for reads (default: 20)'
    )
    
    parser.add_argument(
        '--min-baseq',
        type=int,
        default=20,
        help='Minimum base quality for positions (default: 20)'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='SMN CNV Pipeline V3 Depth Calculator 3.0'
    )
    
    args = parser.parse_args()
    
    # Create depth calculator
    calculator = DepthCalculator(
        input_bam_dir=args.input_bam_dir,
        bed_file=args.bed_file,
        output_dir=args.output_dir,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        verbose=args.verbose
    )
    
    # Run calculation
    success = calculator.run()
    
    if success:
        logger.info("Depth calculation completed successfully!")
        sys.exit(0)
    else:
        logger.error("Depth calculation failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
