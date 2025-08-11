#!/usr/bin/env python3

"""
test_v3_integration.py - Integration tests for SMN CNV Pipeline V3
Usage: python test_v3_integration.py [--create-test-data]
"""

import sys
import os
import tempfile
import shutil
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import argparse

#!/usr/bin/env python3

"""
test_v3_integration.py - Integration tests for SMN CNV Pipeline V3
Usage: python test_v3_integration.py [--create-test-data]
"""

import sys
import os
import tempfile
import shutil
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import argparse
import yaml

def create_test_data(output_dir):
    """Create synthetic test data for V3 pipeline testing"""
    test_dir = Path(output_dir)
    test_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Creating test data in: {test_dir}")
    
    # Create synthetic coverage data
    samples = ['ref001', 'ref002', 'ref003', 'control001', 'normal001', 
              'patient001', 'patient002', 'test001', 'sample_missing_exon7']
    exons = ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']
    
    # Generate realistic coverage data
    np.random.seed(42)
    coverage_data = {}
    
    for sample in samples:
        coverage_data[sample] = {}
        
        if 'ref' in sample or 'control' in sample or 'normal' in sample:
            # Reference samples - normal copy number (CN=2)
            base_coverage = np.random.normal(50, 10, len(exons))
            base_coverage = np.maximum(base_coverage, 5)  # Minimum coverage
            
        elif 'patient001' in sample:
            # SMA affected - homozygous deletion (CN=0) for SMN1
            base_coverage = []
            for exon in exons:
                if 'SMN1' in exon:
                    base_coverage.append(np.random.normal(2, 1))  # Very low coverage
                else:
                    base_coverage.append(np.random.normal(50, 10))  # Normal SMN2
            base_coverage = np.maximum(base_coverage, 0)
            
        elif 'patient002' in sample:
            # SMA carrier - heterozygous deletion (CN=1) for SMN1
            base_coverage = []
            for exon in exons:
                if 'SMN1' in exon:
                    base_coverage.append(np.random.normal(25, 8))  # Half coverage
                else:
                    base_coverage.append(np.random.normal(50, 10))  # Normal SMN2
            base_coverage = np.maximum(base_coverage, 2)
            
        elif 'missing_exon7' in sample:
            # Sample with missing exon 7 data
            base_coverage = []
            for exon in exons:
                if 'exon7' in exon:
                    base_coverage.append(np.nan)  # Missing data
                else:
                    base_coverage.append(np.random.normal(45, 12))
            
        else:
            # Test samples - normal with some variation
            base_coverage = np.random.normal(45, 12, len(exons))
            base_coverage = np.maximum(base_coverage, 5)
        
        for i, exon in enumerate(exons):
            coverage_data[sample][exon] = base_coverage[i]
    
    # Create coverage summary file
    coverage_df = pd.DataFrame(coverage_data).T
    coverage_file = test_dir / 'coverage_summary.txt'
    coverage_df.to_csv(coverage_file, sep='\t')
    
    # Create coverage pivot file
    pivot_file = test_dir / 'coverage_summary_pivot.txt'
    coverage_df.to_csv(pivot_file, sep='\t')
    
    # Create sample info file
    sample_info = []
    for sample in samples:
        if any(keyword in sample.lower() for keyword in ['ref', 'control', 'normal']):
            sample_type = 'reference'
        else:
            sample_type = 'test'
        
        sample_info.append({
            'sample_id': sample,
            'bam_path': f'/fake/path/{sample}.bam',
            'sample_type': sample_type
        })
    
    sample_info_df = pd.DataFrame(sample_info)
    sample_info_file = test_dir / 'sample_info.txt'
    sample_info_df.to_csv(sample_info_file, sep='\t', index=False)
    
    # Create synthetic allele count data
    allele_data = []
    for sample in samples:
        for gene in ['SMN1', 'SMN2']:
            for snp in ['c.840C>T', 'c.888+100G>A']:
                if 'patient001' in sample and gene == 'SMN1':
                    # Homozygous deletion - very low counts
                    ref_count = np.random.poisson(1)
                    alt_count = np.random.poisson(1)
                elif 'patient002' in sample and gene == 'SMN1':
                    # Heterozygous deletion - intermediate counts
                    ref_count = np.random.poisson(15)
                    alt_count = np.random.poisson(15)
                else:
                    # Normal samples
                    ref_count = np.random.poisson(30)
                    alt_count = np.random.poisson(30)
                
                total_depth = ref_count + alt_count
                ref_freq = ref_count / max(1, total_depth)
                alt_freq = alt_count / max(1, total_depth)
                
                allele_data.append({
                    'sample_id': sample,
                    'chrom': 'chr5',
                    'pos': 70925001 if snp == 'c.840C>T' else 70924784,
                    'snp_name': snp,
                    'gene': gene,
                    'ref_allele': 'C' if snp == 'c.840C>T' else 'G',
                    'alt_allele': 'T' if snp == 'c.840C>T' else 'A',
                    'ref_count': ref_count,
                    'alt_count': alt_count,
                    'total_depth': total_depth,
                    'other_count': 0,
                    'ref_freq': ref_freq,
                    'alt_freq': alt_freq,
                    'sample_type': sample_info_df[sample_info_df['sample_id'] == sample]['sample_type'].iloc[0]
                })
    
    allele_df = pd.DataFrame(allele_data)
    allele_file = test_dir / 'allele_counts.txt'
    allele_df.to_csv(allele_file, sep='\t', index=False)
    
    # Create control normalization factors
    control_factors = []
    for sample in samples:
        # Add some realistic variation in normalization factors
        factor = np.random.normal(1.0, 0.1)
        factor = max(0.5, min(2.0, factor))  # Reasonable range
        
        control_factors.append({
            'sample_id': sample,
            'normalization_factor': factor,
            'method': 'geometric_mean'
        })
    
    control_factors_df = pd.DataFrame(control_factors)
    control_factors_file = test_dir / 'normalization_factors.txt'
    control_factors_df.to_csv(control_factors_file, sep='\t', index=False)
    
    print(f"Test data created:")
    print(f"  Samples: {len(samples)}")
    print(f"  Exons: {len(exons)}")
    print(f"  Coverage file: {coverage_file}")
    print(f"  Sample info: {sample_info_file}")
    print(f"  Allele counts: {allele_file}")
    print(f"  Control factors: {control_factors_file}")
    
    return test_dir

def test_control_gene_analysis(test_data_dir, output_dir):
    """Test control gene analysis component"""
    print("\n=== Testing Control Gene Analysis ===")
    
    # Create a minimal control genes BED file for testing
    control_genes_content = """# Test control genes
chr12	6534405	6538374	GAPDH	.	-
chr7	5527148	5530601	ACTB	.	-
chr7	44838788	44846306	PPIA	.	+"""
    
    control_genes_file = test_data_dir / 'control_genes.bed'
    with open(control_genes_file, 'w') as f:
        f.write(control_genes_content)
    
    coverage_file = test_data_dir / 'coverage_summary.txt'
    control_output_dir = output_dir / 'control_analysis_test'
    
    cmd = [
        'python3', 'bin/control_gene_analysis.py',
        str(coverage_file),
        str(control_genes_file), 
        str(control_output_dir),
        '--no-plots'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✓ Control gene analysis completed successfully")
        
        # Check expected output files
        expected_files = [
            'control_selection.txt',
            'normalization_factors.txt', 
            'gene_stability_metrics.txt'
        ]
        
        for file_name in expected_files:
            file_path = control_output_dir / file_name
            if file_path.exists():
                print(f"✓ Output file created: {file_name}")
            else:
                print(f"✗ Missing output file: {file_name}")
                
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Control gene analysis failed: {e}")
        print(f"  stdout: {e.stdout}")
        print(f"  stderr: {e.stderr}")
        return False

def test_cbs_segmentation(test_data_dir, output_dir):
    """Test CBS segmentation component"""
    print("\n=== Testing CBS Segmentation ===")
    
    coverage_file = test_data_dir / 'coverage_summary_pivot.txt'
    cbs_output_dir = output_dir / 'cbs_test'
    
    cmd = [
        'python3', 'bin/cbs_segmentation.py',
        str(coverage_file),
        str(cbs_output_dir),
        '--alpha', '0.05',  # Less stringent for test
        '--min-markers', '2',
        '--nperm', '100',   # Fewer permutations for speed
        '--no-plots'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✓ CBS segmentation completed successfully")
        
        # Check output files
        segments_file = cbs_output_dir / 'cbs_segments.txt'
        params_file = cbs_output_dir / 'cbs_parameters.txt'
        
        if segments_file.exists():
            segments_df = pd.read_csv(segments_file, sep='\t')
            print(f"✓ CBS segments file created with {len(segments_df)} segments")
        else:
            print("✗ CBS segments file not created")
            
        if params_file.exists():
            print("✓ CBS parameters file created")
        else:
            print("✗ CBS parameters file not created")
            
        return segments_file.exists()
        
    except subprocess.CalledProcessError as e:
        print(f"✗ CBS segmentation failed: {e}")
        print(f"  stdout: {e.stdout}")
        print(f"  stderr: {e.stderr}")
        return False

def test_hmm_segmentation(test_data_dir, output_dir):
    """Test HMM segmentation component"""
    print("\n=== Testing HMM Segmentation ===")
    
    coverage_file = test_data_dir / 'coverage_summary_pivot.txt'
    hmm_output_dir = output_dir / 'hmm_test'
    
    cmd = [
        'python3', 'bin/hmm_segmentation.py',
        str(coverage_file),
        str(hmm_output_dir),
        '--states', '5',
        '--max-iter', '20',  # Fewer iterations for speed
        '--interpolate',
        '--no-plots'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✓ HMM segmentation completed successfully")
        
        # Check output files
        segments_file = hmm_output_dir / 'hmm_segments.txt'
        params_file = hmm_output_dir / 'hmm_parameters.txt'
        interp_file = hmm_output_dir / 'interpolated_coverage.txt'
        
        if segments_file.exists():
            segments_df = pd.read_csv(segments_file, sep='\t')
            print(f"✓ HMM segments file created with {len(segments_df)} segments")
        else:
            print("✗ HMM segments file not created")
            
        if params_file.exists():
            print("✓ HMM parameters file created")
        else:
            print("✗ HMM parameters file not created")
            
        if interp_file.exists():
            print("✓ Interpolated coverage file created")
        else:
            print("✗ Interpolated coverage file not created")
            
        return segments_file.exists()
        
    except subprocess.CalledProcessError as e:
        print(f"✗ HMM segmentation failed: {e}")
        print(f"  stdout: {e.stdout}")
        print(f"  stderr: {e.stderr}")
        return False

def test_consensus_calling(output_dir):
    """Test consensus calling component"""
    print("\n=== Testing Consensus Calling ===")
    
    cbs_file = output_dir / 'cbs_test' / 'cbs_segments.txt'
    hmm_file = output_dir / 'hmm_test' / 'hmm_segments.txt'
    consensus_output_dir = output_dir / 'consensus_test'
    
    if not (cbs_file.exists() and hmm_file.exists()):
        print("✗ Missing input files for consensus calling")
        return False
    
    cmd = [
        'python3', 'bin/consensus_calling.py',
        str(cbs_file),
        str(hmm_file),
        str(consensus_output_dir),
        '--min-algorithms', '1',  # Accept single algorithm for test
        '--confidence-threshold', '0.5',
        '--no-plots'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✓ Consensus calling completed successfully")
        
        # Check output files
        consensus_file = consensus_output_dir / 'consensus_calls.txt'
        high_conf_file = consensus_output_dir / 'high_confidence_calls.txt'
        metrics_file = consensus_output_dir / 'quality_metrics.txt'
        
        if consensus_file.exists():
            consensus_df = pd.read_csv(consensus_file, sep='\t')
            print(f"✓ Consensus calls file created with {len(consensus_df)} calls")
        else:
            print("✗ Consensus calls file not created")
            
        if high_conf_file.exists():
            high_conf_df = pd.read_csv(high_conf_file, sep='\t')
            print(f"✓ High confidence calls file created with {len(high_conf_df)} calls")
        else:
            print("✗ High confidence calls file not created")
            
        if metrics_file.exists():
            print("✓ Quality metrics file created")
        else:
            print("✗ Quality metrics file not created")
            
        return consensus_file.exists()
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Consensus calling failed: {e}")
        print(f"  stdout: {e.stdout}")
        print(f"  stderr: {e.stderr}")
        return False

def test_missing_data_handler(test_data_dir, output_dir):
    """Test missing data handler component"""
    print("\n=== Testing Missing Data Handler ===")
    
    coverage_file = test_data_dir / 'coverage_summary_pivot.txt'
    missing_output_dir = output_dir / 'missing_data_test'
    
    cmd = [
        'python3', 'bin/missing_data_handler.py',
        str(coverage_file),
        str(missing_output_dir),
        '--method', 'linear',  # Use faster method for test
        '--validate',
        '--no-plots'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✓ Missing data handler completed successfully")
        
        # Check output files
        interp_file = missing_output_dir / 'interpolated_coverage.txt'
        mask_file = missing_output_dir / 'missing_data_mask.txt'
        analysis_file = missing_output_dir / 'missing_data_analysis.txt'
        quality_file = missing_output_dir / 'interpolation_quality.txt'
        
        files_to_check = [
            ('Interpolated coverage', interp_file),
            ('Missing data mask', mask_file),
            ('Missing data analysis', analysis_file),
            ('Interpolation quality', quality_file)
        ]
        
        all_exist = True
        for name, file_path in files_to_check:
            if file_path.exists():
                print(f"✓ {name} file created")
            else:
                print(f"✗ {name} file not created")
                all_exist = False
                
        return all_exist
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Missing data handler failed: {e}")
        print(f"  stdout: {e.stdout}")
        print(f"  stderr: {e.stderr}")
        return False

def test_v3_integration(test_data_dir):
    """Test complete V3 integration"""
    print("\n=== Testing V3 Integration ===")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        
        # Test individual components
        tests = [
            ('Control Gene Analysis', lambda: test_control_gene_analysis(test_data_dir, output_dir)),
            ('CBS Segmentation', lambda: test_cbs_segmentation(test_data_dir, output_dir)),
            ('HMM Segmentation', lambda: test_hmm_segmentation(test_data_dir, output_dir)),
            ('Missing Data Handler', lambda: test_missing_data_handler(test_data_dir, output_dir)),
            ('Consensus Calling', lambda: test_consensus_calling(output_dir))
        ]
        
        results = {}
        for test_name, test_func in tests:
            try:
                results[test_name] = test_func()
            except Exception as e:
                print(f"✗ {test_name} failed with exception: {e}")
                results[test_name] = False
        
        # Summary
        print(f"\n=== V3 Integration Test Summary ===")
        passed = sum(results.values())
        total = len(results)
        
        for test_name, passed in results.items():
            status = "✓ PASS" if passed else "✗ FAIL"
            print(f"  {status}: {test_name}")
        
        print(f"\nOverall: {passed}/{total} tests passed")
        
        if passed == total:
            print("🎉 All V3 components working correctly!")
            return True
        else:
            print("❌ Some V3 components need attention")
            return False

def validate_v3_dependencies():
    """Validate that all V3 dependencies are available"""
    print("=== Validating V3 Dependencies ===")
    
    # Check Python packages
    required_packages = [
        'pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn', 
        'sklearn', 'yaml'
    ]
    
    missing_packages = []
    for package in required_packages:
        try:
            if package == 'sklearn':
                import sklearn
            elif package == 'yaml':
                import yaml
            else:
                __import__(package)
            print(f"✓ {package}")
        except ImportError:
            print(f"✗ {package} (missing)")
            missing_packages.append(package)
    
    if missing_packages:
        print(f"\n❌ Missing packages: {', '.join(missing_packages)}")
        print("Install with: pip install scikit-learn PyYAML")
        return False
    else:
        print("\n✓ All required Python packages available")
        return True

def main():
    parser = argparse.ArgumentParser(description='V3 Integration Tests')
    parser.add_argument('--create-test-data', action='store_true',
                       help='Create test data and exit')
    parser.add_argument('--test-data-dir', help='Directory with existing test data')
    parser.add_argument('--skip-deps', action='store_true',
                       help='Skip dependency validation')
    
    args = parser.parse_args()
    
    if not args.skip_deps:
        if not validate_v3_dependencies():
            print("❌ Dependency validation failed")
            sys.exit(1)
    
    if args.create_test_data:
        test_data_dir = create_test_data('test_data_v3')
        print(f"\n✓ Test data created in: {test_data_dir}")
        sys.exit(0)
    
    # Use provided test data directory or create temporary one
    if args.test_data_dir:
        test_data_dir = Path(args.test_data_dir)
        if not test_data_dir.exists():
            print(f"❌ Test data directory not found: {test_data_dir}")
            sys.exit(1)
    else:
        with tempfile.TemporaryDirectory() as temp_dir:
            test_data_dir = create_test_data(temp_dir + '/test_data')
    
    # Run integration tests
    success = test_v3_integration(test_data_dir)
    
    if success:
        print("\n🎉 V3 Integration Tests: ALL PASSED")
        print("SMN CNV Detection Pipeline V3 is ready for production use!")
        sys.exit(0)
    else:
        print("\n❌ V3 Integration Tests: SOME FAILED")
        print("Please check the failed components before using V3 in production")
        sys.exit(1)

if __name__ == "__main__":
    main()
