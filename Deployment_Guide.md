# SMN CNV Detection Pipeline V3 - Deployment Guide

## Quick Start - V3 Deployment

### 1. System Requirements

**Hardware Requirements:**
- CPU: 4+ cores recommended (8+ cores for large cohorts)
- RAM: 8GB minimum, 16GB+ recommended
- Storage: 10GB+ free space for pipeline and results
- OS: Linux/Unix (Ubuntu 18.04+, CentOS 7+, macOS)

**Software Requirements:**
- samtools ≥1.10
- Python 3.7+
- Git (for installation)

### 2. Installation

```bash
# Clone the repository
git clone <repository_url>
cd smn_cnv_pipeline

# Install V3 dependencies
pip install -r requirements_v3.txt

# Make scripts executable
chmod +x run_pipeline_v3.sh bin/*.sh bin/*.py

# Verify installation
python test_v3_integration.py --skip-deps
```

### 3. Quick V3 Test

```bash
# Create test data
python test_v3_integration.py --create-test-data

# Run V3 on test data (takes ~5 minutes)
./run_pipeline_v3.sh test_data_v3/bam_files/ --skip-plots

# Check results
ls results_v3/
cat results_v3/v3_pipeline_summary.txt
```

### 4. Production Usage

```bash
# Full V3 analysis with all features
./run_pipeline_v3.sh /path/to/your/bam/files/

# View results
firefox results/reports_v3/SAMPLE_ID/SAMPLE_ID_report.html
```

---

## Detailed V3 Deployment

### Environment Setup

#### Option A: Conda Environment (Recommended)

```bash
# Create conda environment
conda create -n smn_cnv_v3 python=3.9
conda activate smn_cnv_v3

# Install core dependencies
conda install pandas numpy scipy matplotlib seaborn scikit-learn pyyaml

# Install samtools
conda install -c bioconda samtools

# Verify installation
python -c "import pandas, numpy, sklearn, yaml; print('V3 dependencies OK')"
```

#### Option B: System Installation

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install samtools python3 python3-pip

# CentOS/RHEL
sudo yum install samtools python3 python3-pip

# Install Python packages
pip3 install -r requirements_v3.txt
```

#### Option C: Docker Deployment

```dockerfile
# Dockerfile for V3
FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    samtools \
    python3 \
    python3-pip \
    git

COPY requirements_v3.txt .
RUN pip3 install -r requirements_v3.txt

COPY . /pipeline
WORKDIR /pipeline
RUN chmod +x run_pipeline_v3.sh bin/*.sh bin/*.py

CMD ["bash"]
```

Build and run:
```bash
docker build -t smn-cnv-v3 .
docker run -v /path/to/data:/data smn-cnv-v3 ./run_pipeline_v3.sh /data/bam_files/
```

### Configuration Customization

#### V3 Parameters File (`config/segmentation_params.yaml`)

```yaml
# Example customization for clinical lab
cbs:
  alpha: 0.01          # More stringent for clinical use
  min_segment_markers: 3
  nperm: 10000

hmm:
  n_states: 5
  max_iter: 100
  
consensus:
  min_algorithms: 2    # Require both CBS and HMM
  confidence_threshold: 0.8  # Higher confidence for clinical

quality_control:
  coverage:
    min_depth: 20      # Clinical minimum
    max_cv: 1.5        # Stricter variability
```

#### Control Genes Configuration

Add lab-specific control genes to `config/control_genes.bed`:

```bash
# Add your validated control genes
echo -e "chr1\t1234567\t1234890\tCUSTOM_CONTROL1\t.\t+" >> config/control_genes.bed
```

### Batch Processing Setup

#### Large Cohort Processing

```bash
#!/bin/bash
# process_cohort.sh - Process multiple cohorts

cohorts=("cohort1" "cohort2" "cohort3")

for cohort in "${cohorts[@]}"; do
    echo "Processing $cohort..."
    
    ./run_pipeline_v3.sh /data/${cohort}/bam_files/ \
        --results /results/${cohort}/ \
        --sample-type auto \
        --consensus \
        --min-algorithms 2
        
    # Generate cohort summary
    python bin/cohort_summary.py /results/${cohort}/ > /results/${cohort}_summary.txt
done
```

#### Parallel Sample Processing

```bash
# Process samples in parallel (GNU parallel)
find /data/bam_files/ -name "*.bam" | \
parallel -j 4 ./run_pipeline_v3.sh {} --sample-type test --skip-plots
```

### Quality Control Setup

#### Pre-Analysis QC

```bash
#!/bin/bash
# pre_analysis_qc.sh

bam_dir=$1
echo "Pre-analysis QC for: $bam_dir"

# Check BAM files
for bam in $bam_dir/*.bam; do
    if [ ! -f "$bam.bai" ]; then
        echo "WARNING: Missing index for $bam"
        samtools index "$bam"
    fi
    
    # Basic statistics
    samtools flagstat "$bam" > "${bam}_flagstat.txt"
done

# Check coverage in SMN regions
samtools depth -b config/smn_exons.bed $bam_dir/*.bam > smn_region_coverage.txt
python bin/coverage_qc.py smn_region_coverage.txt
```

#### Post-Analysis QC

```bash
# Automated QC checks
python bin/v3_quality_checker.py results/ --thresholds clinical
```

### Clinical Laboratory Integration

#### LIMS Integration Template

```python
# lims_integration.py
import json
import subprocess
from pathlib import Path

def process_clinical_sample(sample_id, bam_path, patient_info):
    """Process single clinical sample with V3 pipeline"""
    
    # Create sample-specific directory
    sample_dir = Path(f"/clinical_results/{sample_id}")
    sample_dir.mkdir(exist_ok=True)
    
    # Copy BAM to processing area
    temp_bam_dir = sample_dir / "input"
    temp_bam_dir.mkdir(exist_ok=True)
    
    # Run V3 pipeline
    cmd = [
        "./run_pipeline_v3.sh",
        str(temp_bam_dir),
        "--results", str(sample_dir / "results"),
        "--consensus",
        "--confidence", "0.8",
        "--sample-type", "test"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        # Parse results for LIMS
        return parse_clinical_results(sample_dir / "results")
    else:
        raise Exception(f"Pipeline failed: {result.stderr}")

def parse_clinical_results(results_dir):
    """Parse V3 results for clinical reporting"""
    consensus_file = results_dir / "consensus" / "high_confidence_calls.txt"
    
    if consensus_file.exists():
        import pandas as pd
        df = pd.read_csv(consensus_file, sep='\t')
        
        # Extract SMN1/SMN2 copy numbers
        clinical_result = {
            'SMN1_copy_number': extract_smn_cn(df, 'SMN1'),
            'SMN2_copy_number
