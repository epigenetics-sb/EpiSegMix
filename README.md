# EpiSegMix Pipeline
---
This repository contains a Nextflow workflow for chromatin segmentation to annotate both coding and non-coding regions of the genome.
EpiSegMix first estimates the parameters of a hidden Markov model, where each state corresponds to a different combination of epigenetic modifications and thus represents a functional role, such as enhancer, transcription start site, active or silent gene. The spatial relations are captured via the transition probabolities. After the parameter estimation, each region in the genome is annotated with the most likely chromatin state. The implementation allows to choose for each histone modification a different distributional assumption (the available distributions are listed below). Similar tools are ChromHMM or EpiCSeg (references [2] and [3]).
The implementation of the HMM is in C++. The parameters are estimated using the Baum-Welch algorithm and for decoding the user can select either the Viterbi algorithm or posterior decoding. For all distributions where the MLE contains no closed form solution, the parameters are updated by numerical optimization using the Migrad minimizer from ROOT (https://root.cern/root/htmldoc/guides/minuit2/Minuit2Letter.pdf). The initial parameters of the HMM are found using k-Means clustering. Each cluster is assumed to correspond to one state and the parameters are initialized using method of moment estimates.


## Quick Start

1. **Install Nextflow** (≥22.10) and **Docker** (or Singularity).
2. **Clone & Run Test:**

```bash
git clone https://github.com/epigenetics-sb/EpiSegMix.git 
cd EpiSegMix 

# Run with test data
nextflow run main.nf -profile test,docker
# or
nextflow run main.nf -profile test,singularity
```

---

## Input: `samplesheet.csv`

Create a file named `samplesheet.csv` with these exact columns:

| Column | Description | Example |
| --- | --- | --- |
| `sample_id` | Unique ID for the sample | `sample_01` |
| `replicate` | Replicate number | `1` |
| `epigenetic_mark` | Mark name (H3K4me3, WGBS, etc.) | `H3K27ac` |
| `file_name` | Absolute path to BAM or BED file | `/data/file.bam` |
| `modality` | `ChIP-seq`, `WGBS`, `ATAC-seq` | `ChIP-seq` |
| `paired_end` | `true` or `false` | `true` |

> **Note:** File paths must exist on your filesystem and match the modality:
> - ChIP-seq / ATAC-seq → `.bam`  
> - WGBS / NOMe-seq → `.bed` or `.bed.gz`  

---

## Usage

### Standard Mode (Default)

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --episegmix_mode standard \
    --states 8 \
    --outdir ./results
```

### Duration Modeling Mode

Uses topology-based HMM to model state duration.

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --episegmix_mode duration \
    --states 8 \
    --outdir ./results_dm
```

---

##  Key Outputs

Results are saved in your output directory:

* **`EpiSegMix/{id}/Plots/{id}.html`**: **Main Summary Report** (start here).  
* **`EpiSegMix/{id}/Segmentation/{id}.bed.gz`**:  Final chromatin segmentation file.  
* **`EpiSegMix/{id}/Models/{id}.model.json`**:  Trained model parameters.

---

## Key Parameters

| Flag | Default | Description |
| --- | --- | --- |
| `--episegmix_mode` | `standard` | Select `standard` or `duration`. |
| `--states` | `8` | Number of chromatin states. |
| `--genome` | `hg38` | Reference genome (e.g., `mm10`, `hg19`). |
| `--merge` | `false` | Merge replicates/modalities before training. |
| `--binsize` | `200` | Genome bin size in bp. |
| `--dist_histone` | `NBI` | Distribution for histone counts (Poisson, BI, NBI, etc.). |
| `--dist_methyl` | `BI` | Distribution for methylation counts (BI or BB only). |
| `--chr_parameter_estimation` | `pilot_hg38` | Chromosome parameter estimation (string or int). |
| `--decoding_algorithm` | `viterbi` | Algorithm for decoding states (`viterbi`). |
| `--outdir` | `./results` | Output directory for all results. |

---

## Notes

* Use **absolute paths** for input files.  
* Docker or Singularity will automatically mount the pipeline directory.  
* For testing, use the `-profile test,docker` profile.  
