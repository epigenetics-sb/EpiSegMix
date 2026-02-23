# **EpiSegMix**

**A Nextflow workflow for chromatin segmentation to annotate coding and non-coding genomic regions.**

## **Introduction**

**EpiSegMix** is a Nextflow pipeline for chromatin segmentation. It uses a hidden Markov model (HMM) to annotate genomic regions with functional states (e.g., enhancers, promoters) based on combinations of epigenetic modifications, capturing spatial relations via transition probabilities.  
The C++ backend allows users to select flexible statistical distributions for histone modifications. Parameters are optimized via Baum-Welch training and ROOT's Migrad minimizer, followed by Viterbi or posterior decoding to assign the most likely chromatin states.

## **Pipeline Summary**

EpiSegMix uses a modular architecture that dynamically executes subworkflows based on your inputs:

* **Setup & Pre-processing** (input\_check.nf, prepare\_genome.nf, generate\_bins.nf): Validates the samplesheet, prepares reference genomes, and creates fixed-size genomic bins.  
* **Data Processing** (process\_histones.nf, process\_methyl.nf): Dynamically processes provided .bam (ChIP/ATAC-seq) and .bed (WGBS/NOMe-seq) files.  
* **Data Merging** (merge\_data.nf): Combines modalities prior to training if \--merge true is set.  
* **Model Training & Decoding**: Runs one of four modules based on the \--episegmix\_mode flag:  
  * **standard** (model\_training\_std.nf): Standard Baum-Welch training and decoding.  
  * **duration** (model\_training\_dm.nf): Topology-based duration modeling HMM.  
  * **DNA** (model\_training\_dna.nf): Segmentation solely on DNA methylation data.  
  * **fitting** (distribution\_fitting.nf): Empirically determines best-fitting distributions and generates an optimized samplesheet for downstream training.

## **Quick Start**

1. Install [Nextflow](https://www.google.com/search?q=https://www.nextflow.io/docs/latest/getstarted.html%23installation) (\>=22.10.1) and either [Docker](https://docs.docker.com/engine/installation/) or [Singularity](https://singularity.lbl.gov/all-releases).  
2. Clone the repository:  
   git clone [https://github.com/epigenetics-sb/EpiSegMix.git](https://github.com/epigenetics-sb/EpiSegMix.git)  
   cd EpiSegMix

3. Run the test profile:  
   nextflow run main.nf \-profile test,docker

4. **Standard Analysis:**  
   nextflow run main.nf \-profile docker \\  
       \--input samplesheet.csv \\  
       \--episegmix\_mode standard \\  
       \--states 8 \\  
       \--outdir ./results

5. **Advanced Run (Duration Mode \+ Merging):**  
   nextflow run main.nf \-profile docker \\  
       \--input samplesheet.csv \\  
       \--episegmix\_mode duration \\  
       \--merge true \\  
       \--states 8 \\  
       \--genome hg38 \\  
       \--binsize 200 \\  
       \--decoding\_algorithm viterbi \\  
       \--dist\_histone NBI \\  
       \--dist\_methyl BI \\  
       \--outdir ./results\_duration

6. **Multiple States Run:** Evaluate multiple models in a single run by passing a comma-separated list:  
   nextflow run main.nf \-profile docker \\  
       \--input samplesheet.csv \\  
       \--states 8,9,10 \\  
       \--outdir ./results\_multistate

## **Documentation**

### **Input Samplesheet**

Provide a comma-separated samplesheet.csv with these exact columns:

| sample\_id | replicate | epigenetic\_mark | file\_name | modality | paired\_end | distribution |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| sample\_01 | 1 | H3K27ac | data/file.bam | ChIP-seq | true | NBI |
| sample\_02 | 1 | WGBS | data/file.bed | WGBS | false | BI |

**Path Recommendations:**

* **Recommended:** Use **relative paths** (e.g., data/file.bam) relative to your project directory for better portability and pipeline mounting.  
* **Alternative:** You may use **absolute paths** if your data is stored in fixed external storage locations.

**Notes:**

* ChIP-seq/ATAC-seq require .bam files.  
* WGBS/NOMe-seq require .bed or .bed.gz files.  
* Use the distribution column to specify a statistical model (e.g., NBI, BI). Unsure which to use? Run the fitting module to auto-generate this column.

### **Input File Formats & Pre-processing**

EpiSegMix minimizes manual pre-processing:

* **BAM Files:** No need to sort or index prior to running. The pipeline automatically cleans alignments, standardizes chromosome prefixes, sorts, and indexes (.bai) the data.  
* **BED Files:** Must follow a standard 11-column format containing methylation percentages and coverage metrics:  
  1  10468  10469  100.00  2  \+  10468  10469  210,0,0  0  100.00

### **Analysis Modules (--episegmix\_mode)**

1. **standard**: Default mode. Standard transition probabilities without explicit duration constraints.  
2. **duration**: Topology-based HMM modeling chromatin state lengths for improved spatial resolution.  
3. **DNA**: Segmentation utilizing purely DNA methylation data.  
4. **fitting**: Identifies best-fitting statistical distributions and outputs an updated fitted\_samplesheet.csv containing the optimal models for your data.

### **Key Parameters**

| Parameter | Default | Description |
| :---- | :---- | :---- |
| **Core Inputs & Outputs** |  |  |
| \--input | ./samplesheet.csv | Path to the samplesheet. |
| \--outdir | ./results | Output directory. |
| **Model Configuration** |  |  |
| \--episegmix\_mode | standard | Module selection (standard, duration, DNA, or fitting). |
| \--states | 8 | Number of states (single int or comma-separated list, e.g., 8,9,10). |
| \--genome | hg38 | Reference genome assembly. |
| \--binsize | 200 | Genome bin size in base pairs (bp). |
| \--merge | false | Merge replicates/modalities before training. |
| **Algorithmic Settings** |  |  |
| \--decoding\_algorithm | viterbi | State decoding method (viterbi or posterior). |
| \--chr\_parameter\_estimation | pilot\_hg38 | Chromosome(s) used for parameter estimation. |
| \--iter | 200 | Max iterations for Baum-Welch training. |
| \--epsilon | 1 | Convergence threshold. |
| \--adjustment | 2 | Parameter estimation adjustment factor. |
| **Statistical Distributions** |  |  |
| \--dist\_histone | NBI | Histone count distribution (e.g., PO, BI, NBI). |
| \--dist\_methyl | BI | Methylation count distribution (BI or BB). |

### **Advanced Execution**

**HPC Clusters (Slurm):** EpiSegMix is highly scalable and includes a built-in slurm profile. When running on an HPC, it is highly recommended to use **Singularity**:  
nextflow run main.nf \-profile singularity,slurm \--input samplesheet.csv \--outdir ./results

### **Outputs**

Results are saved in your defined \--outdir:

* EpiSegMix/{id}/Plots/{id}.html: **Main Summary Report** (interactive emission matrices/state distributions).  
* EpiSegMix/{id}/Segmentation/{id}.bed.gz: Final segmentation file.  
* EpiSegMix/{id}/Models/{id}.model.json: Trained HMM parameters.  
* EpiSegMix/Fitting/fitted\_samplesheet.csv: Generated *only* in fitting mode.

## **Support & Contributions**

For bugs, questions, or feature requests, please open an issue on the [EpiSegMix GitHub Issues page](https://www.google.com/search?q=https://github.com/epigenetics-sb/EpiSegMix/issues).  
**Core Developer:** Aaryan Jaitly (@aaryanjaitly)

Incase of queries please contact Nihit Aggarwal: nihit.aggarwal@uni-saarland.de 

## **Citations**

If you use EpiSegMix for your research, please cite:

1. Schmitz, J. E., et al. (2023). EpiSegMix: a flexible distribution hidden Markov model with duration modeling for chromatin state discovery. *Bioinformatics*. doi: [10.1093/bioinformatics/btae178](https://doi.org/10.1093/bioinformatics/btae178)  
2. Aggarwal, N., et al. (2025). Integrated flexible DNA methylation-chromatin segmentation modeling enhances epigenomic state annotation. *bioRxiv*. doi: [10.1101/2025.07.25.666820](https://doi.org/10.1101/2025.07.25.666820)

