# Genomic Landscape and Preclinical Therapeutic Insights Developed from Xenograft, Organoid, and Cell Line Models in the NCI’s Patient-Derived Models Repository

# NGS PIPELINE (PDMR LANDSCAPE PROJECT)
====================================

# DESCRIPTION
-----------
This directory contains scripts and configuration to process next‑generation sequencing (NGS) data for the PDMR Landscape project. 
The goal is to produce reproducible alignment, quality control (QC), variant calls/annotations, and optional copy‑number/coverage
summaries suitable for downstream analyses and figure generation.

# SCOPE
-----
## Inputs:
  - Raw sequencing reads (typically paired-end FASTQ files)
  - Required references (genome fasta, indices, and any caller resources)
  - Optional capture/target BED files for WES/panel datasets

## Outputs:
  - Aligned, coordinate‑sorted BAM/CRAM with index and alignment metrics
  - Variant calls (VCF/MAF as applicable)
  - Copy‑number and/or coverage tables (if enabled)
  - QC reports (e.g., FastQC/MultiQC, coverage summaries)

# FOLDER LAYOUT
-------------
- ngs_pipeline/  : wrappers, configs, and helper scripts to run the pipeline on HPC/Linux
- DATA/          : input tables and intermediate data products (not all files are tracked in Git)
- Scripts/       : figure‑generation and analysis scripts
- Figures/       : rendered figures (PDF/PNG)

# PREREQUISITES
-------------
- Linux environment (HPC or workstation)
- Standard bioinformatics toolchain installed and on PATH (aligner, samtools/bcftools, variant caller, QC tools)
- Sufficient memory/storage for alignment and variant calling
- For HPC (e.g., NIH Biowulf/Helix): run interactive work or submissions from a login/interactive node; submit long jobs via Slurm.
  Compute nodes may restrict outbound network access—stage references/containers in advance.

# REFERENCES & METADATA
---------------------
- Human reference genome (e.g., GRCh38/hg38) plus indices for your aligner(s)
- Known sites/resources as required by chosen callers (e.g., dbSNP, Mills/1000G)
- Optional: target BED file for capture assays (WES/panel)
- Use consistent sample identifiers and a clear filename convention across all outputs. Keep a sample sheet or manifest for batch runs.

# QUICK START
-----------
1) Prepare inputs
   - Place FASTQs in your working area and verify read pairing (_R1/_R2).
   - Ensure references and indices are readable on shared storage.

2) Dry run / help
   - Check the header of the driver script(s) in this folder for expected arguments and environment setup.
   - Many scripts support a "-h" or "--help" option.

3) Run locally (small test) or submit to Slurm (HPC)
   Local example:
     bash run_pipeline.sh -m samples.tsv -o results --threads 8

   Slurm example:
     sbatch --cpus-per-task=16 --mem=64g --time=24:00:00 \
       --wrap="bash run_pipeline.sh -m samples.tsv -o results --threads 16"

   (Replace the command and flags with those used by your specific wrapper or workflow manager.)

# INPUTS AND OUTPUTS
------------------
Expected paths (examples; actual names may vary by tool/config):

  TYPE     PATH EXAMPLE                              NOTES
  ----     ---------------------------------------   -----------------------------------------------
  FASTQ    data/<sample>_R1.fastq.gz                 Paired-end reads expected as _R1/_R2
  BAM      results/align/<sample>.bam                Coordinate-sorted with .bai and alignment metrics
  VCF/MAF  results/variants/<sample>.vcf.gz          Per-sample calls; joint/multi-sample if configured
  CNV      results/cnv/<sample>.*                    Copy-number outputs (tool-dependent)
  QC       results/qc/*                              FastQC/MultiQC, coverage, contamination etc.

# CONFIGURATION
-------------
- Tool paths, reference locations, and run options are defined in config files or at the top of the driver scripts.
- Common toggles include:
  * Aligner choice (e.g., BWA‑MEM2/Novoalign) and read group handling
  * Variant calling mode (somatic vs germline; tumor/normal pairing)
  * Targeted vs whole‑exome/genome settings (BED vs whole genome)
  * QC thresholds (minimum coverage, mapping quality, on‑target rate)

# REPRODUCIBILITY
---------------
- Log tool versions and exact parameters for every run.
- Prefer environment modules or containers (e.g., Singularity/Apptainer, Docker where appropriate).
- Keep a copy of configuration files used for published results.
- Ensure sample metadata (IDs, assay, platform) are consistent across all outputs.

# TROUBLESHOOTING
---------------
- "File not found": confirm the working directory and update input/reference paths in configs.
- Permissions: verify read access to references and write access to output directories.
- Network access on compute nodes: some clusters restrict outbound connectivity; perform downloads/pulls on login/transfer nodes.
- Slow/memory‑heavy steps: adjust thread/memory settings or split batches.

# CITATION AND LICENSE
--------------------
- Project repository: FNL‑MoCha / PDMR_landscape_paper_public
- License: MIT (see the repository top‑level LICENSE file).

CONTACT
-------
For questions or issues, open a GitHub issue in the repository or contact the maintainers above.