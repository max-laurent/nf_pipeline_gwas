# GWAS Pipeline

This repository contains a Nextflow pipeline and R scripts to perform GWAS (Genome-Wide Association Study) using **FASTLMM** and process the results. The pipeline includes formatting inputs, integrating genotype and phenotype data, and extracting significant markers. This pipelines has been designed to receive maize B73v2 genotype data and converts it to B73v4, it had been adapted to run on other organisms yet.

---

## 📁 Repository Structure

- `main.nf` — Nextflow pipeline coordinating all tasks.
- `bin/`
    - `matrix_formating.R` — Formats phenotype and genotype data into a matrix for FASTLMM.
    - `split_pheno_matrix.sh` — Bash script to split phenotypic matrices into separate traits.
    - `extract_signif_marker.R` — Extracts significant GWAS markers (LogPvalue > 5).
    - `formating_fastlmm_output.R` — Formats FASTLMM GWAS outputs, adds statistics, MAF filtering, allele information, and genome version update.
- `docker/`
  - `Dockerfile` — Docker environment.
- Create a `data/` directory locally that contains the data required to run the pipeline, make sure to update the paths in the main.nf file

---

## 🧪 Workflow Overview

1. **Prepare phenotype and genotype data** with `matrix_formating.R`.
2. **Split the phenotype matrice per column into single phenotype files** using `split_pheno_matrix.sh`.
3. **Run FASTLMM with the genotype, phenotype and kinship data** in a bash command
3. **Format GWAS results** using `formating_fastlmm_output.R`.
4. **Extract significant markers** with `extract_signif_marker.R`.
---

## 🚀 Getting Started

### Prerequisites

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)

### Install dependencies

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Build the Docker image from the docker directory
cd docker
docker build -t image_fastlmm_plink:latest .
```

### Run the pipeline

```bash
nextflow run main.nf
```
