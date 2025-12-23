# Isolate Genome Assembly Pipeline

This pipeline performs **genome assembly of microbial isolates** from paired-end Illumina sequencing data using a **reproducible, HPC-friendly workflow** based on **Snakemake**, **Conda**, and **Apptainer/Singularity**.

It is designed for shared computing environments (no root required) and supports multiple assemblers with clear separation of environments and containerized execution.

---

## ✨ Features

- End-to-end **isolate genome assembly**
- Fully reproducible with Conda + container images
- HPC-safe (no Docker, no root privileges)
- Modular, configurable, and easy to extend
- Clean output organization with sentinel `done` files

---

## 🧬 Supported Assemblers

- **metaSPAdes**
  A widely used short-read assembler for microbial genomes.

- **MetaHipMer2 (MHM2)**
  A high-performance, scalable assembler optimized for HPC systems, executed via an Apptainer/Singularity container.

---

## 🔄 Workflow Overview

For each isolate sample, the pipeline performs:

1. **Quality control & trimming**
   - Adapter trimming and quality filtering using `fastp`
   - Generates HTML and JSON QC reports

2. **Genome assembly**
   - Assembly with **metaSPAdes**
   - Assembly with **MHM2** (containerized)

3. **Completion tracking**
   - Each assembler produces a `done` file to mark successful completion

---


## Software to install before running the pipeline

- snakemake>=9.10.0
- apptainer>=1.3.2


## 📥 Input Requirements

### 1️⃣ Sample list file (TSV)

A tab-separated file with **three columns**:

```text
sample_id    /path/to/read1.fastq.gz    /path/to/read2.fastq.gz
```

### 2️⃣ Configuration file (assembly_config.yaml)


## 🧪 Software Environment

Snakemake will automatically create and manage this environment when run with --use-conda.


## ▶️ Running the Pipeline

The command below runs the pipeline using **40 CPU cores**. The `/vol` filesystem is **bind-mounted into the Apptainer container** and used for input and output paths.
Please adjust the number of cores and bind mounts according to your local or HPC environment.

```shell
snakemake -s assembly.smk \
  --use-conda --use-apptainer \
  -j 40 \
  --apptainer-args "--bind /vol:/vol"
```

