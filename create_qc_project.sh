#!/usr/bin/env bash
set -euo pipefail

echo "Project name: "
read -r PROJECT

if [ -z "${PROJECT}" ]; then
    echo "❌ Project name cannot be empty"
    exit 1
fi

echo "📁 Creating project structure for: ${PROJECT}"

# Directories
mkdir -p "${PROJECT}/data/fastq"
mkdir -p "${PROJECT}/workflow"
mkdir -p "${PROJECT}/results/fastqc/pre"
mkdir -p "${PROJECT}/results/fastqc/post"
mkdir -p "${PROJECT}/results/fastp"
mkdir -p "${PROJECT}/results/clean"
mkdir -p "${PROJECT}/results/multiqc/pre"
mkdir -p "${PROJECT}/results/multiqc/post"

echo "✅ Project structure '${PROJECT}' created successfully"

echo "📄 Creating configuration and workflow files..."

# config.yaml
cat > "${PROJECT}/workflow/config.yaml" << 'EOF'
reads_dir: "data/fastq"
outdir: "results"

pe_naming: "auto"   # auto | R1 | 1

fastp:
    qualified_quality_phred: 20
    length_required: 50
    thread: 8
EOF

# Snakefile
cat > "${PROJECT}/workflow/Snakefile" << 'EOF'
# Snakefile for the quality control pipeline

import glob
import os
import re

configfile: "workflow/config.yaml"

# ---------------------------------------------------------------------
# Main configuration
# ---------------------------------------------------------------------

READS_DIR = config["reads_dir"]
OUTDIR = config["outdir"]

# auto | R1 | 1
PE_NAMING = config.get("pe_naming", "auto")

# Valid FASTQ extensions
FASTQ_EXTENSIONS = ["fastq.gz", "fq.gz"]


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def find_files(patterns):

    files = []

    for p in patterns:
        files.extend(glob.glob(p))

    return sorted(files)


# Automatically detects PE naming:
#
# R1/R2:
#   sample_R1.fastq.gz
#   sample_R2.fastq.gz
#
# 1/2:
#   sample_1.fastq.gz
#   sample_2.fastq.gz
#
def detect_mode():

    r1_patterns = [
        os.path.join(
            READS_DIR,
            f"*_R1*.{ext}"
        )
        for ext in FASTQ_EXTENSIONS
    ]

    one_patterns = [
        os.path.join(
            READS_DIR,
            f"*_1.{ext}"
        )
        for ext in FASTQ_EXTENSIONS
    ]

    r1_files = find_files(r1_patterns)
    one_files = find_files(one_patterns)

    if r1_files:
        return "R1"

    if one_files:
        return "1"

    raise ValueError(
        f"No PE files detected in {READS_DIR}"
    )


# Final mode
MODE = detect_mode() if PE_NAMING == "auto" else PE_NAMING

if MODE not in {"R1", "1"}:
    raise ValueError(
        "pe_naming must be: auto | R1 | 1"
    )


# Returns FASTQ path for sample/mate
def fq_path(sample, mate):

    if MODE == "1":

        patterns = [
            os.path.join(
                READS_DIR,
                f"{sample}_{mate}.{ext}"
            )
            for ext in FASTQ_EXTENSIONS
        ]

    else:

        patterns = [
            os.path.join(
                READS_DIR,
                f"{sample}_R{mate}*.{ext}"
            )
            for ext in FASTQ_EXTENSIONS
        ]

    hits = find_files(patterns)

    if not hits:
        raise ValueError(
            f"Could not find FASTQ for {sample} mate {mate}"
        )

    if len(hits) > 1:
        print(
            f"[WARN] {sample} mate {mate}: "
            f"multiple matches found, using {hits[0]}"
        )

    return hits[0]


# Detects valid paired-end samples
def list_samples():

    samples = set()

    if MODE == "R1":

        for ext in FASTQ_EXTENSIONS:

            files = glob.glob(
                os.path.join(
                    READS_DIR,
                    f"*_R1*.{ext}"
                )
            )

            for f in files:

                base = os.path.basename(f)

                sample = re.sub(
                    r"_R1.*",
                    "",
                    base
                )

                r2_exists = False

                for ext2 in FASTQ_EXTENSIONS:

                    r2_pattern = os.path.join(
                        READS_DIR,
                        f"{sample}_R2*.{ext2}"
                    )

                    if glob.glob(r2_pattern):
                        r2_exists = True
                        break

                if r2_exists:
                    samples.add(sample)

    else:

        for ext in FASTQ_EXTENSIONS:

            files = glob.glob(
                os.path.join(
                    READS_DIR,
                    f"*_1.{ext}"
                )
            )

            for f in files:

                base = os.path.basename(f)

                sample = re.sub(
                    r"_1.*",
                    "",
                    base
                )

                r2 = os.path.join(
                    READS_DIR,
                    f"{sample}_2.{ext}"
                )

                if os.path.exists(r2):
                    samples.add(sample)

    if not samples:

        raise ValueError(
            f"No valid PE pairs detected "
            f"in {READS_DIR} (mode {MODE})"
        )

    return sorted(samples)


# Output naming helper
def out_prefix(sample, mate):

    if MODE == "1":
        return f"{sample}_{mate}"

    return f"{sample}_R{mate}"


# Clean FASTQ path
def clean_path(sample, mate):

    if MODE == "1":

        return os.path.join(
            OUTDIR,
            "clean",
            f"{sample}_{mate}_clean.fastq.gz"
        )

    return os.path.join(
        OUTDIR,
        "clean",
        f"{sample}_R{mate}_clean.fastq.gz"
    )


# ---------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------

SAMPLES = list_samples()

print(
    f"[INFO] Detected PE mode: {MODE}. "
    f"Samples detected: {len(SAMPLES)}"
)


# ---------------------------------------------------------------------
# Global outputs
# ---------------------------------------------------------------------

fastqc_pre_outputs = expand(
    os.path.join(
        OUTDIR,
        "fastqc",
        "pre",
        "{sample}_R{mate}_fastqc.zip"
    ),
    sample=SAMPLES,
    mate=[1, 2]
)

fastqc_post_outputs = expand(
    os.path.join(
        OUTDIR,
        "fastqc",
        "post",
        "{sample}_R{mate}_clean_fastqc.zip"
    ),
    sample=SAMPLES,
    mate=[1, 2]
)

FASTP_JSON = expand(
    os.path.join(
        OUTDIR,
        "fastp",
        "{sample}.json"
    ),
    sample=SAMPLES
)


# ---------------------------------------------------------------------
# Final targets
# ---------------------------------------------------------------------

rule all:
    input:
        os.path.join(
            OUTDIR,
            "multiqc",
            "pre",
            "multiqc_report.html"
        ),
        os.path.join(
            OUTDIR,
            "multiqc",
            "post",
            "multiqc_report.html"
        ),
        expand(
            clean_path("{sample}", 1),
            sample=SAMPLES
        ),
        expand(
            clean_path("{sample}", 2),
            sample=SAMPLES
        )


# ---------------------------------------------------------------------
# FASTQC BEFORE CLEANING
# ---------------------------------------------------------------------

rule fastqc_pre:
    input:
        r1=lambda wc: fq_path(wc.sample, 1),
        r2=lambda wc: fq_path(wc.sample, 2)

    output:
        html1=os.path.join(
            OUTDIR,
            "fastqc",
            "pre",
            "{sample}_R1_fastqc.html"
        ),
        zip1=os.path.join(
            OUTDIR,
            "fastqc",
            "pre",
            "{sample}_R1_fastqc.zip"
        ),
        html2=os.path.join(
            OUTDIR,
            "fastqc",
            "pre",
            "{sample}_R2_fastqc.html"
        ),
        zip2=os.path.join(
            OUTDIR,
            "fastqc",
            "pre",
            "{sample}_R2_fastqc.zip"
        )

    threads: 2

    shell:
        r"""
        mkdir -p {OUTDIR}/fastqc/pre

        fastqc \
            -t {threads} \
            -o {OUTDIR}/fastqc/pre \
            {input.r1} \
            {input.r2}
        """


# ---------------------------------------------------------------------
# MULTIQC BEFORE CLEANING
# ---------------------------------------------------------------------

rule multiqc_pre:
    input:
        fastqc_pre_outputs

    output:
        html=os.path.join(
            OUTDIR,
            "multiqc",
            "pre",
            "multiqc_report.html"
        )

    shell:
        r"""
        mkdir -p {OUTDIR}/multiqc/pre

        multiqc \
            {OUTDIR}/fastqc/pre \
            -o {OUTDIR}/multiqc/pre
        """


# ---------------------------------------------------------------------
# FASTP
# ---------------------------------------------------------------------

rule fastp:
    input:
        r1=lambda wc: fq_path(wc.sample, 1),
        r2=lambda wc: fq_path(wc.sample, 2)

    output:
        o1=clean_path("{sample}", 1),
        o2=clean_path("{sample}", 2),
        html=os.path.join(
            OUTDIR,
            "fastp",
            "{sample}.html"
        ),
        json=os.path.join(
            OUTDIR,
            "fastp",
            "{sample}.json"
        )

    params:
        q=config["fastp"]["qualified_quality_phred"],
        min_len=config["fastp"]["length_required"]

    threads: config["fastp"]["thread"]

    shell:
        r"""
        mkdir -p {OUTDIR}/clean
        mkdir -p {OUTDIR}/fastp

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.o1} \
            -O {output.o2} \
            --qualified_quality_phred {params.q} \
            --length_required {params.min_len} \
            --detect_adapter_for_pe \
            -w {threads} \
            -h {output.html} \
            -j {output.json}
        """


# ---------------------------------------------------------------------
# FASTQC AFTER CLEANING
# ---------------------------------------------------------------------

rule fastqc_post:
    input:
        r1=lambda wc: clean_path(wc.sample, 1),
        r2=lambda wc: clean_path(wc.sample, 2)

    output:
        html1=os.path.join(
            OUTDIR,
            "fastqc",
            "post",
            "{sample}_R1_clean_fastqc.html"
        ),
        zip1=os.path.join(
            OUTDIR,
            "fastqc",
            "post",
            "{sample}_R1_clean_fastqc.zip"
        ),
        html2=os.path.join(
            OUTDIR,
            "fastqc",
            "post",
            "{sample}_R2_clean_fastqc.html"
        ),
        zip2=os.path.join(
            OUTDIR,
            "fastqc",
            "post",
            "{sample}_R2_clean_fastqc.zip"
        )

    threads: 2

    shell:
        r"""
        mkdir -p {OUTDIR}/fastqc/post

        fastqc \
            -t {threads} \
            -o {OUTDIR}/fastqc/post \
            {input.r1} \
            {input.r2}
        """


# ---------------------------------------------------------------------
# MULTIQC AFTER CLEANING
# ---------------------------------------------------------------------

rule multiqc_post:
    input:
        fastqc_post_outputs,
        FASTP_JSON

    output:
        html=os.path.join(
            OUTDIR,
            "multiqc",
            "post",
            "multiqc_report.html"
        )

    shell:
        r"""
        mkdir -p {OUTDIR}/multiqc/post

        multiqc \
            {OUTDIR}/fastqc/post \
            {OUTDIR}/fastp \
            -o {OUTDIR}/multiqc/post
        """
EOF

echo "✅ Configuration and workflow files created successfully"

# Check that conda is available
if ! command -v conda >/dev/null 2>&1; then
    echo "❌ Conda is not installed or not available in PATH"
    echo "   Install Miniconda/Anaconda and run this script again"
    exit 1
fi

# Check whether the environment already exists
if conda env list | awk 'NR > 2 {print $1}' | grep -qx "qc-pipe"; then
    echo "ℹ️ Conda environment 'qc-pipe' already exists. Skipping creation"
else
    echo "📦 Creating Conda environment 'qc-pipe' with required packages..."
    conda create -n qc-pipe -c conda-forge -c bioconda snakemake fastqc multiqc fastp -y
    echo "✅ Conda environment 'qc-pipe' created successfully"
fi

echo ""
echo "🚀 Next steps:"
echo "1) Copy your FASTQ files to: ${PROJECT}/data/fastq/"
echo "2) Activate the environment: conda activate qc-pipe"
echo "3) Run the pipeline: snakemake -j 8 -s workflow/Snakefile --scheduler greedy"



