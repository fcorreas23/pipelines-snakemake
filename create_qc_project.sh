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

configfile: "workflow/config.yaml"

# Main pipeline configuration.
READS_DIR = config["reads_dir"]
OUTDIR = config["outdir"]
PE_NAMING = config.get("pe_naming", "auto")  # auto | R1 | 1


# Automatically detects the PE read naming scheme:
# - R1/R2: sample_R1*.fastq.gz and sample_R2*.fastq.gz
# - 1/2:   sample_1.fastq.gz and sample_2.fastq.gz
def detect_mode():
    r1_files = glob.glob(os.path.join(READS_DIR, "*_R1*.fastq.gz"))
    one_files = glob.glob(os.path.join(READS_DIR, "*_1.fastq.gz"))

    if r1_files:
        return "R1"
    if one_files:
        return "1"

    raise ValueError(
        f"No PE files found with *_R1*.fastq.gz or *_1.fastq.gz in {READS_DIR}"
    )


# Final naming mode defined by config or auto-detection.
MODE = detect_mode() if PE_NAMING == "auto" else PE_NAMING
if MODE not in {"R1", "1"}:
    raise ValueError("pe_naming must be: auto | R1 | 1")


# Returns FASTQ path for a sample/mate based on active mode.
def fq_path(sample, mate):
    if MODE == "1":
        return os.path.join(READS_DIR, f"{sample}_{mate}.fastq.gz")

    pattern = os.path.join(READS_DIR, f"{sample}_R{mate}*.fastq.gz")
    hits = sorted(glob.glob(pattern))
    if not hits:
        raise ValueError(f"Could not find {pattern}")
    if len(hits) > 1:
        print(f"[WARN] {sample} mate {mate}: multiple matches, using: {hits[0]}")
    return hits[0]


# Lists valid samples with complete PE pairs (R1/R2 or 1/2).
def list_samples():
    samples = set()

    if MODE == "1":
        for r1 in glob.glob(os.path.join(READS_DIR, "*_1.fastq.gz")):
            base = os.path.basename(r1)
            sample = base[: -len("_1.fastq.gz")]
            r2 = os.path.join(READS_DIR, f"{sample}_2.fastq.gz")
            if os.path.exists(r2):
                samples.add(sample)
    else:
        for r1 in glob.glob(os.path.join(READS_DIR, "*_R1*.fastq.gz")):
            base = os.path.basename(r1)
            sample = base.split("_R1")[0]
            r2 = os.path.join(READS_DIR, base.replace("_R1", "_R2", 1))
            if os.path.exists(r2):
                samples.add(sample)

    samples = sorted(samples)
    if not samples:
        raise ValueError(
            f"No valid PE pairs detected in {READS_DIR} (mode {MODE})"
        )
    return samples


# Builds output prefix per sample/mate for FastQC reuse.
def out_prefix(sample, mate):
    return f"{sample}_{mate}" if MODE == "1" else f"{sample}_R{mate}"


# Path to cleaned FASTQ produced by fastp (_clean suffix).
def clean_path(sample, mate):
    if MODE == "1":
        return os.path.join(OUTDIR, "clean", f"{sample}_{mate}_clean.fastq.gz")
    return os.path.join(OUTDIR, "clean", f"{sample}_R{mate}_clean.fastq.gz")


SAMPLES = list_samples()
print(f"[INFO] Detected PE mode: {MODE}. Samples: {len(SAMPLES)}")

FASTQC_PRE_ZIPS = expand(
    os.path.join(OUTDIR, "fastqc", "pre", "{p}_fastqc.zip"),
    p=[out_prefix(s, 1) for s in SAMPLES] + [out_prefix(s, 2) for s in SAMPLES],
)

FASTQC_POST_ZIPS = expand(
    os.path.join(OUTDIR, "fastqc", "post", "{p}_fastqc.zip"),
    p=[f"{out_prefix(s, 1)}_clean" for s in SAMPLES]
    + [f"{out_prefix(s, 2)}_clean" for s in SAMPLES],
)

FASTP_JSON = expand(os.path.join(OUTDIR, "fastp", "{sample}.json"), sample=SAMPLES)


# Target rule: defines all expected final outputs.
rule all:
    input:
        os.path.join(OUTDIR, "multiqc", "pre", "multiqc_report.html"),
        os.path.join(OUTDIR, "multiqc", "post", "multiqc_report.html"),
        expand(clean_path("{sample}", 1), sample=SAMPLES),
        expand(clean_path("{sample}", 2), sample=SAMPLES)


# Runs FastQC on raw reads (R1 and R2) per sample.
rule fastqc_pre:
    input:
        r1=lambda wc: fq_path(wc.sample, 1),
        r2=lambda wc: fq_path(wc.sample, 2)
    output:
        html1=os.path.join(OUTDIR, "fastqc", "pre", f"{out_prefix('{sample}', 1)}_fastqc.html"),
        zip1=os.path.join(OUTDIR, "fastqc", "pre", f"{out_prefix('{sample}', 1)}_fastqc.zip"),
        html2=os.path.join(OUTDIR, "fastqc", "pre", f"{out_prefix('{sample}', 2)}_fastqc.html"),
        zip2=os.path.join(OUTDIR, "fastqc", "pre", f"{out_prefix('{sample}', 2)}_fastqc.zip")
    threads: 2
    shell:
        r"""
        mkdir -p {OUTDIR}/fastqc/pre
        fastqc -t {threads} -o {OUTDIR}/fastqc/pre {input.r1} {input.r2}
        """


# Consolidates pre-cleaning FastQC reports into MultiQC.
rule multiqc_pre:
    input:
        FASTQC_PRE_ZIPS
    output:
        report=os.path.join(OUTDIR, "multiqc", "pre", "multiqc_report.html")
    shell:
        r"""
        mkdir -p {OUTDIR}/multiqc/pre
        multiqc -o {OUTDIR}/multiqc/pre {OUTDIR}/fastqc/pre
        """


# Runs PE trimming/filtering with fastp and generates HTML/JSON reports.
rule fastp:
    input:
        r1=lambda wc: fq_path(wc.sample, 1),
        r2=lambda wc: fq_path(wc.sample, 2)
    output:
        o1=clean_path("{sample}", 1),
        o2=clean_path("{sample}", 2),
        html=os.path.join(OUTDIR, "fastp", "{sample}.html"),
        json=os.path.join(OUTDIR, "fastp", "{sample}.json")
    params:
        q=config["fastp"]["qualified_quality_phred"],
        min_len=config["fastp"]["length_required"]
    threads: config["fastp"]["thread"]
    shell:
        r"""
        mkdir -p {OUTDIR}/clean {OUTDIR}/fastp
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.o1} -O {output.o2} \
          --qualified_quality_phred {params.q} \
          --length_required {params.min_len} \
          --detect_adapter_for_pe \
          -w {threads} \
          -h {output.html} -j {output.json}
        """


# Runs FastQC on cleaned reads produced by fastp.
rule fastqc_post:
    input:
        r1=lambda wc: clean_path(wc.sample, 1),
        r2=lambda wc: clean_path(wc.sample, 2)
    output:
        # FastQC uses the input filename as base, so the _clean suffix appears.
        html1=os.path.join(OUTDIR, "fastqc", "post", f"{out_prefix('{sample}', 1)}_clean_fastqc.html"),
        zip1=os.path.join(OUTDIR, "fastqc", "post", f"{out_prefix('{sample}', 1)}_clean_fastqc.zip"),
        html2=os.path.join(OUTDIR, "fastqc", "post", f"{out_prefix('{sample}', 2)}_clean_fastqc.html"),
        zip2=os.path.join(OUTDIR, "fastqc", "post", f"{out_prefix('{sample}', 2)}_clean_fastqc.zip")
    threads: 2
    shell:
        r"""
        mkdir -p {OUTDIR}/fastqc/post
        fastqc -t {threads} -o {OUTDIR}/fastqc/post {input.r1} {input.r2}
        """


# Consolidates post-cleaning FastQC + fastp metrics into final MultiQC.
rule multiqc_post:
    input:
        FASTQC_POST_ZIPS + FASTP_JSON
    output:
        report=os.path.join(OUTDIR, "multiqc", "post", "multiqc_report.html")
    shell:
        r"""
        mkdir -p {OUTDIR}/multiqc/post
        multiqc -o {OUTDIR}/multiqc/post {OUTDIR}/fastqc/post {OUTDIR}/fastp
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



