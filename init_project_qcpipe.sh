#!/usr/bin/env bash
set -euo pipefail

echo "📌 Nombre del proyecto:"
read -r PROJECT

if [ -z "${PROJECT}" ]; then
  echo "❌ El nombre del proyecto no puede estar vacío"
  exit 1
fi

echo "📁 Creando estructura para el proyecto: ${PROJECT}"

# Directorios
mkdir -p "${PROJECT}/data/fastq"
mkdir -p "${PROJECT}/workflow"

mkdir -p "${PROJECT}/results/fastqc/pre"
mkdir -p "${PROJECT}/results/fastqc/post"
mkdir -p "${PROJECT}/results/fastp"
mkdir -p "${PROJECT}/results/clean"
mkdir -p "${PROJECT}/results/multiqc/pre"
mkdir -p "${PROJECT}/results/multiqc/post"

# README
cat > "${PROJECT}/README.md" << 'EOF'
# Pipeline FASTQ (PE) — FastQC / MultiQC / fastp (auto _1/_2 o _R1/_R2)

## Inputs soportados (paired-end)
- sample_1.fastq.gz / sample_2.fastq.gz
- sample_R1.fastq.gz / sample_R2.fastq.gz
- también soporta sufijos: sample_R1_001.fastq.gz / sample_R2_001.fastq.gz

## Outputs fastp
Se guardan con sufijo `_clean`:
- sample_1_clean.fastq.gz / sample_2_clean.fastq.gz
o
- sample_R1_clean.fastq.gz / sample_R2_clean.fastq.gz

## Uso
1) Copia tus FASTQ a `data/fastq/`
2) Instala herramientas:
   conda create -n qcpipe -c conda-forge -c bioconda snakemake fastqc multiqc fastp -y
   conda activate qcpipe
3) Ejecuta:
   snakemake -j 8 -s workflow/Snakefile --scheduler greedy
EOF

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

# Snakefile (AUTO, PE, outputs con sufijo _clean)
cat > "${PROJECT}/workflow/Snakefile" << 'EOF'
import os, glob

configfile: "workflow/config.yaml"

READS_DIR = config["reads_dir"]
OUTDIR = config["outdir"]
PE_NAMING = config.get("pe_naming", "auto")  # auto | R1 | 1

def detect_mode():
    r1 = glob.glob(os.path.join(READS_DIR, "*_R1*.fastq.gz"))
    u1 = glob.glob(os.path.join(READS_DIR, "*_1.fastq.gz"))
    if r1:
        return "R1"
    if u1:
        return "1"
    raise ValueError(f"No encuentro PE con *_R1*.fastq.gz ni *_1.fastq.gz en {READS_DIR}")

MODE = detect_mode() if PE_NAMING == "auto" else PE_NAMING
if MODE not in {"R1", "1"}:
    raise ValueError("pe_naming debe ser: auto | R1 | 1")

def fq_path(sample, mate):
    if MODE == "1":
        return os.path.join(READS_DIR, f"{sample}_{mate}.fastq.gz")

    pattern = os.path.join(READS_DIR, f"{sample}_R{mate}*.fastq.gz")
    hits = sorted(glob.glob(pattern))
    if not hits:
        raise ValueError(f"No encontré {pattern}")
    if len(hits) > 1:
        print(f"[WARN] {sample} mate {mate}: múltiples matches, usando: {hits[0]}")
    return hits[0]

def list_samples():
    samples = set()
    if MODE == "1":
        for r1 in glob.glob(os.path.join(READS_DIR, "*_1.fastq.gz")):
            base = os.path.basename(r1)
            sample = base[:-len("_1.fastq.gz")]
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
        raise ValueError(f"No se detectaron pares PE válidos en {READS_DIR} (modo {MODE})")
    return samples

SAMPLES = list_samples()
print(f"[INFO] PE mode detectado: {MODE}. Muestras: {len(SAMPLES)}")

def out_prefix(sample, mate):
    return f"{sample}_{mate}" if MODE == "1" else f"{sample}_R{mate}"

def clean_path(sample, mate):
    # salida con sufijo _clean
    if MODE == "1":
        return os.path.join(OUTDIR, "clean", f"{sample}_{mate}_clean.fastq.gz")
    return os.path.join(OUTDIR, "clean", f"{sample}_R{mate}_clean.fastq.gz")

FASTQC_PRE_ZIPS = expand(os.path.join(OUTDIR, "fastqc", "pre", "{p}_fastqc.zip"),
                         p=[out_prefix(s,1) for s in SAMPLES] + [out_prefix(s,2) for s in SAMPLES])

FASTQC_POST_ZIPS = expand(os.path.join(OUTDIR, "fastqc", "post", "{p}_fastqc.zip"),
                          p=[f"{out_prefix(s,1)}_clean" for s in SAMPLES] + [f"{out_prefix(s,2)}_clean" for s in SAMPLES])

FASTP_JSON = expand(os.path.join(OUTDIR, "fastp", "{sample}.json"), sample=SAMPLES)

rule all:
    input:
        os.path.join(OUTDIR, "multiqc", "pre", "multiqc_report.html"),
        os.path.join(OUTDIR, "multiqc", "post", "multiqc_report.html"),
        expand(clean_path("{sample}", 1), sample=SAMPLES),
        expand(clean_path("{sample}", 2), sample=SAMPLES)

rule fastqc_pre:
    input:
        r1=lambda wc: fq_path(wc.sample, 1),
        r2=lambda wc: fq_path(wc.sample, 2)
    output:
        html1=os.path.join(OUTDIR, "fastqc", "pre", lambda wc: f"{out_prefix(wc.sample,1)}_fastqc.html"),
        zip1=os.path.join(OUTDIR, "fastqc", "pre", lambda wc: f"{out_prefix(wc.sample,1)}_fastqc.zip"),
        html2=os.path.join(OUTDIR, "fastqc", "pre", lambda wc: f"{out_prefix(wc.sample,2)}_fastqc.html"),
        zip2=os.path.join(OUTDIR, "fastqc", "pre", lambda wc: f"{out_prefix(wc.sample,2)}_fastqc.zip")
    threads: 2
    shell:
        r"""
        mkdir -p {OUTDIR}/fastqc/pre
        fastqc -t {threads} -o {OUTDIR}/fastqc/pre {input.r1} {input.r2}
        """

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

rule fastp:
    input:
        r1=lambda wc: fq_path(wc.sample, 1),
        r2=lambda wc: fq_path(wc.sample, 2)
    output:
        o1=lambda wc: clean_path(wc.sample, 1),
        o2=lambda wc: clean_path(wc.sample, 2),
        html=os.path.join(OUTDIR, "fastp", "{sample}.html"),
        json=os.path.join(OUTDIR, "fastp", "{sample}.json")
    threads: config["fastp"]["thread"]
    shell:
        r"""
        mkdir -p {OUTDIR}/clean {OUTDIR}/fastp
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.o1} -O {output.o2} \
          --qualified_quality_phred {config[fastp][qualified_quality_phred]} \
          --length_required {config[fastp][length_required]} \
          --detect_adapter_for_pe \
          -w {threads} \
          -h {output.html} -j {output.json}
        """

rule fastqc_post:
    input:
        r1=lambda wc: clean_path(wc.sample, 1),
        r2=lambda wc: clean_path(wc.sample, 2)
    output:
        # FastQC nombra outputs según el nombre del archivo de entrada, por eso usamos *_clean aquí
        html1=os.path.join(OUTDIR, "fastqc", "post", lambda wc: f"{out_prefix(wc.sample,1)}_clean_fastqc.html"),
        zip1=os.path.join(OUTDIR, "fastqc", "post", lambda wc: f"{out_prefix(wc.sample,1)}_clean_fastqc.zip"),
        html2=os.path.join(OUTDIR, "fastqc", "post", lambda wc: f"{out_prefix(wc.sample,2)}_clean_fastqc.html"),
        zip2=os.path.join(OUTDIR, "fastqc", "post", lambda wc: f"{out_prefix(wc.sample,2)}_clean_fastqc.zip")
    threads: 2
    shell:
        r"""
        mkdir -p {OUTDIR}/fastqc/post
        fastqc -t {threads} -o {OUTDIR}/fastqc/post {input.r1} {input.r2}
        """

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

echo "✅ Proyecto '${PROJECT}' creado (auto _1/_2 o _R1/_R2) + fastp outputs *_clean.fastq.gz"
echo "➡️ Siguiente:"
echo "   1) Copia tus FASTQ a: ${PROJECT}/data/fastq/"
echo "   2) Entra: cd ${PROJECT}"
echo "   3) Ejecuta: snakemake -j 8 -s workflow/Snakefile --scheduler greedy"

