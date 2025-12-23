from os import path

# Use workflow.source_path to find config relative to the Snakefile location
configfile: workflow.source_path("assembly_config.yaml")

# ---- Setup Paths & Variables ----
dataset   = config["dataset"]
# Ensure output dir is absolute to avoid confusion
BASE_OUT  = path.abspath(config["out_dir"])
proj_dir  = path.join(BASE_OUT, dataset)

threads_default = int(config.get("threads", 20))

out_data_dir = path.join(proj_dir, "data")
qc_dir       = path.join(out_data_dir, "qc_fq")
logs_dir     = path.join(proj_dir, "logs")
results_dir  = path.join(proj_dir, "results")

# ---- Parse Samples ----
fq_file_dict = {}
sample_list_path = config["sample_list_file"]

if not path.exists(sample_list_path):
    raise FileNotFoundError(f"Sample list not found at: {sample_list_path}")

with open(sample_list_path, "r") as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 3:
            raise ValueError(f"Invalid line (need 3 cols): {line!r}")
        # Store as [R1, R2]
        fq_file_dict[cols[0]] = cols[1:3]

samples = list(fq_file_dict.keys())

# ---- Helpers & Constraints ----
def get_fq(wc):
    return fq_file_dict[wc.sample]

wildcard_constraints:
    assembler="[a-zA-Z0-9]+",
    sample="[A-Za-z-_0-9]+"

assemblers = ["metaspades", "mhm2"]

# ---- Env & Config Loading ----
ASSEMBLY_ENV = config.get("assembly_env", "assembly_env.yaml")
FASTP_EXTRA  = config.get("fastp", {}).get("extra", "")

# ---- Rules ----

rule all:
    input:
        expand(results_dir + "/{assembler}/{sample}/done", sample=samples, assembler=assemblers)

rule fastp:
    input:
        reads = get_fq
    output:
        or1  = path.join(qc_dir, "{sample}.qc.r1.fq.gz"),
        or2  = path.join(qc_dir, "{sample}.qc.r2.fq.gz"),
        html = path.join(qc_dir, "{sample}.qc.report.html"),
        json = path.join(qc_dir, "{sample}.qc.report.json")
    threads: 8
    conda:
        ASSEMBLY_ENV
    log:
        path.join(logs_dir, "fastp/{sample}.log")
    shell:
        r"""
        fastp -i {input.reads[0]} -I {input.reads[1]} \
              -o {output.or1} -O {output.or2} \
              {FASTP_EXTRA} \
              -h {output.html} -j {output.json} \
              -w {threads} \
              2> {log}
        """

rule metaspades:
    input:
        r1 = rules.fastp.output.or1,
        r2 = rules.fastp.output.or2,
    output:
        done = touch(path.join(results_dir, "metaspades/{sample}/done")),
    threads: threads_default
    resources:
        mem_mb = 250000 # Request 250GB RAM (adjust as needed)
    conda:
        ASSEMBLY_ENV
    log:
        path.join(logs_dir, "metaspades/{sample}.log")
    params:
        outdir = path.join(results_dir, "metaspades/{sample}"),
    shell:
        r"""
        # Spades fails if dir exists and is not empty. Clean start.
        if [ -d "{params.outdir}" ]; then
            rm -rf "{params.outdir}"
        fi
        
        metaspades.py -k auto \
                      -1 {input.r1} -2 {input.r2} \
                      -o {params.outdir} \
                      -t {threads} > {log} 2>&1
        """

rule ungz:
    """
    Decompresses QC reads for MHM2.
    Marked as temp() to save disk space after MHM2 is done.
    """
    input:
        r1 = rules.fastp.output.or1,
        r2 = rules.fastp.output.or2,
    output:
        ungz_r1 = path.join(qc_dir, "{sample}.qc.r1.fq"),
        ungz_r2 = path.join(qc_dir, "{sample}.qc.r2.fq"),
    threads: 10
    conda:
        ASSEMBLY_ENV
    shell:
        r"""
        # Explicitly stream to output file for safety
        pigz -dc -p {threads} {input.r1} > {output.ungz_r1}
        pigz -dc -p {threads} {input.r2} > {output.ungz_r2}
        """

rule mhm2:
    input:
        r1 = rules.ungz.output.ungz_r1,
        r2 = rules.ungz.output.ungz_r2,
    output:
        done = touch(path.join(results_dir, "mhm2/{sample}/done")),
    threads: threads_default
    resources:
        mem_mb = 250000 
    log:
        path.join(logs_dir, "mhm2/{sample}.log")
    params:
        outdir = path.join(results_dir, "mhm2/{sample}"),
        extra  = config.get("mhm2", {}).get("extra_mhm2", "").strip(),
        procs  = lambda wc: int(config.get("mhm2", {}).get("ranks", threads_default)),
        nodes  = lambda wc: int(config.get("mhm2", {}).get("nodes", 1)),
    singularity:
        config["mhm2_image"]
    shell:
        """
        mhm2.py --procs {params.procs} --nodes {params.nodes} {params.extra} \
                -p {input.r1} {input.r2} \
                -o {params.outdir} > {log} 2>&1
        """
