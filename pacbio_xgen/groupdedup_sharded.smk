"""Sharded isoseq groupdedup for all tissues in the hifi_locations CSV.

Plan → split (single pass) → N×groupdedup (parallel) → merge → pbindex.
Drop-in replacement for the main pipeline's `rule isoseq_groupdedup`:
  - Inputs:  isoseq/subsample/{tissue}/{tissue}.fltnc.corrected.sorted.subsampled.bam
             isoseq/bcstats/{tissue}/{tissue}.bcstats.tsv
  - Outputs: isoseq/dedup/{tissue}/{tissue}.dedup.bam
             isoseq/dedup/{tissue}/{tissue}.dedup.fasta

The final outputs land at the standard pipeline path so downstream rules
(pbmm2_align, isoseq_collapse, ...) can be triggered by the main snakefile
afterwards. WARNING: if a non-sharded groupdedup job for any of these tissues
is still running, scancel it first to avoid both processes racing on the same
output paths.

Config overrides on the command line:
  --config hifi_locations_csv=…/hifi_locations_xgen.csv n_shards=20

Intermediate per-shard files live under isoseq/dedup_shards/{tissue}/.
"""
configfile: "config.yaml"
import csv
import os

WORKDIR = config["work_directory"]
N_SHARDS = int(config.get("n_shards", 20))
SHARD_IDS = [f"{i:03d}" for i in range(N_SHARDS)]

SHARDING_SCRIPTS = os.path.join(
    os.path.dirname(workflow.snakefile), "sharding_scripts"
)

def _load_tissues():
    csv_path = config["hifi_locations_csv"]
    tissues = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            t = (row.get("tissue") or "").strip()
            if t:
                tissues.append(t)
    if not tissues:
        raise ValueError(f"No tissues found in {csv_path}")
    return sorted(set(tissues))


TISSUES = _load_tissues()


DEFAULT_THREADS = config.get("threads", 32)
DEFAULT_MEMORY = config.get("default_memory", "64G")
DEFAULT_PARTITION = config.get("partition", "owners")
DEFAULT_TIME = config.get("time", "2-0:00:00")
DEFAULT_SLURM_EXTRA = config.get("slurm_extra", "")


def get_resource(rule_name, resource_type, default_value):
    if rule_name in config:
        return config[rule_name].get(resource_type, default_value)
    return default_value


rule all:
    input:
        expand(os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.bam"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.bam.pbi"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.fasta"), tissue=TISSUES),


rule plan_shards:
    input:
        bcstats = os.path.join(WORKDIR, "isoseq/bcstats/{tissue}/{tissue}.bcstats.tsv"),
    output:
        cbs = expand(
            os.path.join(WORKDIR, "isoseq/dedup_shards/{{tissue}}/plan/{shard}.cbs.txt"),
            shard=SHARD_IDS,
        ),
    log:
        os.path.join(WORKDIR, "logs/plan_dedup_shards/{tissue}.log"),
    params:
        plan_dir = lambda wc: os.path.join(WORKDIR, f"isoseq/dedup_shards/{wc.tissue}/plan"),
        n_shards = N_SHARDS,
        max_reads = config.get("max_reads_per_barcode", 0),
        script = os.path.join(SHARDING_SCRIPTS, "plan_dedup_shards.py"),
    threads: get_resource("plan_shards", "threads", 1)
    resources:
        mem = get_resource("plan_shards", "mem", "4G"),
        partition = get_resource("plan_shards", "partition", DEFAULT_PARTITION),
        time = get_resource("plan_shards", "time", "0:30:00"),
        slurm_extra = get_resource("plan_shards", "slurm_extra", DEFAULT_SLURM_EXTRA),
    shell:
        """
        python {params.script} \
            --bcstats {input.bcstats} \
            --n-shards {params.n_shards} \
            --max-reads-per-barcode {params.max_reads} \
            --out-dir {params.plan_dir} 2> {log}
        """


rule split_bam_shards:
    input:
        bam = os.path.join(WORKDIR, "isoseq/subsample/{tissue}/{tissue}.fltnc.corrected.sorted.subsampled.bam"),
        cbs = expand(
            os.path.join(WORKDIR, "isoseq/dedup_shards/{{tissue}}/plan/{shard}.cbs.txt"),
            shard=SHARD_IDS,
        ),
    output:
        bams = expand(
            os.path.join(WORKDIR, "isoseq/dedup_shards/{{tissue}}/bam/{shard}.bam"),
            shard=SHARD_IDS,
        ),
        pbis = expand(
            os.path.join(WORKDIR, "isoseq/dedup_shards/{{tissue}}/bam/{shard}.bam.pbi"),
            shard=SHARD_IDS,
        ),
    log:
        os.path.join(WORKDIR, "logs/split_bam_shards/{tissue}.log"),
    params:
        plan_dir = lambda wc: os.path.join(WORKDIR, f"isoseq/dedup_shards/{wc.tissue}/plan"),
        bam_dir = lambda wc: os.path.join(WORKDIR, f"isoseq/dedup_shards/{wc.tissue}/bam"),
        n_shards = N_SHARDS,
        script = os.path.join(SHARDING_SCRIPTS, "split_bam_by_cb_shards.py"),
    threads: get_resource("split_bam_shards", "threads", 8)
    resources:
        mem = get_resource("split_bam_shards", "mem", "32G"),
        partition = get_resource("split_bam_shards", "partition", DEFAULT_PARTITION),
        time = get_resource("split_bam_shards", "time", "6:00:00"),
        slurm_extra = get_resource("split_bam_shards", "slurm_extra", DEFAULT_SLURM_EXTRA),
    shell:
        """
        python {params.script} \
            --input-bam {input.bam} \
            --shard-cb-dir {params.plan_dir} \
            --out-dir {params.bam_dir} \
            --n-shards {params.n_shards} \
            --threads {threads} 2> {log}
        """


rule isoseq_groupdedup_shard:
    input:
        bam = os.path.join(WORKDIR, "isoseq/dedup_shards/{tissue}/bam/{shard}.bam"),
        pbi = os.path.join(WORKDIR, "isoseq/dedup_shards/{tissue}/bam/{shard}.bam.pbi"),
    output:
        bam = os.path.join(WORKDIR, "isoseq/dedup_shards/{tissue}/dedup/{shard}.dedup.bam"),
        fasta = os.path.join(WORKDIR, "isoseq/dedup_shards/{tissue}/dedup/{shard}.dedup.fasta"),
    log:
        os.path.join(WORKDIR, "logs/isoseq_groupdedup_shard/{tissue}.{shard}.log"),
    threads: get_resource("isoseq_groupdedup_shard", "threads", 16)
    resources:
        mem = get_resource("isoseq_groupdedup_shard", "mem", "64G"),
        partition = get_resource("isoseq_groupdedup_shard", "partition", DEFAULT_PARTITION),
        time = get_resource("isoseq_groupdedup_shard", "time", "12:00:00"),
        slurm_extra = get_resource("isoseq_groupdedup_shard", "slurm_extra", DEFAULT_SLURM_EXTRA),
    shell:
        """
        isoseq groupdedup {input.bam} {output.bam} \
            --num-threads {threads} 2> {log}
        """


rule merge_dedup_shards:
    input:
        bams = expand(
            os.path.join(WORKDIR, "isoseq/dedup_shards/{{tissue}}/dedup/{shard}.dedup.bam"),
            shard=SHARD_IDS,
        ),
        fastas = expand(
            os.path.join(WORKDIR, "isoseq/dedup_shards/{{tissue}}/dedup/{shard}.dedup.fasta"),
            shard=SHARD_IDS,
        ),
    output:
        bam = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.bam"),
        pbi = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.bam.pbi"),
        fasta = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.fasta"),
    log:
        os.path.join(WORKDIR, "logs/merge_dedup_shards/{tissue}.log"),
    params:
        rename_dir = lambda wc: os.path.join(WORKDIR, f"isoseq/dedup_shards/{wc.tissue}/dedup_renamed"),
    threads: get_resource("merge_dedup_shards", "threads", 8)
    resources:
        mem = get_resource("merge_dedup_shards", "mem", "16G"),
        partition = get_resource("merge_dedup_shards", "partition", DEFAULT_PARTITION),
        time = get_resource("merge_dedup_shards", "time", "2:00:00"),
        slurm_extra = get_resource("merge_dedup_shards", "slurm_extra", DEFAULT_SLURM_EXTRA),
    shell:
        """
        # `isoseq groupdedup` re-numbers every output read as `molecule/N` starting
        # at N=0 in each shard AND tags every shard with the same `@RG PU:molecule`,
        # so pbmerge sees identical (movie, ZMW) pairs across all 100 shards and
        # fatally errors with "cannot sort CCS/transcripts that share both movie
        # name & ZMW hole number". pbmerge reads the movie name from `@RG PU:`,
        # so we rewrite both the read names AND the PU field to be shard-unique.
        set -euo pipefail
        mkdir -p {params.rename_dir}
        : > {output.fasta}
        renamed_bams=""
        for bam in {input.bams}; do
            shard=$(basename "$bam" .dedup.bam)
            renamed="{params.rename_dir}/${{shard}}.renamed.bam"
            samtools view -h "$bam" \
                | sed -E "/^@RG/ s|PU:molecule\\b|PU:molecule_${{shard}}|; /^@/!s|^molecule/|molecule_${{shard}}/|" \
                | samtools view -bS - > "$renamed"
            renamed_bams="$renamed_bams $renamed"
        done
        pbmerge -j {threads} -o {output.bam} $renamed_bams 2> {log}
        pbindex {output.bam} 2>> {log}
        for fa in {input.fastas}; do
            shard=$(basename "$fa" .dedup.fasta)
            sed "s|^>molecule/|>molecule_${{shard}}/|" "$fa" >> {output.fasta}
        done
        rm -rf {params.rename_dir}
        """
