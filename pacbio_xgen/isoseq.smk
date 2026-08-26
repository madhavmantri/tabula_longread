configfile: "config.yaml"
import os
import sys
from pathlib import Path
import csv
import glob
import math

sys.path = [p for p in sys.path if not p.startswith("/share/software/user/")]

WORKDIR = config["work_directory"]

def _load_tissue_dirs():
    """Build {tissue: hifi_dir} from config['hifi_locations_csv']
    (columns: tissue,hifi_dir). `hifi_dir` must point at the per-tissue
    directory containing the .bam files."""
    csv_path = config["hifi_locations_csv"]
    tissue_dirs = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            tissue = (row.get("tissue") or "").strip()
            hifi_dir = (row.get("hifi_dir") or "").strip()
            if not tissue or not hifi_dir:
                continue
            if tissue in tissue_dirs:
                raise ValueError(
                    f"[isoseq.smk] duplicate tissue {tissue!r} in {csv_path}"
                )
            tissue_dirs[tissue] = hifi_dir
    return tissue_dirs

TISSUE_DIRS = _load_tissue_dirs()
TISSUES = sorted(TISSUE_DIRS.keys())

# Default resource parameters
DEFAULT_THREADS = config.get("threads", 32)
DEFAULT_MEMORY = config.get("default_memory", "64G")
DEFAULT_PARTITION = config.get("partition", "owners")
DEFAULT_TIME = config.get("time", "2-0:00:00")
DEFAULT_SLURM_EXTRA = config.get("slurm_extra", "")

def get_bam_files(tissue):
    """Return sorted list of HiFi .bam files in a tissue's HiFi directory.
    Matches `*.hifi_reads*.bam` to exclude fail_reads BAMs, then drops
    `*.hifi_reads.unassigned.bam` (reads that failed demultiplexing) so
    only barcode-assigned HiFi reads enter the pipeline."""
    bams = glob.glob(os.path.join(TISSUE_DIRS[tissue], "*.hifi_reads*.bam"))
    bams = [b for b in bams if not os.path.basename(b).endswith(".hifi_reads.unassigned.bam")]
    # Multiplexed Signios delivery folders pack several tissues into one directory,
    # embedding the tissue id in each barcoded BAM (e.g. '..._T43.bam'). When such
    # tags are present, keep only this tissue's BAM. Legacy single-tissue dirs whose
    # BAMs carry no tissue tag (e.g. '...bcM0003.bam') match nothing here and fall
    # through to the full list unchanged. '_T4.bam' won't match '_T43.bam' (exact
    # suffix), so adjacent tissue ids never collide.
    tagged = [b for b in bams if os.path.basename(b).endswith(f"_{tissue}.bam")]
    if tagged:
        bams = tagged
    return sorted(bams)

def get_bam_input(wildcards):
    """Map {tissue} + {bam_idx} back to the actual input .bam file."""
    bams = get_bam_files(wildcards.tissue)
    idx = int(wildcards.bam_idx) - 1
    return bams[idx]

def get_all_skera_outputs():
    """Build the full list of expected skera output files."""
    outputs = []
    for tissue in TISSUES:
        bams = get_bam_files(tissue)
        for i in range(1, len(bams) + 1):
            outputs.append(
                os.path.join(WORKDIR, f"skera/{tissue}/segmented-{i}.bam")
            )
    return outputs

def get_resource(rule_name, resource_type, default_value):
    """Get resource value with rule-specific override support from config.yaml"""
    if rule_name in config:
        return config[rule_name].get(resource_type, default_value)
    return default_value

def get_segmented_bams(wildcards):
    """Return expected segmented BAM paths for a tissue based on input BAM count."""
    bams = get_bam_files(wildcards.tissue)
    return [
        os.path.join(WORKDIR, f"skera/{wildcards.tissue}/segmented-{i}.bam")
        for i in range(1, len(bams) + 1)
    ]

rule all:
    input:
        #get_all_skera_outputs(),
        #expand(os.path.join(WORKDIR, "skera/{tissue}/segmented.bam"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "pigeon/seurat/{tissue}/genes_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "pigeon/seurat/{tissue}/isoforms_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "pigeon/seurat_pbids/{tissue}/genes_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "pigeon/seurat_pbids/{tissue}/isoforms_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "isoseq/bcstats/{tissue}/{tissue}.bcstats.tsv"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "pigeon/classify/{tissue}/IsoAnnotLite_annotated_transcripts.txt"), tissue=TISSUES),
        #expand(os.path.join(WORKDIR, "isoseq/dedup_combined/{tissue}/{tissue}.dedup.bam"), tissue=TISSUES),
        #expand(os.path.join(WORKDIR, "isoseq/dedup_combined/{tissue}/{tissue}.dedup.fasta"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "isoseq/saturation/{tissue}/{tissue}.rarefaction.tsv"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, "isoseq/saturation/{tissue}/{tissue}.saturation_rarefaction.pdf"), tissue=TISSUES)

rule skera_split:
    input:
        bam=get_bam_input
    output:
        bam=os.path.join(WORKDIR, "skera/{tissue}/segmented-{bam_idx}.bam")
    params:
        primers=config["primer_fasta"],
    log:
        os.path.join(WORKDIR, "logs/skera/{tissue}_{bam_idx}.log")
    benchmark:
        "benchmarks/skera/{tissue}_{bam_idx}.tsv"
    threads: get_resource("skera_split", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("skera_split", "mem", DEFAULT_MEMORY),
        partition=get_resource("skera_split", "partition", DEFAULT_PARTITION),
        time=get_resource("skera_split", "time", DEFAULT_TIME),
        slurm_extra=get_resource("skera_split", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        skera split \
            --num-threads {threads} \
            {input.bam} \
            {params.primers} \
            {output.bam} \
            > {log} 2>&1
        """

rule merge_segmented:
    input:
        bams = get_segmented_bams
    output:
        merged_bam = os.path.join(WORKDIR, "raw_bams/{tissue}/segmented.merged.bam")
    log:
        os.path.join(WORKDIR, "logs/merge_segmented/{tissue}.log")
    threads: get_resource("merge_segmented", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("merge_segmented", "mem", DEFAULT_MEMORY),
        partition=get_resource("merge_segmented", "partition", DEFAULT_PARTITION),
        time=get_resource("merge_segmented", "time", DEFAULT_TIME),
        slurm_extra=get_resource("merge_segmented", "slurm_extra", DEFAULT_SLURM_EXTRA)
    run:
        if len(input.bams) == 1:
            shell("ln -s {input.bams[0]} {output.merged_bam} 2>> {log}")
        else:
            shell("pbmerge -j {threads} -o {output.merged_bam} {input.bams} 2>> {log}")

rule lima_remove_primers:
    input:
        segmented_bam = os.path.join(WORKDIR, "raw_bams/{tissue}/segmented.merged.bam"),
        primers = os.path.join(WORKDIR, "primers/primers.fasta")
    output: 
        fl_bam = os.path.join(WORKDIR, "lima/{tissue}/{tissue}.fl.5p--3p.bam")
    log:
        os.path.join(WORKDIR, "logs/lima_remove_primers/{tissue}.log")
    params:
        prefix = os.path.join(WORKDIR, "lima/{tissue}/{tissue}.fl.bam")
    threads: get_resource("lima_remove_primers", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("lima_remove_primers", "mem", DEFAULT_MEMORY),
        partition=get_resource("lima_remove_primers", "partition", DEFAULT_PARTITION),
        time=get_resource("lima_remove_primers", "time", DEFAULT_TIME),
        slurm_extra=get_resource("lima_remove_primers", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        lima {input.segmented_bam} {input.primers} {params.prefix} \
            --isoseq \
            -j {threads} 2> {log}
        """

rule isoseq_tag:
    input:
        fl_bam = os.path.join(WORKDIR, "lima/{tissue}/{tissue}.fl.5p--3p.bam")
    output:
        flt_bam = os.path.join(WORKDIR, "isoseq/tag/{tissue}/{tissue}.flt.bam")
    log:
        os.path.join(WORKDIR, "logs/isoseq_tag/{tissue}.log")
    params:
        design = "T-12U-16B"
    threads: get_resource("isoseq_tag", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("isoseq_tag", "mem", DEFAULT_MEMORY),
        partition=get_resource("isoseq_tag", "partition", DEFAULT_PARTITION),
        time=get_resource("isoseq_tag", "time", DEFAULT_TIME),
        slurm_extra=get_resource("isoseq_tag", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        isoseq tag {input.fl_bam} {output.flt_bam} \
            --design {params.design} \
            -j {threads} 2> {log}
        """

rule isoseq_refine:
    input:
        flt_bam = os.path.join(WORKDIR, "isoseq/tag/{tissue}/{tissue}.flt.bam"),
        primers = os.path.join(WORKDIR, "primers/primers.fasta")
    output:
        fltnc_bam = os.path.join(WORKDIR, "isoseq/refine/{tissue}/{tissue}.fltnc.bam")
    log:
        os.path.join(WORKDIR, "logs/isoseq_refine/{tissue}.log")
    threads: get_resource("isoseq_refine", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("isoseq_refine", "mem", DEFAULT_MEMORY),
        partition=get_resource("isoseq_refine", "partition", DEFAULT_PARTITION),
        time=get_resource("isoseq_refine", "time", DEFAULT_TIME),
        slurm_extra=get_resource("isoseq_refine", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        isoseq refine {input.flt_bam} {input.primers} {output.fltnc_bam} \
            --require-polya \
            -j {threads} 2> {log}
        """

rule isoseq_correct:
    input:
        fltnc_bam = os.path.join(WORKDIR, "isoseq/refine/{tissue}/{tissue}.fltnc.bam"),
        barcodes = os.path.join(WORKDIR, "barcodes/3M-february-2018-REVERSE-COMPLEMENTED.txt")
    output:
        corrected_bam = os.path.join(WORKDIR, "isoseq/correct/{tissue}/{tissue}.fltnc.corrected.bam")
    log:
        os.path.join(WORKDIR, "logs/isoseq_correct/{tissue}.log")
    params:
        percentile = config.get("percentile_cutoff", 90)
    threads: get_resource("isoseq_correct", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("isoseq_correct", "mem", DEFAULT_MEMORY),
        partition=get_resource("isoseq_correct", "partition", DEFAULT_PARTITION),
        time=get_resource("isoseq_correct", "time", DEFAULT_TIME),
        slurm_extra=get_resource("isoseq_correct", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        isoseq correct --barcodes {input.barcodes} {input.fltnc_bam} {output.corrected_bam} \
            --method percentile --percentile {params.percentile} -j {threads} 2> {log}
        """

rule samtools_sort_by_barcode:
    input:
        corrected_bam = os.path.join(WORKDIR, "isoseq/correct/{tissue}/{tissue}.fltnc.corrected.bam")
    output:
        sorted_bam = os.path.join(WORKDIR, "isoseq/sort/{tissue}/{tissue}.fltnc.corrected.sorted.bam")
    log:
        os.path.join(WORKDIR, "logs/samtools_sort_by_barcode/{tissue}.log")
    threads: get_resource("samtools_sort_by_barcode", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("samtools_sort_by_barcode", "mem", DEFAULT_MEMORY),
        partition=get_resource("samtools_sort_by_barcode", "partition", DEFAULT_PARTITION),
        time=get_resource("samtools_sort_by_barcode", "time", DEFAULT_TIME),
        slurm_extra=get_resource("samtools_sort_by_barcode", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        samtools sort -t CB -@ {threads} {input.corrected_bam} -o {output.sorted_bam} 2> {log}
        """

rule isoseq_bcstats:
    input:
        sorted_bam = os.path.join(WORKDIR, "isoseq/sort/{tissue}/{tissue}.fltnc.corrected.sorted.bam")
    output:
        json = os.path.join(WORKDIR, "isoseq/bcstats/{tissue}/{tissue}.bcstats.json"),
        tsv = os.path.join(WORKDIR, "isoseq/bcstats/{tissue}/{tissue}.bcstats.tsv")
    log:
        os.path.join(WORKDIR, "logs/isoseq_bcstats/{tissue}.log")
    params:
        percentile = config.get("percentile_cutoff", 90)
    threads: get_resource("isoseq_bcstats", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("isoseq_bcstats", "mem", DEFAULT_MEMORY),
        partition=get_resource("isoseq_bcstats", "partition", DEFAULT_PARTITION),
        time=get_resource("isoseq_bcstats", "time", DEFAULT_TIME),
        slurm_extra=get_resource("isoseq_bcstats", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        isoseq bcstats --json {output.json} -o {output.tsv} {input.sorted_bam} \
            -j {threads} --method percentile --percentile {params.percentile} 2> {log}
        """

rule saturation_rarefaction:
    input:
        sorted_bam = os.path.join(WORKDIR, "isoseq/sort/{tissue}/{tissue}.fltnc.corrected.sorted.bam")
    output:
        rarefaction_tsv = os.path.join(WORKDIR, "isoseq/saturation/{tissue}/{tissue}.rarefaction.tsv"),
        per_cell_tsv = os.path.join(WORKDIR, "isoseq/saturation/{tissue}/{tissue}.per_cell_saturation.tsv"),
        pdf = os.path.join(WORKDIR, "isoseq/saturation/{tissue}/{tissue}.saturation_rarefaction.pdf")
    params:
        out_prefix = os.path.join(WORKDIR, "isoseq/saturation/{tissue}/{tissue}"),
        barcode_tag = config.get("barcode_tag", "CB")
    log:
        os.path.join(WORKDIR, "logs/saturation_rarefaction/{tissue}.log")
    threads: get_resource("saturation_rarefaction", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("saturation_rarefaction", "mem", DEFAULT_MEMORY),
        partition=get_resource("saturation_rarefaction", "partition", DEFAULT_PARTITION),
        time=get_resource("saturation_rarefaction", "time", DEFAULT_TIME),
        slurm_extra=get_resource("saturation_rarefaction", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        python scripts/saturation_rarefaction.py \
            --bam {input.sorted_bam} \
            --out-prefix {params.out_prefix} \
            --barcode-tag {params.barcode_tag} \
            --threads {threads} \
            2> {log}
        """

rule index_sorted_bam:
    input:
        sorted_bam = os.path.join(WORKDIR, "isoseq/sort/{tissue}/{tissue}.fltnc.corrected.sorted.bam")
    output:
        bai = os.path.join(WORKDIR, "isoseq/sort/{tissue}/{tissue}.fltnc.corrected.sorted.bam.bai")
    log:
        os.path.join(WORKDIR, "logs/index_bam/{tissue}.log")
    threads: get_resource("index_sorted_bam", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("index_sorted_bam", "mem", DEFAULT_MEMORY),
        partition=get_resource("index_sorted_bam", "partition", DEFAULT_PARTITION),
        time=get_resource("index_sorted_bam", "time", DEFAULT_TIME),
        slurm_extra=get_resource("index_sorted_bam", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        samtools index {input.sorted_bam} -@ {threads} 2> {log}
        """

rule pbindex_sorted_bam:
    input:
        sorted_bam = os.path.join(WORKDIR, "isoseq/sort/{tissue}/{tissue}.fltnc.corrected.sorted.bam")
    output:
        pbi = os.path.join(WORKDIR, "isoseq/sort/{tissue}/{tissue}.fltnc.corrected.sorted.bam.pbi")
    log:
        os.path.join(WORKDIR, "logs/pbindex_sorted_bam/{tissue}.log")
    threads: get_resource("pbindex_sorted_bam", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pbindex_sorted_bam", "mem", DEFAULT_MEMORY),
        partition=get_resource("pbindex_sorted_bam", "partition", DEFAULT_PARTITION),
        time=get_resource("pbindex_sorted_bam", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pbindex_sorted_bam", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        pbindex \
            -j {threads} \
            {input.sorted_bam} \
            2> {log}
        """

rule subsample_cell_barcode:
    input:
        sorted_bam = os.path.join(WORKDIR, "isoseq/sort/{tissue}/{tissue}.fltnc.corrected.sorted.bam")
    output:
        subsampled_bam = os.path.join(WORKDIR, "isoseq/subsample/{tissue}/{tissue}.fltnc.corrected.sorted.subsampled.bam"),
        subsampled_pbi = os.path.join(WORKDIR, "isoseq/subsample/{tissue}/{tissue}.fltnc.corrected.sorted.subsampled.bam.pbi")
    params:
        max_reads = config.get("max_reads_per_barcode", 300000),
        tag = config.get("barcode_tag", "CB")
    log:
        os.path.join(WORKDIR, "logs/subsample_cell_barcode/{tissue}.log")
    threads: get_resource("subsample_cell_barcode", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("subsample_cell_barcode", "mem", DEFAULT_MEMORY),
        partition=get_resource("subsample_cell_barcode", "partition", DEFAULT_PARTITION),
        time=get_resource("subsample_cell_barcode", "time", DEFAULT_TIME),
        slurm_extra=get_resource("subsample_cell_barcode", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        python scripts/subsample_cell_barcode.py \
            --bam {input.sorted_bam} \
            --maxreads {params.max_reads} \
            --tag {params.tag} \
            --threads {threads} \
            --output {output.subsampled_bam} \
            2> {log}
        pbindex -j {threads} {output.subsampled_bam} 2>> {log}
        """

rule isoseq_groupdedup:
    input:
        sorted_bam = os.path.join(WORKDIR, "isoseq/subsample/{tissue}/{tissue}.fltnc.corrected.sorted.subsampled.bam"),
        sorted_pbi = os.path.join(WORKDIR, "isoseq/subsample/{tissue}/{tissue}.fltnc.corrected.sorted.subsampled.bam.pbi")
    output:
        dedup_bam = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.bam"),
        dedup_fasta = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.fasta")
    log:
        os.path.join(WORKDIR, "logs/isoseq_groupdedup/{tissue}.log")
    threads: get_resource("isoseq_groupdedup", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("isoseq_groupdedup", "mem", DEFAULT_MEMORY),
        partition=get_resource("isoseq_groupdedup", "partition", DEFAULT_PARTITION),
        time=get_resource("isoseq_groupdedup", "time", DEFAULT_TIME),
        slurm_extra=get_resource("isoseq_groupdedup", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        isoseq groupdedup {input.sorted_bam} {output.dedup_bam} \
         --num-threads {threads} 2> {log}
        """

rule pbmm2_align:
    input:
        dedup_bam = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.bam"),
        reference = config["reference_genome"]  
    output:
        aligned_bam = os.path.join(WORKDIR, "pbmm2/aligned/{tissue}/{tissue}.aligned.bam")
    log:
        os.path.join(WORKDIR, "logs/pbmm2_align/{tissue}.log")
    threads: get_resource("pbmm2_align", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pbmm2_align", "mem", DEFAULT_MEMORY),
        partition=get_resource("pbmm2_align", "partition", DEFAULT_PARTITION),
        time=get_resource("pbmm2_align", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pbmm2_align", "slurm_extra", DEFAULT_SLURM_EXTRA)
    params:
        preset = "ISOSEQ"
    shell:
        """
        pbmm2 align {input.reference} {input.dedup_bam} {output.aligned_bam} \
            --preset {params.preset} \
            --sort \
            -j {threads} 2> {log}
        """

rule pbindex:
    input:
        bam = os.path.join(WORKDIR, "pbmm2/aligned/{tissue}/{tissue}.aligned.bam")
    output:
        pbi = os.path.join(WORKDIR, "pbmm2/aligned/{tissue}/{tissue}.aligned.bam.pbi")
    log:
        os.path.join(WORKDIR, "logs/pbindex/{tissue}.log")
    threads: get_resource("pbindex", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pbindex", "mem", DEFAULT_MEMORY),
        partition=get_resource("pbindex", "partition", DEFAULT_PARTITION),
        time=get_resource("pbindex", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pbindex", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        pbindex \
            -j {threads} \
            {input.bam} \
            2> {log}
        """

rule isoseq_collapse:
    input:
        bam = os.path.join(WORKDIR, "pbmm2/aligned/{tissue}/{tissue}.aligned.bam"),
    output:
        gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.gff"),
        group = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.group.txt"),
        abundance = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.abundance.txt")        
    log:
        os.path.join(WORKDIR, "logs/isoseq_collapse/{tissue}.log")
    threads: get_resource("isoseq_collapse", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("isoseq_collapse", "mem", DEFAULT_MEMORY),
        partition=get_resource("isoseq_collapse", "partition", DEFAULT_PARTITION),
        time=get_resource("isoseq_collapse", "time", DEFAULT_TIME),
        slurm_extra=get_resource("isoseq_collapse", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p $(dirname {output.gff})
        isoseq collapse \
            {input.bam} \
            {output.gff} \
            -j {threads} \
            --do-not-collapse-extra-5exons \
            --max-5p-diff 1000 \
            --max-3p-diff 1000 \
            2> {log}
        """

rule pigeon_index_gtf:
    input:
        gtf = config["reference_gtf"]
    output:
        gtf_index = config["reference_gtf"] + ".pgi"
    log:
        os.path.join(WORKDIR, "logs/pigeon_index_gtf/index_gtf.log")
    threads: get_resource("pigeon_index_gtf", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_index_gtf", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_index_gtf", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_index_gtf", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_index_gtf", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        pigeon index \
            {input.gtf} \
            -j {threads} \
            2> {log}
        """

rule pigeon_index_cage_peak:
    input:
        cage_peak = config["cage_peak"]
    output:
        cage_peak_index = config["cage_peak"] + ".pgi"
    log:
        os.path.join(WORKDIR, "logs/pigeon_index_cage_peak/index_cage_peak.log")
    resources:
        mem=get_resource("pigeon_index_cage_peak", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_index_cage_peak", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_index_cage_peak", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_index_cage_peak", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        pigeon index \
            {input.cage_peak} \
            2> {log}
        """

rule index_fasta:
    input:
        genome = config["reference_genome"]
    output:
        fai = config["reference_genome"] + ".fai"
    log:
        os.path.join(WORKDIR, "logs/index_fasta/faidx.log")
    resources:
        mem=get_resource("index_fasta", "mem", DEFAULT_MEMORY),
        partition=get_resource("index_fasta", "partition", DEFAULT_PARTITION),
        time=get_resource("index_fasta", "time", DEFAULT_TIME),
        slurm_extra=get_resource("index_fasta", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        samtools faidx {input.genome} 2> {log}
        """

rule pigeon_sort:
    input:
        gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.gff")
    output:
        sorted_gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.gff.sorted")
    log:
        os.path.join(WORKDIR, "logs/pigeon_sort/{tissue}.log")
    threads: get_resource("pigeon_sort", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_sort", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_sort", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_sort", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_sort", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        pigeon sort \
            {input.gff} \
            {output.sorted_gff} \
            2> {log}
        """

rule pigeon_classify:    
    input:
        gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.gff.sorted"),
        gtf = config["reference_gtf"],
        genome = config["reference_genome"],
        abundance = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.abundance.txt"),
        gtf_index = config["reference_gtf"] + ".pgi",
        fai = config["reference_genome"] + ".fai",
        cage_peak = config["cage_peak"],
        cage_peak_index = config["cage_peak"] + ".pgi",
        polyA_motif = config["polyA_motif_list"]
    output:
        classified = os.path.join(WORKDIR, "pigeon/classify/{tissue}/collapse_classification.txt")
    params:
        outdir = os.path.join(WORKDIR, "pigeon/classify/{tissue}/")
    log:
        os.path.join(WORKDIR, "logs/pigeon_classify/{tissue}.log")
    threads: get_resource("pigeon_classify", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_classify", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_classify", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_classify", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_classify", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p {params.outdir}
        pigeon classify \
            {input.gff} \
            {input.gtf} \
            {input.genome} \
            --out-dir {params.outdir} \
            --num-threads {threads} \
            --fl {input.abundance} \
            --cage-peak {input.cage_peak} \
            --poly-a {input.polyA_motif} \
            2> {log}
        """

rule pigeon_filter:
    input:
        classification = os.path.join(WORKDIR, "pigeon/classify/{tissue}/collapse_classification.txt"),
        sorted_gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.gff.sorted")
    output:
        sorted_gff_renamed = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.sorted.gff"),
        filtered_gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.sorted.filtered_lite.gff"),
        filtered_classification = os.path.join(WORKDIR, "pigeon/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        filtered_junctions = os.path.join(WORKDIR, "pigeon/classify/{tissue}/collapse_classification.filtered_lite_junctions.txt")
    log:
        os.path.join(WORKDIR, "logs/pigeon_filter/{tissue}.log")
    threads: get_resource("pigeon_filter", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_filter", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_filter", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_filter", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_filter", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        cp {input.sorted_gff} {output.sorted_gff_renamed}
        pigeon filter \
            {input.classification} \
            --isoforms {output.sorted_gff_renamed} \
            -j {threads} \
            2> {log}
        """
        
rule run_IsoAnnotLite:
    input:
        filtered_gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.sorted.filtered_lite.gff"),
        filtered_classification = os.path.join(WORKDIR, "pigeon/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        filtered_junctions = os.path.join(WORKDIR, "pigeon/classify/{tissue}/collapse_classification.filtered_lite_junctions.txt")
    output:
        annotated_gff = os.path.join(WORKDIR, "pigeon/classify/{tissue}/IsoAnnotLite_annotated.gff3"),
        annotated_transcript_gff = os.path.join(WORKDIR, "pigeon/classify/{tissue}/IsoAnnotLite_annotated_transcripts.txt")
    params:
        outprefix = os.path.join(WORKDIR, "pigeon/classify/{tissue}/IsoAnnotLite_annotated"),
        isoannotlite_script = config["isoannotlite_script"]
    log:
        os.path.join(WORKDIR, "logs/run_IsoAnnotLite/{tissue}.log")
    resources:
        mem=get_resource("run_IsoAnnotLite", "mem", DEFAULT_MEMORY),
        partition=get_resource("run_IsoAnnotLite", "partition", DEFAULT_PARTITION),
        time=get_resource("run_IsoAnnotLite", "time", DEFAULT_TIME),
        slurm_extra=get_resource("run_IsoAnnotLite", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        python {params.isoannotlite_script} \
            {input.filtered_gff} \
            {input.filtered_classification} \
            {input.filtered_junctions} \
            -o {params.outprefix} \
            2> {log}
        awk '$3=="transcript" {{
            split($9,a,";");
            for(i in a) {{
                if(a[i]~/^ID=/) {{
                    id=substr(a[i],4);
                    print $1,(id=="novel"?$1:id);
                }}
            }}
        }}' {output.annotated_gff} > {output.annotated_transcript_gff}
        """

rule pigeon_make_seurat:
    input:
        classification = os.path.join(WORKDIR, "pigeon/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        dedup = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.fasta"),
        group = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.group.txt")
    output:
        gene_matrix = os.path.join(WORKDIR, "pigeon/seurat/{tissue}/genes_seurat/matrix.mtx"),
        isoform_matrix = os.path.join(WORKDIR, "pigeon/seurat/{tissue}/isoforms_seurat/matrix.mtx")
    params:
        outdir = os.path.join(WORKDIR, "pigeon/seurat/{tissue}"),
        keep_novel = "--keep-novel-genes" if config.get("keep_novel_genes", False) else "",
        keep_ribo_mito = "--keep-ribo-mito-genes" if config.get("keep_ribo_mito_genes", False) else ""
    threads: get_resource("pigeon_make_seurat", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_make_seurat", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_make_seurat", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_make_seurat", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_make_seurat", "slurm_extra", DEFAULT_SLURM_EXTRA)
    log:
        os.path.join(WORKDIR, "logs/pigeon_make_seurat/{tissue}.log")
    shell:
        """
        mkdir -p {params.outdir}
        pigeon make-seurat \
            {params.keep_novel} \
            {params.keep_ribo_mito} \
            -o {wildcards.tissue} \
            -j {threads} \
            --dedup {input.dedup} \
            -g {input.group} \
            -d {params.outdir} \
            {input.classification} \
            2> {log}
        """


rule pigeon_make_seurat_with_pbids:
    # Same as pigeon_make_seurat but keeps PBIDs in the matrix (--keep-pbids) and writes
    # to the "pigeon/seurat_pbids" folder instead of "pigeon/seurat".
    input:
        classification = os.path.join(WORKDIR, "pigeon/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        dedup = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.fasta"),
        group = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.group.txt")
    output:
        gene_matrix = os.path.join(WORKDIR, "pigeon/seurat_pbids/{tissue}/genes_seurat/matrix.mtx"),
        isoform_matrix = os.path.join(WORKDIR, "pigeon/seurat_pbids/{tissue}/isoforms_seurat/matrix.mtx")
    params:
        outdir = os.path.join(WORKDIR, "pigeon/seurat_pbids/{tissue}"),
        keep_novel = "--keep-novel-genes" if config.get("keep_novel_genes", False) else "",
        keep_ribo_mito = "--keep-ribo-mito-genes" if config.get("keep_ribo_mito_genes", False) else ""
    threads: get_resource("pigeon_make_seurat", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_make_seurat", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_make_seurat", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_make_seurat", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_make_seurat", "slurm_extra", DEFAULT_SLURM_EXTRA)
    log:
        os.path.join(WORKDIR, "logs/pigeon_make_seurat_pbids/{tissue}.log")
    shell:
        """
        mkdir -p {params.outdir}
        pigeon make-seurat \
            {params.keep_novel} \
            {params.keep_ribo_mito} \
            --keep-pbids \
            -o {wildcards.tissue} \
            -j {threads} \
            --dedup {input.dedup} \
            -g {input.group} \
            -d {params.outdir} \
            {input.classification} \
            2> {log}
        """