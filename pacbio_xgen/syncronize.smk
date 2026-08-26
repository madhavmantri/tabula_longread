configfile: "config.yaml"

import os
import csv
import shlex
from pathlib import Path

# Get directories from config and list of tissues
for i, item in enumerate(sys.path):
    if item.startswith("/share/software/user/"):
        sys.path.pop(i)
        break

WORKDIR = config["work_directory"]
DATADIR = config["raw_data_directory"]
RECOLLAPSED = config.get("recollapsed_directory", "recollapsed")
SYNCED = config.get("synced_directory", "synced")

# Donor(s) whose tissues are excluded from this branch. Accepts a YAML list
# (exclude_donor: [TSP33.1, TSP31]) or a single/comma-separated string from the
# CLI (--config exclude_donor=TSP33.1,TSP31).
_exclude_donor_cfg = config.get("exclude_donor", "TSP31")
if isinstance(_exclude_donor_cfg, str):
    EXCLUDE_DONORS = [d.strip() for d in _exclude_donor_cfg.split(",") if d.strip()]
else:
    EXCLUDE_DONORS = [str(d).strip() for d in _exclude_donor_cfg if str(d).strip()]

def _load_tissues_from_csv():
    """Tissue list from config['hifi_locations_csv'] (column: tissue),
    matching how isoseq.smk discovers tissues so the sync pipeline runs on
    exactly the tissues the upstream Iso-Seq pipeline produced."""
    csv_path = config["hifi_locations_csv"]
    tissues = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            tissue = (row.get("tissue") or "").strip()
            if tissue and tissue not in tissues:
                tissues.append(tissue)
    return tissues

def _tissues_for_donors(donors):
    """Set of tissue (tube_id) IDs belonging to any of `donors`, read from
    config['sample_metadata_csv'] (defaults to csvs/sample_metadata.csv). For
    the xGen pipeline this points at csvs/xgen_metadata.csv, since the xGen
    samples (X33, T70-T73) are absent from the main sample_metadata.csv."""
    donors = set(donors)
    metadata_csv = config.get(
        "sample_metadata_csv",
        os.path.join(WORKDIR, "../csvs/sample_metadata.csv"),
    )
    with open(metadata_csv) as f:
        return {
            (row.get("tube_id") or "").strip()
            for row in csv.DictReader(f)
            if (row.get("donor") or "").strip() in donors
        }

# All HiFi tissues except those from the excluded donor(s)
EXCLUDE_TISSUES = _tissues_for_donors(EXCLUDE_DONORS)
TISSUES = sorted(t for t in _load_tissues_from_csv() if t not in EXCLUDE_TISSUES)

# Configurable threshold
BARCODE_GROUP_SIZE = config.get("barcode_group_size", 100000)
SATURATION_THRESHOLD = config.get("saturation_threshold", 0.75)

# Default resource parameters
DEFAULT_THREADS = config.get("threads", 32)
DEFAULT_MEMORY = config.get("default_memory", "64G")
DEFAULT_PARTITION = config.get("partition", "owners")
DEFAULT_TIME = config.get("time", "2-0:00:00")
DEFAULT_SLURM_EXTRA = config.get("slurm_extra", "")

def get_resource(rule_name, resource_type, default_value):
    """Get resource value with rule-specific override support from config.yaml"""
    if rule_name in config:
        return config[rule_name].get(resource_type, default_value)
    return default_value

rule all:
    input:
        #expand(os.path.join(WORKDIR, SYNCED + "/seurat/{tissue}/genes_seurat/matrix.mtx"), tissue=TISSUES),
        #expand(os.path.join(WORKDIR, SYNCED + "/seurat/{tissue}/isoforms_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, RECOLLAPSED + "/seurat/{tissue}/genes_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, RECOLLAPSED + "/seurat/{tissue}/isoforms_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, RECOLLAPSED + "/seurat_pbids/{tissue}/genes_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, RECOLLAPSED + "/seurat_pbids/{tissue}/isoforms_seurat/matrix.mtx"), tissue=TISSUES),
        expand(os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/IsoAnnotLite_annotated_transcripts.txt"), tissue=TISSUES)

# =============================================================================
# Sync approach: build a global ID map, then rewrite per-tissue files
#
# Output tree:  synced/
#                 ├── id_map/pacbio_id_map.tsv
#                 ├── collapse/{tissue}/collapse.synced.*
#                 ├── classify/{tissue}/collapse_classification.*
#                 └── seurat/{tissue}/genes_seurat/ & isoforms_seurat/
# =============================================================================

rule sync_pacbio_id_map:
    input:
        gffs = expand(os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.gff"), tissue=TISSUES)
    output:
        mapping = os.path.join(WORKDIR, SYNCED + "/id_map/pacbio_id_map.tsv")
    params:
        script = "scripts/sync_pacbio_ids.py",
        gff_args = lambda wildcards, input: " ".join(
            f"--tissue {shlex.quote(tissue)} --gff {shlex.quote(gff)}"
            for tissue, gff in zip(TISSUES, input.gffs)
        )
    log:
        os.path.join(WORKDIR, "logs/" + SYNCED + "/pacbio_id_map.log")
    threads: get_resource("sync_pacbio_id_map", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("sync_pacbio_id_map", "mem", DEFAULT_MEMORY),
        partition=get_resource("sync_pacbio_id_map", "partition", DEFAULT_PARTITION),
        time=get_resource("sync_pacbio_id_map", "time", DEFAULT_TIME),
        slurm_extra=get_resource("sync_pacbio_id_map", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p $(dirname {output.mapping}) $(dirname {log})
        python3 {params.script} build-map \
            {params.gff_args} \
            --map-out {output.mapping} \
            2> {log}
        """


rule sync_pacbio_outputs:
    input:
        mapping = rules.sync_pacbio_id_map.output.mapping,
        gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.gff"),
        group = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.group.txt"),
        abundance = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.abundance.txt")
    output:
        gff = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.gff"),
        group = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.group.txt"),
        abundance = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.abundance.txt")
    params:
        script = "scripts/sync_pacbio_ids.py"
    log:
        os.path.join(WORKDIR, "logs/" + SYNCED + "/{tissue}_rewrite.log")
    threads: get_resource("sync_pacbio_outputs", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("sync_pacbio_outputs", "mem", DEFAULT_MEMORY),
        partition=get_resource("sync_pacbio_outputs", "partition", DEFAULT_PARTITION),
        time=get_resource("sync_pacbio_outputs", "time", DEFAULT_TIME),
        slurm_extra=get_resource("sync_pacbio_outputs", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p $(dirname {output.gff}) $(dirname {log})
        python3 {params.script} rewrite \
            --mapping {input.mapping} \
            --tissue {wildcards.tissue} \
            --gff {input.gff} \
            --group {input.group} \
            --abundance {input.abundance} \
            --gff-out {output.gff} \
            --group-out {output.group} \
            --abundance-out {output.abundance} \
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

rule pigeon_sort_synced:
   input:
       gff = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.gff")
   output:
       sorted_gff = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.gff.sorted")
   log:
       os.path.join(WORKDIR, "logs/" + SYNCED + "/{tissue}_sort.log")
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


rule pigeon_classify_synced:
    input:
        gff = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.gff.sorted"),
        gtf = config["reference_gtf"],
        genome = config["reference_genome"],
        abundance = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.abundance.txt"),
        gtf_index = config["reference_gtf"] + ".pgi",
        fai = config["reference_genome"] + ".fai",
        cage_peak = config["cage_peak"],
        cage_peak_index = config["cage_peak"] + ".pgi",
        polyA_motif = config["polyA_motif_list"]
    output:
        classified = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/collapse_classification.txt")
    params:
        outdir = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/")
    log:
        os.path.join(WORKDIR, "logs/" + SYNCED + "/{tissue}_classify.log")
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


rule pigeon_filter_synced:
    input:
        classification = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/collapse_classification.txt"),
        sorted_gff = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.gff.sorted")
    output:
        sorted_gff_renamed = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.sorted.gff"),
        filtered_gff = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.sorted.filtered_lite.gff"),
        filtered_classification = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        filtered_junctions = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/collapse_classification.filtered_lite_junctions.txt")
    log:
        os.path.join(WORKDIR, "logs/" + SYNCED + "/{tissue}_filter.log")
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


rule run_IsoAnnotLite_synced:
    input:
        filtered_gff = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.sorted.filtered_lite.gff"),
        filtered_classification = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        filtered_junctions = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/collapse_classification.filtered_lite_junctions.txt")
    output:
        annotated_gff = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/IsoAnnotLite_annotated.gff3"),
        annotated_transcript_gff = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/IsoAnnotLite_annotated_transcripts.txt")
    params:
        outprefix = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/IsoAnnotLite_annotated"),
        isoannotlite_script = config["isoannotlite_script"]
    log:
        os.path.join(WORKDIR, "logs/" + SYNCED + "/{tissue}_IsoAnnotLite.log")
    threads: get_resource("run_IsoAnnotLite", "threads", DEFAULT_THREADS)
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


rule pigeon_make_seurat_synced:
    input:
        classification = os.path.join(WORKDIR, SYNCED + "/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        dedup = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.fasta"),
        group = os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.group.txt")
    output:
        gene_matrix = os.path.join(WORKDIR, SYNCED + "/seurat/{tissue}/genes_seurat/matrix.mtx"),
        isoform_matrix = os.path.join(WORKDIR, SYNCED + "/seurat/{tissue}/isoforms_seurat/matrix.mtx")
    params:
        outdir = os.path.join(WORKDIR, SYNCED + "/seurat/{tissue}"),
        keep_novel = "--keep-novel-genes" if config.get("keep_novel_genes", False) else "",
        keep_ribo_mito = "--keep-ribo-mito-genes" if config.get("keep_ribo_mito_genes", False) else ""
    threads: get_resource("pigeon_make_seurat", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_make_seurat", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_make_seurat", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_make_seurat", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_make_seurat", "slurm_extra", DEFAULT_SLURM_EXTRA)
    log:
        os.path.join(WORKDIR, "logs/" + SYNCED + "/{tissue}_make_seurat.log")
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


rule all_synced:
    input:
        synced_gffs = expand(os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.gff"), tissue=TISSUES),
        synced_groups = expand(os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.group.txt"), tissue=TISSUES),
        synced_abundances = expand(os.path.join(WORKDIR, SYNCED + "/collapse/{tissue}/collapse.synced.abundance.txt"), tissue=TISSUES),
        seurat_gene_outputs = expand(os.path.join(WORKDIR, SYNCED + "/seurat/{tissue}/genes_seurat/matrix.mtx"), tissue=TISSUES),
        seurat_isoform_outputs = expand(os.path.join(WORKDIR, SYNCED + "/seurat/{tissue}/isoforms_seurat/matrix.mtx"), tissue=TISSUES)

# =============================================================================
# Recollapse approach: prefix fastas, pool, re-align, re-collapse, rewrite
#
# Output tree:  recollapsed/
#                 ├── prefix/{tissue}/collapse.prefixed.fasta
#                 ├── pooled/collapse.prefixed.{fasta,aligned.bam,gff,group,abundance}
#                 ├── collapse/{tissue}/collapse.recollapsed.*
#                 ├── classify/{tissue}/collapse_classification.*
#                 └── seurat/{tissue}/genes_seurat/ & isoforms_seurat/
# =============================================================================

rule prefix_collapsed_fasta:
    input:
        fasta = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.fasta")
    output:
        fasta = os.path.join(WORKDIR, RECOLLAPSED + "/prefix/{tissue}/collapse.prefixed.fasta")
    params:
        script = "scripts/recollapse_from_pooled_fasta.py"
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/{tissue}_prefix_fasta.log")
    threads: get_resource("prefix_collapsed_fasta", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("prefix_collapsed_fasta", "mem", DEFAULT_MEMORY),
        partition=get_resource("prefix_collapsed_fasta", "partition", DEFAULT_PARTITION),
        time=get_resource("prefix_collapsed_fasta", "time", DEFAULT_TIME),
        slurm_extra=get_resource("prefix_collapsed_fasta", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p $(dirname {output.fasta}) $(dirname {log})
        python3 {params.script} prefix-fasta \
            --tissue {wildcards.tissue} \
            --fasta {input.fasta} \
            --fasta-out {output.fasta} \
            2> {log}
        """


rule pool_prefixed_fastas:
    input:
        fastas = expand(os.path.join(WORKDIR, RECOLLAPSED + "/prefix/{tissue}/collapse.prefixed.fasta"), tissue=TISSUES)
    output:
        fasta = os.path.join(WORKDIR, RECOLLAPSED + "/pooled/collapse.prefixed.fasta")
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/pool_prefixed_fastas.log")
    threads: get_resource("pool_prefixed_fastas", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pool_prefixed_fastas", "mem", DEFAULT_MEMORY),
        partition=get_resource("pool_prefixed_fastas", "partition", DEFAULT_PARTITION),
        time=get_resource("pool_prefixed_fastas", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pool_prefixed_fastas", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p $(dirname {output.fasta}) $(dirname {log})
        cat {input.fastas} > {output.fasta} 2> {log}
        """


rule pbmm2_align_prefixed_pool:
    input:
        genome = config["reference_genome"],
        fasta = rules.pool_prefixed_fastas.output.fasta
    output:
        bam = os.path.join(WORKDIR, RECOLLAPSED + "/pooled/collapse.prefixed.aligned.bam")
    threads: get_resource("pbmm2_align", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pbmm2_align", "mem", DEFAULT_MEMORY),
        partition=get_resource("pbmm2_align", "partition", DEFAULT_PARTITION),
        time=get_resource("pbmm2_align", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pbmm2_align", "slurm_extra", DEFAULT_SLURM_EXTRA)
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/pbmm2_align_prefixed_pool.log")
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        pbmm2 align \
            --preset ISOSEQ \
            --sort \
            -j {threads} \
            {input.genome} \
            {input.fasta} \
            {output.bam} \
            2> {log}
        """


rule isoseq_collapse_prefixed_pool:
    input:
        bam = rules.pbmm2_align_prefixed_pool.output.bam
    output:
        gff = os.path.join(WORKDIR, RECOLLAPSED + "/pooled/collapse.gff"),
        group = os.path.join(WORKDIR, RECOLLAPSED + "/pooled/collapse.group.txt"),
        abundance = os.path.join(WORKDIR, RECOLLAPSED + "/pooled/collapse.abundance.txt")
    threads: get_resource("isoseq_collapse", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("isoseq_collapse", "mem", DEFAULT_MEMORY),
        partition=get_resource("isoseq_collapse", "partition", DEFAULT_PARTITION),
        time=get_resource("isoseq_collapse", "time", DEFAULT_TIME),
        slurm_extra=get_resource("isoseq_collapse", "slurm_extra", DEFAULT_SLURM_EXTRA)
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/isoseq_collapse_prefixed_pool.log")
    shell:
        """
        mkdir -p $(dirname {output.gff}) $(dirname {log})
        isoseq collapse \
            {input.bam} \
            {output.gff} \
            -j {threads} \
            --do-not-collapse-extra-5exons \
            --max-5p-diff 1000 \
            --max-3p-diff 1000 \
            2> {log}
        """


rule rewrite_recollapsed_outputs:
    input:
        pooled_group = rules.isoseq_collapse_prefixed_pool.output.group,
        gff = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.gff"),
        group = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.group.txt"),
        abundance = os.path.join(WORKDIR, "isoseq/collapse/{tissue}/collapse.abundance.txt")
    output:
        mapping = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.map.tsv"),
        gff = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.gff"),
        group = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.group.txt"),
        abundance = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.abundance.txt")
    params:
        script = "scripts/recollapse_from_pooled_fasta.py"
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/{tissue}_rewrite.log")
    threads: get_resource("rewrite_recollapsed_outputs", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("rewrite_recollapsed_outputs", "mem", DEFAULT_MEMORY),
        partition=get_resource("rewrite_recollapsed_outputs", "partition", DEFAULT_PARTITION),
        time=get_resource("rewrite_recollapsed_outputs", "time", DEFAULT_TIME),
        slurm_extra=get_resource("rewrite_recollapsed_outputs", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p $(dirname {output.mapping}) $(dirname {log})
        python3 {params.script} rewrite \
            --tissue {wildcards.tissue} \
            --pooled-group {input.pooled_group} \
            --gff {input.gff} \
            --group {input.group} \
            --abundance {input.abundance} \
            --mapping-out {output.mapping} \
            --gff-out {output.gff} \
            --group-out {output.group} \
            --abundance-out {output.abundance} \
            2> {log}
        """


rule pigeon_sort_recollapsed:
   input:
       gff = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.gff")
   output:
       sorted_gff = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.gff.sorted")
   log:
       os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/{tissue}_sort.log")
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


rule pigeon_classify_recollapsed:
    input:
        gff = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.gff.sorted"),
        gtf = config["reference_gtf"],
        genome = config["reference_genome"],
        abundance = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.abundance.txt"),
        gtf_index = config["reference_gtf"] + ".pgi",
        fai = config["reference_genome"] + ".fai",
        cage_peak = config["cage_peak"],
        cage_peak_index = config["cage_peak"] + ".pgi",
        polyA_motif = config["polyA_motif_list"]
    output:
        classified = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/collapse_classification.txt")
    params:
        outdir = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/")
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/{tissue}_classify.log")
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


rule pigeon_filter_recollapsed:
    input:
        classification = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/collapse_classification.txt"),
        sorted_gff = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.gff.sorted")
    output:
        sorted_gff_renamed = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.sorted.gff"),
        filtered_gff = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.sorted.filtered_lite.gff"),
        filtered_classification = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        filtered_junctions = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/collapse_classification.filtered_lite_junctions.txt")
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/{tissue}_filter.log")
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


rule run_IsoAnnotLite_recollapsed:
    input:
        filtered_gff = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.sorted.filtered_lite.gff"),
        filtered_classification = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        filtered_junctions = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/collapse_classification.filtered_lite_junctions.txt")
    output:
        annotated_gff = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/IsoAnnotLite_annotated.gff3"),
        annotated_transcript_gff = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/IsoAnnotLite_annotated_transcripts.txt")
    params:
        outprefix = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/IsoAnnotLite_annotated"),
        isoannotlite_script = config["isoannotlite_script"]
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/{tissue}_IsoAnnotLite.log")
    threads: get_resource("run_IsoAnnotLite", "threads", DEFAULT_THREADS)
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


rule pigeon_make_seurat_recollapsed:
    input:
        classification = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        dedup = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.fasta"),
        group = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.group.txt")
    output:
        gene_matrix = os.path.join(WORKDIR, RECOLLAPSED + "/seurat/{tissue}/genes_seurat/matrix.mtx"),
        isoform_matrix = os.path.join(WORKDIR, RECOLLAPSED + "/seurat/{tissue}/isoforms_seurat/matrix.mtx")
    params:
        outdir = os.path.join(WORKDIR, RECOLLAPSED + "/seurat/{tissue}"),
        keep_novel = "--keep-novel-genes" if config.get("keep_novel_genes", False) else "",
        keep_ribo_mito = "--keep-ribo-mito-genes" if config.get("keep_ribo_mito_genes", False) else ""
    threads: get_resource("pigeon_make_seurat", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_make_seurat", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_make_seurat", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_make_seurat", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_make_seurat", "slurm_extra", DEFAULT_SLURM_EXTRA)
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/{tissue}_make_seurat.log")
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


rule pigeon_make_seurat_recollapsed_with_pbids:
    # Same as pigeon_make_seurat_recollapsed but keeps PBIDs in the matrix (--keep-pbids)
    # and writes to the "seurat_pbids" folder instead of "seurat".
    input:
        classification = os.path.join(WORKDIR, RECOLLAPSED + "/classify/{tissue}/collapse_classification.filtered_lite_classification.txt"),
        dedup = os.path.join(WORKDIR, "isoseq/dedup/{tissue}/{tissue}.dedup.fasta"),
        group = os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.group.txt")
    output:
        gene_matrix = os.path.join(WORKDIR, RECOLLAPSED + "/seurat_pbids/{tissue}/genes_seurat/matrix.mtx"),
        isoform_matrix = os.path.join(WORKDIR, RECOLLAPSED + "/seurat_pbids/{tissue}/isoforms_seurat/matrix.mtx")
    params:
        outdir = os.path.join(WORKDIR, RECOLLAPSED + "/seurat_pbids/{tissue}"),
        keep_novel = "--keep-novel-genes" if config.get("keep_novel_genes", False) else "",
        keep_ribo_mito = "--keep-ribo-mito-genes" if config.get("keep_ribo_mito_genes", False) else ""
    threads: get_resource("pigeon_make_seurat", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("pigeon_make_seurat", "mem", DEFAULT_MEMORY),
        partition=get_resource("pigeon_make_seurat", "partition", DEFAULT_PARTITION),
        time=get_resource("pigeon_make_seurat", "time", DEFAULT_TIME),
        slurm_extra=get_resource("pigeon_make_seurat", "slurm_extra", DEFAULT_SLURM_EXTRA)
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED + "/{tissue}_make_seurat_pbids.log")
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


rule all_recollapsed:
    input:
        pooled_alignment = rules.pbmm2_align_prefixed_pool.output.bam,
        pooled_group = rules.isoseq_collapse_prefixed_pool.output.group,
        recollapsed_gffs = expand(os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.gff"), tissue=TISSUES),
        recollapsed_groups = expand(os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.group.txt"), tissue=TISSUES),
        recollapsed_abundances = expand(os.path.join(WORKDIR, RECOLLAPSED + "/collapse/{tissue}/collapse.recollapsed.abundance.txt"), tissue=TISSUES),
        seurat_gene_outputs = expand(os.path.join(WORKDIR, RECOLLAPSED + "/seurat/{tissue}/genes_seurat/matrix.mtx"), tissue=TISSUES),
        seurat_isoform_outputs = expand(os.path.join(WORKDIR, RECOLLAPSED + "/seurat/{tissue}/isoforms_seurat/matrix.mtx"), tissue=TISSUES),
        seurat_pbids_gene_outputs = expand(os.path.join(WORKDIR, RECOLLAPSED + "/seurat_pbids/{tissue}/genes_seurat/matrix.mtx"), tissue=TISSUES),
        seurat_pbids_isoform_outputs = expand(os.path.join(WORKDIR, RECOLLAPSED + "/seurat_pbids/{tissue}/isoforms_seurat/matrix.mtx"), tissue=TISSUES)
