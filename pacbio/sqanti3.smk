configfile: "config.yaml"

import os
import sys
from pathlib import Path

# Remove problematic shared software paths (same as syncronize.smk)
for i, item in enumerate(sys.path):
    if item.startswith("/share/software/user/"):
        sys.path.pop(i)
        break

WORKDIR = config["work_directory"]
RECOLLAPSED_DIR = config.get("recollapsed_directory", "recollapsed")
RECOLLAPSED_COLLAPSE = os.path.join(WORKDIR, RECOLLAPSED_DIR, "collapse")
SQANTI3_DIR = os.path.join(WORKDIR, RECOLLAPSED_DIR, "sqanti3")
TISSUES = [d.name for d in Path(RECOLLAPSED_COLLAPSE).iterdir() if d.is_dir()]

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
        expand(os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_classification.txt"), tissue=TISSUES),
        expand(os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_RulesFilter_classification.txt"), tissue=TISSUES),
        expand(os.path.join(SQANTI3_DIR, "{tissue}/IsoAnnotLite_annotated_transcripts.txt"), tissue=TISSUES)


rule sqanti3_qc:
    input:
        gff = os.path.join(RECOLLAPSED_COLLAPSE, "{tissue}/collapse.recollapsed.sorted.filtered_lite.gff"),
        gtf = config["reference_gtf"],
        genome = config["reference_genome"],
        abundance = os.path.join(RECOLLAPSED_COLLAPSE, "{tissue}/collapse.recollapsed.abundance.txt")
    output:
        classification = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_classification.txt"),
        junctions = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_junctions.txt"),
        corrected_fasta = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_corrected.fasta"),
        corrected_gtf = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_corrected.gtf"),
        params_file = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}.qc_params.txt")
    params:
        sqanti3_qc_script = config["sqanti3_qc_script"],
        outdir = os.path.join(SQANTI3_DIR, "{tissue}"),
        cage_peak = config["cage_peak"],
        polyA_motif_list = config["polyA_motif_list"],
        # Optional PolyA-site reference (BED). Populates within_polyA_site /
        # dist_to_polyA_site. Empty string (no config['polyA_peak']) = flag omitted.
        polyA_peak_arg = ("--polyA_peak " + config["polyA_peak"]) if config.get("polyA_peak") else ""
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED_DIR + "/sqanti3/{tissue}_qc.log")
    threads: get_resource("sqanti3_qc", "threads", 16)
    resources:
        mem=get_resource("sqanti3_qc", "mem", "64G"),
        partition=get_resource("sqanti3_qc", "partition", DEFAULT_PARTITION),
        time=get_resource("sqanti3_qc", "time", "8:00:00"),
        slurm_extra=get_resource("sqanti3_qc", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        # FL count = deduplicated MOLECULES only (count_fl, col 2). Previously this passed
        # both count_fl AND fl_assoc (col 3, raw reads); SQANTI3 sums numeric columns as
        # separate samples, so FL came out = molecules + reads. Keep col 2 only -> FL = molecules.
        grep -v '^#' {input.abundance} | cut -f1-2 | awk -F'\\t' 'NR==1{{print; next}} {{count[$1]+=$2}} END{{for(k in count) print k"\\t"count[k]}}' > $(dirname {input.abundance})/abundance_nocomments.tsv
        python {params.sqanti3_qc_script} \
            --isoforms {input.gff} \
            --refGTF {input.gtf} \
            --refFasta {input.genome} \
            -o {wildcards.tissue} \
            -d {params.outdir} \
            -t {threads} \
            --fl_count $(dirname {input.abundance})/abundance_nocomments.tsv \
            --CAGE_peak {params.cage_peak} \
            --polyA_motif_list {params.polyA_motif_list} \
            {params.polyA_peak_arg} \
            --include_ORF \
            --report skip \
            2> {log}
        """


rule sqanti3_filter:
    input:
        classification = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_classification.txt"),
        corrected_fasta = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_corrected.fasta"),
        corrected_gtf = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_corrected.gtf")
    output:
        filtered_classification = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_RulesFilter_classification.txt"),
        filtered_fasta = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}.filtered.fasta"),
        filtered_gtf = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}.filtered.gtf")
    params:
        sqanti3_filter_script = config["sqanti3_filter_script"],
        outdir = os.path.join(SQANTI3_DIR, "{tissue}")
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED_DIR + "/sqanti3/{tissue}_filter.log")
    threads: get_resource("sqanti3_filter", "threads", 4)
    resources:
        mem=get_resource("sqanti3_filter", "mem", "32G"),
        partition=get_resource("sqanti3_filter", "partition", DEFAULT_PARTITION),
        time=get_resource("sqanti3_filter", "time", "2:00:00"),
        slurm_extra=get_resource("sqanti3_filter", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        mkdir -p $(dirname {log})
        python {params.sqanti3_filter_script} rules \
            --sqanti_class {input.classification} \
            --filter_isoforms {input.corrected_fasta} \
            --filter_gtf {input.corrected_gtf} \
            -o {wildcards.tissue} \
            -d {params.outdir} \
            2> {log}
        """


rule run_IsoAnnotLite:
    input:
        filtered_gtf = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}.filtered.gtf"),
        filtered_classification = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_RulesFilter_classification.txt"),
        junctions = os.path.join(SQANTI3_DIR, "{tissue}/{tissue}_junctions.txt")
    output:
        annotated_gff = os.path.join(SQANTI3_DIR, "{tissue}/IsoAnnotLite_annotated.gff3"),
        annotated_transcript_gff = os.path.join(SQANTI3_DIR, "{tissue}/IsoAnnotLite_annotated_transcripts.txt")
    params:
        outprefix = os.path.join(SQANTI3_DIR, "{tissue}/IsoAnnotLite_annotated"),
        isoannotlite_script = config["isoannotlite_script"]
    log:
        os.path.join(WORKDIR, "logs/" + RECOLLAPSED_DIR + "/sqanti3/{tissue}_IsoAnnotLite.log")
    threads: get_resource("run_IsoAnnotLite", "threads", DEFAULT_THREADS)
    resources:
        mem=get_resource("run_IsoAnnotLite", "mem", DEFAULT_MEMORY),
        partition=get_resource("run_IsoAnnotLite", "partition", DEFAULT_PARTITION),
        time=get_resource("run_IsoAnnotLite", "time", DEFAULT_TIME),
        slurm_extra=get_resource("run_IsoAnnotLite", "slurm_extra", DEFAULT_SLURM_EXTRA)
    shell:
        """
        python {params.isoannotlite_script} \
            {input.filtered_gtf} \
            {input.filtered_classification} \
            {input.junctions} \
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
