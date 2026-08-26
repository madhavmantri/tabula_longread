#!/usr/bin/env python3

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


def parse_attributes(raw):
    attributes = OrderedDict()
    for field in raw.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
            value = value.strip('"')
        else:
            key, value = field, ""
        attributes[key] = value
    return attributes


def format_attributes(attributes):
    return " ".join(f'{key} "{value}";' if value != "" else key for key, value in attributes.items())


def remap_token(value, mapping):
    if value in mapping:
        return mapping[value]
    if "," in value:
        return ",".join(mapping.get(part, part) for part in value.split(","))
    return value


def rewrite_gff(mapping, input_path, output_path):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(input_path) as src, open(output_path, "w") as dst:
        for line in src:
            if line.startswith("#") or not line.strip():
                dst.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                dst.write(line)
                continue
            attributes = parse_attributes(fields[8])
            for key, value in list(attributes.items()):
                attributes[key] = remap_token(value, mapping)
            fields[8] = format_attributes(attributes)
            dst.write("\t".join(fields) + "\n")


def rewrite_group(mapping, input_path, output_path):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(input_path) as src, open(output_path, "w") as dst:
        for line in src:
            if not line.strip():
                dst.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if fields:
                fields[0] = mapping.get(fields[0], fields[0])
            dst.write("\t".join(fields) + "\n")


def rewrite_abundance(mapping, input_path, output_path):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(input_path) as src, open(output_path, "w") as dst:
        header = src.readline()
        if not header:
            return
        dst.write(header)
        columns = header.rstrip("\n").split("\t")
        pbid_index = columns.index("pbid") if "pbid" in columns else 0
        for line in src:
            if not line.strip():
                dst.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > pbid_index:
                fields[pbid_index] = mapping.get(fields[pbid_index], fields[pbid_index])
            dst.write("\t".join(fields) + "\n")


def transcript_id_from_attributes(attributes):
    for key in ("ID", "transcript_id"):
        value = attributes.get(key)
        if value:
            return value
    parent = attributes.get("Parent")
    if parent:
        return parent.split(",")[0]
    raise ValueError("Could not determine transcript ID from GFF attributes")


def prefix_command(args):
    prefix = args.tissue + "|"
    mapping = {}
    with open(args.gff) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            tx_id = transcript_id_from_attributes(parse_attributes(fields[8]))
            mapping.setdefault(tx_id, prefix + tx_id)

    rewrite_gff(mapping, args.gff, args.gff_out)
    rewrite_group(mapping, args.group, args.group_out)
    rewrite_abundance(mapping, args.abundance, args.abundance_out)


def load_transcripts(gff_paths):
    transcripts = OrderedDict()
    exons = {}
    template_lines = {}

    for gff_path in gff_paths:
        with open(gff_path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                seqid, source, feature, start, end, score, strand, phase, raw_attrs = fields
                attributes = parse_attributes(raw_attrs)
                tx_id = transcript_id_from_attributes(attributes)
                start = int(start)
                end = int(end)
                if feature == "transcript":
                    transcripts[tx_id] = {
                        "seqid": seqid,
                        "source": source,
                        "feature": feature,
                        "start": start,
                        "end": end,
                        "score": score,
                        "strand": strand,
                        "phase": phase,
                        "attributes": attributes,
                    }
                    template_lines[tx_id] = [line.rstrip("\n")]
                elif feature == "exon":
                    exons.setdefault(tx_id, []).append((start, end))
                    template_lines.setdefault(tx_id, []).append(line.rstrip("\n"))

    records = []
    for tx_id, tx in transcripts.items():
        exon_coords = tuple(sorted(exons.get(tx_id, [(tx["start"], tx["end"])])))
        signature = (tx["seqid"], tx["strand"], exon_coords)
        records.append((tx_id, signature, tx, exon_coords, template_lines.get(tx_id, [])))
    return records


def build_catalog_command(args):
    signature_to_global = OrderedDict()
    rows = []
    catalog_lines = []

    for tx_id, signature, tx, exon_coords, template_lines in load_transcripts(args.gff):
        if signature not in signature_to_global:
            global_id = "PB.pool.{0}".format(len(signature_to_global) + 1)
            signature_to_global[signature] = global_id

            tx_attrs = OrderedDict(tx["attributes"])
            tx_attrs["ID"] = global_id
            tx_fields = [
                tx["seqid"], tx["source"], tx["feature"], str(tx["start"]), str(tx["end"]),
                tx["score"], tx["strand"], tx["phase"], format_attributes(tx_attrs)
            ]
            catalog_lines.append("\t".join(tx_fields))
            for exon_start, exon_end in exon_coords:
                exon_attrs = OrderedDict([("Parent", global_id)])
                exon_fields = [
                    tx["seqid"], tx["source"], "exon", str(exon_start), str(exon_end),
                    ".", tx["strand"], ".", format_attributes(exon_attrs)
                ]
                catalog_lines.append("\t".join(exon_fields))

        rows.append((tx_id, signature_to_global[signature], tx["seqid"], tx["strand"], len(exon_coords)))

    Path(args.map_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.map_out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["old_id", "new_id", "seqid", "strand", "exon_count"])
        writer.writerows(rows)

    Path(args.catalog_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.catalog_out, "w") as handle:
        for line in catalog_lines:
            handle.write(line + "\n")


def load_mapping(mapping_path):
    mapping = {}
    with open(mapping_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            mapping[row["old_id"]] = row["new_id"]
    return mapping


def rewrite_command(args):
    mapping = load_mapping(args.mapping)
    rewrite_gff(mapping, args.gff, args.gff_out)
    rewrite_group(mapping, args.group, args.group_out)
    rewrite_abundance(mapping, args.abundance, args.abundance_out)


def main():
    parser = argparse.ArgumentParser(description="Prefix tissue-specific PacBio IDs and build a pooled global catalog.")
    subparsers = parser.add_subparsers(dest="command")

    prefix_parser = subparsers.add_parser("prefix")
    prefix_parser.add_argument("--tissue", required=True)
    prefix_parser.add_argument("--gff", required=True)
    prefix_parser.add_argument("--group", required=True)
    prefix_parser.add_argument("--abundance", required=True)
    prefix_parser.add_argument("--gff-out", required=True)
    prefix_parser.add_argument("--group-out", required=True)
    prefix_parser.add_argument("--abundance-out", required=True)
    prefix_parser.set_defaults(func=prefix_command)

    build_catalog_parser = subparsers.add_parser("build-catalog")
    build_catalog_parser.add_argument("--gff", action="append", required=True)
    build_catalog_parser.add_argument("--map-out", required=True)
    build_catalog_parser.add_argument("--catalog-out", required=True)
    build_catalog_parser.set_defaults(func=build_catalog_command)

    rewrite_parser = subparsers.add_parser("rewrite")
    rewrite_parser.add_argument("--mapping", required=True)
    rewrite_parser.add_argument("--gff", required=True)
    rewrite_parser.add_argument("--group", required=True)
    rewrite_parser.add_argument("--abundance", required=True)
    rewrite_parser.add_argument("--gff-out", required=True)
    rewrite_parser.add_argument("--group-out", required=True)
    rewrite_parser.add_argument("--abundance-out", required=True)
    rewrite_parser.set_defaults(func=rewrite_command)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        parser.exit(2)
    args.func(args)


if __name__ == "__main__":
    main()
