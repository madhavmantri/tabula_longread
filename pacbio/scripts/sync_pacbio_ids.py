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


def transcript_id_from_attributes(attributes):
    for key in ("ID", "transcript_id"):
        value = attributes.get(key)
        if value:
            return value
    parent = attributes.get("Parent")
    if parent:
        return parent.split(",")[0]
    raise ValueError("Could not determine transcript ID from GFF attributes")


def load_transcript_signatures(gff_path):
    transcripts = {}
    exons = {}
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
                transcripts[tx_id] = (seqid, strand, start, end)
            elif feature == "exon":
                exons.setdefault(tx_id, []).append((start, end))

    signatures = {}
    for tx_id, (seqid, strand, tx_start, tx_end) in transcripts.items():
        exon_coords = tuple(sorted(exons.get(tx_id, [(tx_start, tx_end)])))
        tx_bounds = (exon_coords[0][0], exon_coords[-1][1]) if exon_coords else (tx_start, tx_end)
        signatures[tx_id] = (seqid, strand, tx_bounds, exon_coords)
    return signatures


def build_map(args):
    if len(args.tissue) != len(args.gff):
        raise ValueError("Each --tissue requires a matching --gff")

    rows = []
    signature_to_new_id = OrderedDict()

    for tissue, gff_path in sorted(zip(args.tissue, args.gff), key=lambda item: item[0]):
        signatures = load_transcript_signatures(gff_path)
        for old_id, signature in sorted(signatures.items()):
            if signature not in signature_to_new_id:
                signature_to_new_id[signature] = f"PB.sync.{len(signature_to_new_id) + 1}"
            new_id = signature_to_new_id[signature]
            seqid, strand, (start, end), exon_coords = signature
            rows.append(
                {
                    "tissue": tissue,
                    "old_id": old_id,
                    "new_id": new_id,
                    "seqid": seqid,
                    "strand": strand,
                    "start": start,
                    "end": end,
                    "exon_count": len(exon_coords),
                    "exons": ",".join(f"{exon_start}-{exon_end}" for exon_start, exon_end in exon_coords),
                }
            )

    Path(args.map_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.map_out, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["tissue", "old_id", "new_id", "seqid", "strand", "start", "end", "exon_count", "exons"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_mapping(mapping_path, tissue):
    mapping = {}
    with open(mapping_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["tissue"] == tissue:
                mapping[row["old_id"]] = row["new_id"]
    if not mapping:
        raise ValueError(f"No mapping entries found for tissue '{tissue}'")
    return mapping


def remap_value(value, mapping):
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
                attributes[key] = remap_value(value, mapping)
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


def rewrite(args):
    mapping = load_mapping(args.mapping, args.tissue)
    rewrite_gff(mapping, args.gff, args.gff_out)
    rewrite_group(mapping, args.group, args.group_out)
    rewrite_abundance(mapping, args.abundance, args.abundance_out)


def main():
    parser = argparse.ArgumentParser(description="Synchronize PacBio transcript IDs across tissues.")
    subparsers = parser.add_subparsers(dest="command")

    build_map_parser = subparsers.add_parser("build-map")
    build_map_parser.add_argument("--tissue", action="append", required=True)
    build_map_parser.add_argument("--gff", action="append", required=True)
    build_map_parser.add_argument("--map-out", required=True)
    build_map_parser.set_defaults(func=build_map)

    rewrite_parser = subparsers.add_parser("rewrite")
    rewrite_parser.add_argument("--mapping", required=True)
    rewrite_parser.add_argument("--tissue", required=True)
    rewrite_parser.add_argument("--gff", required=True)
    rewrite_parser.add_argument("--group", required=True)
    rewrite_parser.add_argument("--abundance", required=True)
    rewrite_parser.add_argument("--gff-out", required=True)
    rewrite_parser.add_argument("--group-out", required=True)
    rewrite_parser.add_argument("--abundance-out", required=True)
    rewrite_parser.set_defaults(func=rewrite)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        parser.exit(2)
    args.func(args)


if __name__ == "__main__":
    main()
