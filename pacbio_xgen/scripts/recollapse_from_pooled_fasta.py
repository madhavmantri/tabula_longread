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


def remap_value(value, mapping):
    if value in mapping:
        return mapping[value]
    if "," in value:
        return ",".join(mapping.get(part, part) for part in value.split(","))
    return value


def transcript_id_from_attributes(attributes):
    for key in ("ID", "transcript_id"):
        value = attributes.get(key)
        if value:
            return value
    parent = attributes.get("Parent")
    if parent:
        return parent.split(",")[0]
    raise ValueError("Could not determine transcript ID from GFF attributes")


def prefix_header(header, tissue):
    prefix = tissue + "|"
    if "|" in header:
        first, rest = header.split("|", 1)
        if rest.startswith(first + ":"):
            suffix = rest[len(first):]
            return "{0}{1}|{0}{1}{2}".format(prefix, first, suffix)
        return "{0}{1}|{2}".format(prefix, first, rest)
    return "{0}{1}".format(prefix, header)


def prefix_fasta_command(args):
    Path(args.fasta_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.fasta) as src, open(args.fasta_out, "w") as dst:
        for line in src:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                dst.write(">{0}\n".format(prefix_header(header, args.tissue)))
            else:
                dst.write(line)


def load_pooled_mapping(pooled_group_path, tissue):
    prefix = tissue + "|"
    mapping = {}
    with open(pooled_group_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            pooled_id, members = line.rstrip("\n").split("\t", 1)
            for member in members.split(","):
                if not member.startswith(prefix):
                    continue
                original_id = member[len(prefix):].split("|")[0]
                mapping[original_id] = pooled_id
    if not mapping:
        raise ValueError("No pooled group assignments found for tissue '{0}'".format(tissue))
    return mapping


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


def rewrite_command(args):
    mapping = load_pooled_mapping(args.pooled_group, args.tissue)

    Path(args.mapping_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.mapping_out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tissue", "old_id", "new_id"])
        for old_id, new_id in sorted(mapping.items()):
            writer.writerow([args.tissue, old_id, new_id])

    rewrite_gff(mapping, args.gff, args.gff_out)
    rewrite_group(mapping, args.group, args.group_out)
    rewrite_abundance(mapping, args.abundance, args.abundance_out)


def main():
    parser = argparse.ArgumentParser(description="Recollapse pooled collapse.fasta inputs and project pooled IDs back to tissues.")
    subparsers = parser.add_subparsers(dest="command")

    prefix_parser = subparsers.add_parser("prefix-fasta")
    prefix_parser.add_argument("--tissue", required=True)
    prefix_parser.add_argument("--fasta", required=True)
    prefix_parser.add_argument("--fasta-out", required=True)
    prefix_parser.set_defaults(func=prefix_fasta_command)

    rewrite_parser = subparsers.add_parser("rewrite")
    rewrite_parser.add_argument("--tissue", required=True)
    rewrite_parser.add_argument("--pooled-group", required=True)
    rewrite_parser.add_argument("--gff", required=True)
    rewrite_parser.add_argument("--group", required=True)
    rewrite_parser.add_argument("--abundance", required=True)
    rewrite_parser.add_argument("--mapping-out", required=True)
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
