#!/usr/bin/env python

import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
import parse_protocol
import pysam
import utils
from __init__ import ASSAY


def get_tag_barcode_mismatch_dict(barcode_dict):
    mismatch_dict = {}
    for seq_id, seq in barcode_dict.items():
        for mismatch_seq in parse_protocol.findall_mismatch(seq, n_mismatch=2):
            mismatch_dict[mismatch_seq] = seq_id
    return mismatch_dict


def read_fasta(fasta_file, equal=False):
    """
    Args:
        equal: if True, seq in fasta must have equal length
    Returns:
        {seq_id: seq} dict
    """
    fa_dict = {}
    length = None
    with pysam.FastxFile(fasta_file) as infile:
        for index, record in enumerate(infile):
            seq = record.sequence
            if index == 0:
                length = len(seq)
            if equal:
                if length != len(seq):
                    raise Exception(f"{fasta_file} have different seq length")
            fa_dict[record.name] = seq
    return fa_dict, length


class TagBarcode:
    def __init__(self, args):
        self.args = args

    def run(self):
        barcode_dict, barcode_length = read_fasta(self.args.tag_barcode_fasta, equal=True)
        pattern_dict = parse_protocol.parse_pattern(self.args.r2_pattern)
        feature_bc_slice = pattern_dict["C"][0]
        bc_pattern_len = feature_bc_slice.stop - feature_bc_slice.start
        if bc_pattern_len != barcode_length:
            raise ValueError(f"""The length of tag barcode in r2_pattern({bc_pattern_len}) != 
                length of feature barcode in tag_barcode_fasta({barcode_length})""")
        mismatch_dict = get_tag_barcode_mismatch_dict(barcode_dict)
        bcs = utils.read_one_col(self.args.match_barcode)
        bcs = set(bcs)

        total_reads = 0
        tag_reads = 0
        incell_reads = 0
        bc_antibody_umi = utils.nested_defaultdict(dim=3)
        with pysam.FastxFile(self.args.fq) as infile:
            for record in infile:
                total_reads += 1
                seq = record.sequence
                bc, umi = record.name.split(":")[0:2]
                tag_barcode = seq[feature_bc_slice]
                if tag_barcode in mismatch_dict:
                    tag_reads += 1
                    if bc in bcs:
                        antibody = mismatch_dict[tag_barcode]
                        bc_antibody_umi[bc][antibody][umi] += 1
                        incell_reads += 1

        # median umi per cell
        bc_umi = defaultdict(int)
        bc_antibody = utils.nested_defaultdict(dim=2)
        for bc in bc_antibody_umi:
            for antibody in bc_antibody_umi[bc]:
                umi_count = len(bc_antibody_umi[bc][antibody])
                bc_umi[bc] += umi_count
                bc_antibody[bc][antibody] = umi_count
        median_umi_per_cell = np.median(list(bc_umi.values()))

        # metrics
        stats = {}
        stats["Number of Cells"] = len(bcs)
        stats["Fraction Tag Reads"] = utils.get_frac(tag_reads / total_reads)
        stats["Fraction Tag Reads in Cell"] = utils.get_frac(incell_reads / tag_reads)
        stats["Median UMI per Cell"] = median_umi_per_cell
        utils.write_multiqc(stats, self.args.sample, ASSAY, "tag_barcode" + ".stats")

        # write csv
        df = pd.DataFrame.from_dict(bc_antibody)
        df.fillna(0, inplace=True)
        df = df.astype(int)
        df.to_csv(f"{self.args.sample}.tag_barcode.csv.gz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--r2_pattern", required=True)
    parser.add_argument("--tag_barcode_fasta", required=True)
    parser.add_argument("--fq", help="R2 read fastq.", required=True)
    parser.add_argument("--match_barcode", help="match barcodes.", required=True)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    TagBarcode(args).run()
