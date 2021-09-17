#!/usr/bin/env python3

#
# Copyright (C) 2021 Jerry Hoogenboom
#
# This file is part of FDSTools, data analysis tools for Massively
# Parallel Sequencing of forensic DNA markers.
#
# FDSTools is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FDSTools is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FDSTools.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Name microhaplotypes.

Adds a new column 'microhaplotype' to the output, which contains the
base calls on given microhaplotype positions. The POSITIONS file should
contain two tab-separated values on each line: a marker name and a
genomic position number (without chromosome number). If there is no
[genome_position] set in the library file, the first nucleotide in the
target fragment is position 1.
"""
import argparse
import re
import sys

from errno import EPIPE
from strnaming import libsequence

from ..lib.cli import add_sequence_format_args, add_input_output_args, get_input_output_files
from ..lib.io import get_column_ids
from ..lib.seq import SEQ_SPECIAL_VALUES, ensure_sequence_format

__version__ = "1.0.0"


# Pattern to match a variant in a microhaplotype-defining position.
PAT_VARIANT = re.compile("^N\d+([ACGT])$")


def get_n_refseqs(infile, library):
    refseqs = {}
    for line in infile:
        marker, position = line.rstrip("\r\n").split("\t")
        markerlib = library.get_range(marker)
        position = int(position)
        offset = libsequence.get_genome_pos(markerlib.location, position, invert=True)
        if marker not in refseqs:
            refseqs[marker] = list(markerlib.refseq)
        refseqs[marker][offset] = "N"
    return {marker: "".join(refseq) for marker, refseq in refseqs.items()}
#get_n_refseqs


def name_microhaplotypes(infile, outfile, refseqs, library):
    column_names = infile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return  # Empty file.
    colid_marker, colid_sequence = get_column_ids(column_names, "marker", "sequence")
    colid_microhaplotype = get_column_ids(column_names, "microhaplotype", optional=True)
    if colid_microhaplotype is None:
        column_names.append("microhaplotype")
        colid_microhaplotype = -1
    outfile.write("\t".join(column_names) + "\n")
    for line in infile:
        cols = line.rstrip("\r\n").split("\t")
        if colid_microhaplotype == -1:
            cols.append("")
        marker = cols[colid_marker]
        sequence = ensure_sequence_format(
            cols[colid_sequence], "raw", library=library, marker=marker)
        if marker in refseqs and sequence not in SEQ_SPECIAL_VALUES:
            cols[colid_microhaplotype] = "".join(variant.group(1) for variant in
                map(PAT_VARIANT.match,
                    libsequence.call_variants(refseqs[marker], sequence, location=("M", 1),
                    match_score=1, mismatch_score=-3, gap_open_score=-7, gap_extend_score=-2))
                if variant is not None)
        outfile.write("\t".join(cols) + "\n")
#name_microhaplotypes


def add_arguments(parser):
    add_sequence_format_args(parser, default_format="raw", force=True, require_library=True)
    parser.add_argument("positions", metavar="POSITIONS", type=argparse.FileType("rt", encoding="UTF-8"),
        help="file containing a list of expected variant positions")
    add_input_output_args(parser, single_in=True, batch_support=True, report_out=False)
#add_arguments


def run(args):
    gen = get_input_output_files(args, single_in=True, batch_support=True)
    if not gen:
        raise ValueError("please specify an input file, or pipe in the output of another program")

    # Read list of expected variant positions.
    refseqs = get_n_refseqs(args.positions, args.library)

    for tag, infiles, outfile in gen:
        # TODO: Aggregate data from all infiles of each sample.
        if len(infiles) > 1:
            raise ValueError("multiple input files for sample '%s' specified " % tag)
        try:
            infile = sys.stdin if infiles[0] == "-" else open(infiles[0], "rt", encoding="UTF-8")
            name_microhaplotypes(infile, outfile, refseqs, args.library)
            if infile != sys.stdin:
                infile.close()
        except IOError as e:
            if e.errno == EPIPE:
                continue
            raise
#run
