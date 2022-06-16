"""
 Copyright (c) 2022 Antonio Fernandez-Guerra

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 """


import logging
from dmg_reads.utils import get_arguments, create_output_files
from dmg_reads.lib import load_mdmg_results, filter_damaged_taxa
import pandas as pd
import os
from pathlib import Path
import pysam
import numpy as np
from Bio import SeqIO, Seq, SeqRecord
import gzip
from mimetypes import guess_type
from functools import partial
import tqdm

log = logging.getLogger("my_logger")


def main():

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = get_arguments()
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    log.info("Loading metaDMG results...")
    mdmg_results = load_mdmg_results(args.metaDMG_results)
    # find which taxon are damaged

    damaged_taxa = filter_damaged_taxa(
        df=mdmg_results,
        filter_conditions=args.metaDMG_filter,
    )

    refs = damaged_taxa["reference"].to_list()
    logging.info(f"Loading BAM file")
    save = pysam.set_verbosity(0)
    bam = args.bam
    samfile = pysam.AlignmentFile(bam, "rb")

    references = samfile.references

    chr_lengths = []
    for chrom in samfile.references:
        chr_lengths.append(samfile.get_reference_length(chrom))
    max_chr_length = np.max(chr_lengths)
    pysam.set_verbosity(save)
    ref_lengths = None
    if not samfile.has_index():
        logging.info(f"BAM index not found. Indexing...")
        if max_chr_length > 536870912:
            logging.info(f"A reference is longer than 2^29, indexing with csi")
            pysam.index(bam, "-c")
        else:
            pysam.index(bam)
    reads_damaged_seen = {}
    reads_nondamaged_seen = {}

    out_files = create_output_files(prefix=args.prefix, bam=args.bam)
    encoding = guess_type(out_files["fastq_damaged"])[1]
    _open = partial(gzip.open, mode="wt") if encoding == "gzip" else open

    with _open(out_files["fastq_damaged"]) as dmg_fh, _open(
        out_files["fastq_nondamaged"]
    ) as non_dmg_fh:
        for aln in tqdm.tqdm(
            samfile.fetch(until_eof=True),
            total=samfile.mapped,
            leave=False,
            ncols=80,
            desc=f"Alignments processed",
        ):
            if aln.reference_name in refs:
                if aln.qname not in reads_damaged_seen:
                    reads_damaged_seen[aln.qname] = aln.qname
                    seq = Seq.Seq(aln.seq)
                    qual = aln.query_qualities
                    if aln.is_reverse:
                        seq = seq.reverse_complement()
                        qual = qual[::-1]
                    rec = SeqRecord.SeqRecord(seq, aln.qname, "", "")
                    rec.letter_annotations["phred_quality"] = qual
                    SeqIO.write(rec, dmg_fh, "fastq")
            else:
                if aln.qname not in reads_damaged_seen:
                    reads_nondamaged_seen[aln.qname] = aln.qname
                    seq = Seq.Seq(aln.seq)
                    qual = aln.query_qualities
                    if aln.is_reverse:
                        seq = seq.reverse_complement()
                        qual = qual[::-1]
                    rec = SeqRecord.SeqRecord(seq, aln.qname, "", "")
                    rec.letter_annotations["phred_quality"] = qual
                    SeqIO.write(rec, non_dmg_fh, "fastq")

    logging.info(
        f"Extracted {len(reads_damaged_seen):,} damaged reads to {out_files['fastq_damaged']}"
    )
    logging.info(
        f"Extracted {len(reads_nondamaged_seen):,} non-damaged reads to {out_files['fastq_nondamaged']}"
    )
    logging.info(f"Done!")


if __name__ == "__main__":
    main()
