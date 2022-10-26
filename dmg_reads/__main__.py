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

from dmg_reads.utils import get_arguments, create_output_files, splitkeep, fast_flatten
from dmg_reads.lib import load_mdmg_results, filter_damaged_taxa
from dmg_reads.defaults import valid_ranks
from dmg_reads.extract import get_read_by_taxa
import pandas as pd
import os
import pysam
import numpy as np
from Bio import SeqIO, Seq, SeqRecord
import gzip
from mimetypes import guess_type
from functools import partial
import tqdm
from collections import defaultdict
import re

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

    # Check that rank and taxonomy file are both present
    if args.rank and not args.taxonomy_file:
        log.error("--rank requires --taxonomy")
        exit(1)
    if args.taxonomy_file and not args.rank:
        log("Error: --taxonomy-file requires --rank")
        exit(1)

    log.info("Loading metaDMG results...")
    mdmg_results = load_mdmg_results(args.metaDMG_results)
    # find which taxon are damaged

    damaged_taxa = filter_damaged_taxa(
        df=mdmg_results,
        filter_conditions=args.metaDMG_filter,
    )

    if args.taxonomy_file:
        log.info("Loading taxonomy data...")
        refs_tax = defaultdict()
        taxonomy = pd.read_csv(
            args.taxonomy_file,
            sep="\t",
            header=None,
            names=["reference", "taxonomy"],
        )

        # get get_ranks
        ranks = {valid_ranks[x]: args.rank[x] for x in args.rank if x in valid_ranks}
        # get refs that belong to this taxonomy
        # loop over all rows
        for i, row in tqdm.tqdm(
            taxonomy.iterrows(),
            total=taxonomy.shape[0],
            desc="Taxa processed",
            leave=False,
        ):
            taxs = row[1].split(";")
            for r in taxs:
                r1 = splitkeep(r, "__")
                if r1[0] in ranks and r1[1] in ranks[r1[0]]:
                    v = re.sub("[^0-9a-zA-Z]+", "_", r1[1])
                    refs_tax[row[0]] = f"{r1[0]}{v}"
                    # refs[row[0]] = f
    else:
        ranks = None

    logging.info("Loading BAM file...")
    save = pysam.set_verbosity(0)
    bam = args.bam
    samfile = pysam.AlignmentFile(bam, "rb")

    chr_lengths = []
    for chrom in samfile.references:
        chr_lengths.append(samfile.get_reference_length(chrom))
    max_chr_length = np.max(chr_lengths)

    # Check if BAM files is not sorted by coordinates, sort it by coordinates
    if not samfile.header["HD"]["SO"] == "coordinate":
        log.info("BAM file is not sorted by coordinates, sorting it...")
        sorted_bam = bam.replace(".bam", ".dr-sorted.bam")
        pysam.sort(
            "-@", str(args.threads), "-m", str(args.sort_memory), "-o", sorted_bam, bam
        )
        bam = sorted_bam
        samfile = pysam.AlignmentFile(bam, "rb")

    if not samfile.has_index():
        logging.info("BAM index not found. Indexing...")
        if max_chr_length > 536870912:
            logging.info("A reference is longer than 2^29, indexing with csi")
            pysam.index(bam, "-c", "-@", str(args.threads))
        else:
            pysam.index(
                bam, "-@", str(args.threads)
            )  # Need to reload the samfile after creating index
            log.info("Re-loading BAM file")
            samfile = pysam.AlignmentFile(bam, "rb")
    pysam.set_verbosity(save)

    refs_bam = [
        chrom.contig for chrom in samfile.get_index_statistics() if chrom.mapped > 0
    ]

    if args.taxonomy_file:
        refs_damaged = set(refs_tax.keys()).intersection(
            set(damaged_taxa["reference"].to_list())
        )
        refs_non_damaged = (
            set(refs_tax.keys()).intersection(set(refs_bam)) - refs_damaged
        )
        if not refs_damaged and not refs_non_damaged:
            log.info("No references found in BAM file")
            exit(0)
    else:
        refs_damaged = damaged_taxa["reference"].to_list()
        refs_non_damaged = set(refs_bam) - set(refs_damaged)
        refs_tax = {ref: "all" for ref in refs_bam}

    if args.only_damaged:
        refs = refs_damaged
    else:
        refs = fast_flatten([refs_non_damaged, refs_damaged])

    out_files = create_output_files(prefix=args.prefix, bam=args.bam, taxon=ranks)

    for file in out_files:
        # file exists deleted
        if os.path.exists(out_files[file]):
            os.remove(out_files[file])

    log.info("Processing reads...")

    reads = get_read_by_taxa(
        samfile=samfile, refs=refs, refs_tax=refs_tax, refs_damaged=refs_damaged
    )
    # write reads

    count = defaultdict(int)
    logging.info("Saving reads...")
    if args.taxonomy_file:
        desc = "Taxa processed"
    else:
        desc = "References processed"

    for tax in tqdm.tqdm(reads, ncols=80, desc=desc, leave=False, total=len(reads)):
        if args.taxonomy_file:
            fastq_damaged = out_files[f"fastq_damaged_{tax}"]
            fastq_nondamaged = out_files[f"fastq_nondamaged_{tax}"]
            fastq_multi = out_files[f"fastq_multi_{tax}"]
            fastq_combined = out_files[f"fastq_combined_{tax}"]
        else:
            fastq_damaged = out_files["fastq_damaged"]
            fastq_nondamaged = out_files["fastq_nondamaged"]
            fastq_multi = out_files["fastq_multi"]
            fastq_combined = out_files["fastq_combined"]

        encoding = guess_type(fastq_damaged)[1]
        _open = partial(gzip.open, mode="at") if encoding == "gzip" else open
        # TODO: clean this up
        with _open(fastq_damaged) as f_damaged, _open(
            fastq_nondamaged
        ) as f_nondamaged, _open(fastq_multi) as f_multi, _open(
            fastq_combined
        ) as f_combined:
            for read in tqdm.tqdm(
                reads[tax],
                ncols=80,
                desc="Reads written",
                leave=False,
                total=len(reads[tax]),
                ascii="░▒█",
            ):
                rec = SeqRecord.SeqRecord(reads[tax][read]["seq"], read, "", "")
                rec.letter_annotations["phred_quality"] = reads[tax][read]["qual"]
                if args.combine:
                    SeqIO.write(rec, f_combined, "fastq")
                    if args.taxonomy_file:
                        count[out_files[f"fastq_combined_{tax}"]] += 1
                    else:
                        count[out_files["fastq_combined"]] += 1
                else:
                    if reads[tax][read]["is_damaged"] == "damaged":
                        SeqIO.write(rec, f_damaged, "fastq")
                        if args.taxonomy_file:
                            count[out_files[f"fastq_damaged_{tax}"]] += 1
                        else:
                            count[out_files["fastq_damaged"]] += 1
                    elif reads[tax][read]["is_damaged"] == "non-damaged":
                        SeqIO.write(rec, f_nondamaged, "fastq")
                        if args.taxonomy_file:
                            count[out_files[f"fastq_nondamaged_{tax}"]] += 1
                        else:
                            count[out_files["fastq_nondamaged"]] += 1
                    elif reads[tax][read]["is_damaged"] == "multi":
                        SeqIO.write(rec, f_multi, "fastq")
                        if args.taxonomy_file:
                            count[out_files[f"fastq_multi_{tax}"]] += 1
                        else:
                            count[out_files["fastq_multi"]] += 1

    for file in out_files:
        if count[out_files[file]] and count[out_files[file]] > 0:
            pass
        else:
            os.remove(out_files[file])

    logging.info("Done!")


if __name__ == "__main__":
    main()
