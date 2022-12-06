from wsgiref.util import request_uri
import tqdm
import pandas as pd
import pysam
from Bio import SeqIO, Seq, SeqRecord
import logging
from multiprocessing import Pool
from multiprocessing.managers import BaseManager, DictProxy
from collections import defaultdict
from functools import partial

from dmg_reads.utils import is_debug, calc_chunksize, initializer

import cProfile as profile
import pstats

log = logging.getLogger("my_logger")


class MyManager(BaseManager):
    pass


def ddict():
    return defaultdict(ddict)


MyManager.register("ddict", ddict, DictProxy)


def get_alns(params, bam, reads, refs_tax, refs_damaged, threads=1):
    bam, references = params
    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)
    for reference in references:
        for aln in samfile.fetch(
            contig=reference, multiple_iterators=False, until_eof=True
        ):
            # create read
            # Check if reference is damaged
            aln_reference_name = reference
            aln_qname = aln.qname
            is_damaged = "non-damaged"
            if aln_reference_name in refs_damaged:
                is_damaged = "damaged"

            if reads[refs_tax[aln_reference_name]][aln_qname]:
                dmg = reads[refs_tax[aln_reference_name]][aln_qname]["is_damaged"]
                if dmg == is_damaged:
                    continue
                else:
                    reads[refs_tax[aln_reference_name]][aln_qname][
                        "is_damaged"
                    ] = "multi"
            else:
                seq = Seq.Seq(aln.seq)
                qual = aln.query_qualities
                if aln.is_reverse:
                    seq = seq.reverse_complement()
                    qual = qual[::-1]
                reads[refs_tax[aln_reference_name]][aln_qname] = {
                    "seq": seq,
                    "qual": qual,
                    "is_damaged": is_damaged,
                }
    samfile.close()


def get_read_by_taxa(
    bam, refs_tax, refs, refs_damaged, ref_bam_dict, chunksize=None, threads=1
):
    prof = profile.Profile()
    prof.enable()

    if (chunksize is not None) and ((len(refs) // chunksize) > threads):
        c_size = chunksize
    else:
        c_size = calc_chunksize(n_workers=threads, len_iterable=len(refs), factor=4)

    ref_chunks = [refs[i : i + c_size] for i in range(0, len(refs), c_size)]
    params = zip([bam] * len(ref_chunks), ref_chunks)

    mgr = MyManager()
    mgr.start()
    reads = mgr.ddict()

    if is_debug():
        data = list(
            map(
                partial(
                    get_alns,
                    bam=bam,
                    reads=reads,
                    refs_tax=refs_tax,
                    refs_damaged=refs_damaged,
                    threads=threads,
                ),
                params,
            )
        )
    else:
        p = Pool(
            threads,
            initializer=initializer,
            initargs=([params, refs_tax, refs_damaged],),
        )

        data = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(
                        get_alns,
                        bam=bam,
                        reads=reads,
                        refs_tax=refs_tax,
                        refs_damaged=refs_damaged,
                        threads=threads,
                    ),
                    params,
                    chunksize=1,
                ),
                total=len(ref_chunks),
                leave=False,
                ncols=80,
                desc="References processed",
            )
        )

    pool.close()
    pool.join()

    print(reads)
    prof.disable()
    # print profiling output
    stats = pstats.Stats(prof).sort_stats("tottime")
    stats.print_stats(10)  # top 10 rows
    return reads
