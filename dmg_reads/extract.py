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

from dmg_reads.utils import is_debug, initializer

# import cProfile as profile
# import pstats

log = logging.getLogger("my_logger")


class MyManager(BaseManager):
    pass


def ddict():
    return defaultdict(ddict)


MyManager.register("ddict", ddict, DictProxy)


def get_alns(reference, bam, reads, refs_tax, refs_damaged):

    samfile = pysam.AlignmentFile(bam, "rb")
    for aln in samfile.fetch(
        reference=reference, multiple_iterators=False, until_eof=False
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
                reads[refs_tax[aln_reference_name]][aln_qname]["is_damaged"] = "multi"
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


def get_read_by_taxa(bam, refs_tax, refs, refs_damaged, ref_bam_dict, threads=1):
    # prof = profile.Profile()
    # prof.enable()
    # reads = defaultdict(lambda: defaultdict(dict))

    pool = Pool(processes=threads)
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
                ),
                refs,
            )
        )
    else:

        p = Pool(
            threads,
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
                    ),
                    refs,
                    chunksize=1,
                ),
                total=len(refs),
                leave=False,
                ncols=80,
                desc="References processed",
            )
        )

    pool.close()
    pool.join()

    # prof.disable()
    # # print profiling output
    # stats = pstats.Stats(prof).sort_stats("tottime")
    # stats.print_stats(10)  # top 10 rows
    return reads
