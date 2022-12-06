import tqdm
import pysam
from Bio import Seq
import logging
from multiprocessing import Pool
from collections import defaultdict
from functools import partial
from dmg_reads.utils import is_debug, calc_chunksize, initializer

# import cProfile as profile
# import pstats

log = logging.getLogger("my_logger")


def get_alns(params, refs_tax, refs_damaged, threads=1):

    reads = defaultdict(lambda: defaultdict(dict))
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
    return dict(reads)


def merge_dicts(dicts):

    reads = defaultdict(lambda: defaultdict(dict))

    for d in dicts:
        for tax, tax_reads in d.items():
            for read, read_info in tax_reads.items():
                if reads[tax][read]:
                    dmg = reads[tax][read]["is_damaged"]
                    if dmg == read_info["is_damaged"]:
                        continue
                    else:
                        reads[tax][read]["is_damaged"] = "multi"
                else:
                    reads[tax][read] = read_info
    return dict(reads)


def get_read_by_taxa(
    bam, refs_tax, refs, refs_damaged, ref_bam_dict, chunksize=None, threads=1
):
    # prof = profile.Profile()
    # prof.enable()

    if (chunksize is not None) and ((len(refs) // chunksize) > threads):
        c_size = chunksize
    else:
        c_size = calc_chunksize(n_workers=threads, len_iterable=len(refs), factor=4)

    ref_chunks = [refs[i : i + c_size] for i in range(0, len(refs), c_size)]

    params = zip([bam] * len(ref_chunks), ref_chunks)

    if is_debug():
        data = list(
            map(
                partial(
                    get_alns,
                    refs_tax=refs_tax,
                    refs_damaged=refs_damaged,
                    threads=threads,
                ),
                params,
            )
        )
    else:
        logging.info(
            f"Processing {len(ref_chunks):,} chunks of {c_size:,} references each..."
        )
        p = Pool(
            threads,
            initializer=initializer,
            initargs=([params, refs_damaged, refs_tax],),
        )

        data = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(
                        get_alns,
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

    p.close()
    p.join()
    # prof.disable()
    # # print profiling output
    # stats = pstats.Stats(prof).sort_stats("tottime")
    # stats.print_stats(10)
    log.info(f"Merging {len(ref_chunks)} chunks...")
    data = merge_dicts(data)  # top 10 rows
    return data
