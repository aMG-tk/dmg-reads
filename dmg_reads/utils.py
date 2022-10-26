import argparse
import sys
import gzip
import os
import logging
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull
from dmg_reads import __version__
import time
import json
from dmg_reads.defaults import mdmg_header, valid_ranks
from collections import defaultdict
import re
from itertools import chain

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)
timestr = time.strftime("%Y%m%d-%H%M%S")


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


# From: https://note.nkmk.me/en/python-check-int-float/
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


# function to check if the input value has K, M or G suffix in it
def check_suffix(val, parser, var):
    if var == "--scale":
        units = ["K", "M"]
    else:
        units = ["K", "M", "G"]
    unit = val[-1]
    value = int(val[:-1])

    if is_integer(value) & (unit in units) & (value > 0):
        if var == "--scale":
            if unit == "K":
                val = value * 1000
            elif unit == "M":
                val = value * 1000000
            elif unit == "G":
                val = value * 1000000000
            return val
        else:
            return val
    else:
        parser.error(
            "argument %s: Invalid value %s. Has to be an integer larger than 0 with the following suffix K, M or G"
            % (var, val)
        )


# From https://stackoverflow.com/a/59617044/15704171
def convert_list_to_str(lst):
    n = len(lst)
    if not n:
        return ""
    if n == 1:
        return lst[0]
    return ", ".join(lst[:-1]) + f" or {lst[-1]}"


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


def check_values(val, minval, maxval, parser, var):
    value = float(val)
    if value < minval or value > maxval:
        parser.error(
            "argument %s: Invalid value %s. Range has to be between %s and %s!"
            % (
                var,
                value,
                minval,
                maxval,
            )
        )
    return value


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg, var):
    if not os.path.exists(arg):
        parser.error("argument %s: The file %s does not exist!" % (var, arg))
    else:
        return arg


def is_valid_filter(parser, arg, var):
    arg = json.loads(arg)
    # check if the dictionary keys are in the mdmg header list
    for key in arg.keys():
        if key not in mdmg_header:
            parser.error(
                f"argument {var}: Invalid value {key}.\n"
                f"Valid values are: {convert_list_to_str(mdmg_header)}"
            )

    return arg


def is_valid_rank(parser, arg, var):
    arg = json.loads(arg)
    for key in arg.keys():
        if key not in valid_ranks.keys():
            parser.error(
                f"argument {var}: Invalid value {key}.\n"
                f"Valid values are: {convert_list_to_str(valid_ranks)}"
            )
    return arg


def get_ranks(parser, ranks, var):
    valid_ranks = [
        "domain",
        "kingdom",
        "lineage",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    ranks = ranks.split(",")
    # check if ranks are valid
    for rank in ranks:
        if rank not in valid_ranks:
            parser.error(
                f"argument {var}: Invalid value {rank}.\Rank has to be one of {convert_list_to_str(valid_ranks)}"
            )
        if rank == "all":
            ranks = valid_ranks[1:]
    return ranks


defaults = {
    "metaDMG_filter": {"damage": 0.1, "significance": 2},
    "prefix": None,
    "sort_memory": "1G",
    "threads": 1,
}

help_msg = {
    "metaDMG_results": "A file from metaDMG ran in local mode",
    "metaDMG_filter": "Which filter to use for metaDMG results",
    "bam": "The BAM file used to generate the metaDMG results",
    "prefix": "Prefix used for the output files",
    "taxonomy_file": "A file containing the taxonomy of the BAM references in the format d__;p__;c__;o__;f__;g__;s__.",
    "rank": "Which taxonomic group and rank we want to get the reads extracted.",
    "combine": "If set, the reads damaged and non-damaged will be combined in one fastq file",
    "only_damaged": "If set, only the reads damaged will be extracted",
    "sort_memory": "Set maximum memory per thread for sorting; suffix K/M/G recognized",
    "threads": "Number of threads",
    "debug": "Print debug messages",
    "version": "Print program version",
}


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="A simple tool to extract damaged reads from BAM files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-m",
        "--metaDMG-results",
        type=lambda x: is_valid_file(parser, x, "--metaDMG-results"),
        dest="metaDMG_results",
        help=help_msg["metaDMG_results"],
        required=True,
    )
    parser.add_argument(
        "-f",
        "--metaDMG-filter",
        type=lambda x: is_valid_filter(parser, x, "--metaDMG-filter"),
        dest="metaDMG_filter",
        help=help_msg["metaDMG_filter"],
        required=True,
    )
    parser.add_argument(
        "-b",
        "--bam",
        type=lambda x: is_valid_file(parser, x, "--bam"),
        dest="bam",
        help=help_msg["bam"],
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=defaults["prefix"],
        dest="prefix",
        help=help_msg["prefix"],
    )
    parser.add_argument(
        "--combine", dest="combine", action="store_true", help=help_msg["combine"]
    )
    parser.add_argument(
        "--only-damaged",
        dest="only_damaged",
        action="store_true",
        help=help_msg["only_damaged"],
    )
    parser.add_argument(
        "-T",
        "--taxonomy-file",
        type=lambda x: is_valid_file(parser, x, "---taxonomy-file"),
        dest="taxonomy_file",
        help=help_msg["taxonomy_file"],
    )
    parser.add_argument(
        "-r",
        "--rank",
        type=lambda x: is_valid_rank(parser, x, "--rank"),
        dest="rank",
        help=help_msg["rank"],
    )
    parser.add_argument(
        "-M",
        "--sort-memory",
        type=lambda x: check_suffix(x, parser=parser, var="--sort-memory"),
        default=defaults["sort_memory"],
        dest="sort_memory",
        help=help_msg["sort_memory"],
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--threads")
        ),
        dest="threads",
        default=1,
        help=help_msg["threads"],
    )
    parser.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )
    args = parser.parse_args(None if sys.argv[1:] else ["-h"])
    return args


@contextmanager
def suppress_stdout():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


def create_output_files(prefix, bam, taxon=None, combined=False):
    if prefix is None:
        prefix = bam.split(".")[0]
    # create output files
    if taxon is None:
        if combined:
            out_files = {
                "fastq_combined": f"{prefix}.fastq.gz",
            }
        else:
            out_files = {
                "fastq_damaged": f"{prefix}.damaged.fastq.gz",
                "fastq_nondamaged": f"{prefix}.non-damaged.fastq.gz",
                "fastq_multi": f"{prefix}.multi.fastq.gz",
                "fastq_combined": f"{prefix}.fastq.gz",
            }
    else:
        if combined:
            out_files = defaultdict()
            for k, v in taxon.items():
                for i in v:
                    i = re.sub("[^0-9a-zA-Z]+", "_", i)
                    out_files[f"fastq_{k}{i}"] = f"{prefix}.{k}{i}.fastq.gz"
        else:
            out_files = defaultdict()
            for k, v in taxon.items():
                for i in v:
                    i = re.sub("[^0-9a-zA-Z]+", "_", i)
                    out_files[f"fastq_damaged_{k}{i}"] = f"{prefix}.{k}{i}.damaged.fastq.gz"
                    out_files[
                        f"fastq_nondamaged_{k}{i}"
                    ] = f"{prefix}.{k}{i}.non-damaged.fastq.gz"
                    out_files[f"fastq_multi_{k}{i}"] = f"{prefix}.{k}{i}.multi.fastq.gz"
                    out_files[f"fastq_combined_{k}{i}"] = f"{prefix}.{k}{i}.fastq.gz"
    return out_files


# From https://stackoverflow.com/a/61436083
def splitkeep(s, delimiter):
    split = s.split(delimiter)
    return [substr + delimiter for substr in split[:-1]] + [split[-1]]


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))
