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
from dmg_reads.defaults import mdmg_header

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)
timestr = time.strftime("%Y%m%d-%H%M%S")


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


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


def get_ranks(parser, ranks, var):
    valid_ranks = ["all", "phylum", "class", "order", "family", "genus"]
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
    "metaDMG_filter": {"Bayesian_D_max": 0.1, "Bayesian_z": 2.5},
    "prefix": None,
}

help_msg = {
    "metaDMG_results": f"A file from metaDMG ran in local mode",
    "metaDMG_filter": f"Which filter to use for metaDMG results",
    "bam": f"The BAM file used to generate the metaDMG results",
    "prefix": "Prefix used for the output files",
    "threads": f"Number of threads",
    "debug": f"Print debug messages",
    "version": f"Print program version",
}


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="A simple tool to extract damaged reads from BAM files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-r",
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


def create_output_files(prefix, bam):
    if prefix is None:
        prefix = bam.split(".")[0]
    # create output files
    out_files = {
        "stats": f"{prefix}_stats.tsv.gz",
        "fastq_damaged": f"{prefix}.damaged.fastq.gz",
        "fastq_nondamaged": f"{prefix}.nondamaged.fastq.gz",
    }
    return out_files
