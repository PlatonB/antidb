# pydx
## Super quick start
```
from pydx import Idx, Prs

dbsnp_vcf_path = '/mnt/Storage/databases/dbSNP_platon/GCF_000001405.40.vcf'
idx_prefix = 'allrsids'

idx = Idx(dbsnp_vcf_path,
          idx_prefix)


@idx.idx
def get_id(dbsnp_zst_line):
    return dbsnp_zst_line.split('\t')[2]


get_id()

prs = Prs(dbsnp_vcf_path,
          idx_prefix)

for rs_id in ['rs1009150',
              'rs12044852',
              'rs4902496']:
    for dbsnp_zst_line in prs.prs(rs_id):
        print(dbsnp_zst_line)
```

## Quick start
### your_tool_cli.py
```
from argparse import ArgumentParser


def add_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('-D', '--db-file-path', required=True, metavar='str', dest='db_file_path', type=str,
                            help='Path to indexable multi-line text file (Uncompressed or compressed via Seekable zstd)')
    arg_parser.add_argument('-P', '--idx-prefix', required=True, metavar='str', dest='idx_prefix', type=str,
                            help='Unique supplement to index name that allows to distinguish it from other indexes names')
    arg_parser.add_argument('-T', '--trg-dir-path', required=True, metavar='str', dest='trg_dir_path', type=str,
                            help='Path to directory for results')
    arg_parser.add_argument('-b', '--bench-file-path', metavar='None', dest='bench_file_path', type=str,
                            help='Path for index creation performance measurement results file')
    return arg_parser.parse_args()
```
### your_tool.py
```
import json
import os
from your_tool_cli import add_args
from pydx import Idx, Prs

# Getting from CLI path to
# DB-file, path to directory
# for target files, and prefix
# that makes index names unique.
args = add_args()

# Initialization of a
# class containing methods,
# whose sequential call can
# prepare DB-file indexes.
idx = Idx(args.db_file_path,
          args.idx_prefix)


@idx.idx
def rsmerged_line_parser(db_zst_line):
    '''
    An exemplary algorithm that pulls indexable rsIDs from each line
    of the refsnp-merged.json table. This table can be downloaded from
    https://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-merged.json.bz2.
    Don't forget to decompress it: bzip2 -d /path/to/refsnp-merged.json.bz2
    '''
    db_zst_obj = json.loads(db_zst_line)
    rsmerged_line_parser_out = [db_zst_obj['refsnp_id']] + \
        db_zst_obj['merged_snapshot_data']['merged_into']
    return rsmerged_line_parser_out


# Compressing DB-file and
# creating two indexes. If
# either of these files already
# exists, the corresponding
# step is skipped. Full-index
# contains DB-file lines elements
# found by the decorated function
# and start positions of these
# lines. Mem-index contains every
# N-th element previously indexed
# and the start position of the
# corresponding full-index line.
# Full-index is usually very
# big, while mem-index does
# not overload RAM at all.
rsmerged_line_parser()

# Optional. Writing into separate file the results
# of index creation performance measurement.
if args.bench_file_path and not os.path.exists(args.bench_file_path):
    with open(args.bench_file_path, 'w') as bench_file_opened:
        for file_time in idx.bench.items():
            bench_file_opened.write(f'{file_time[0]}\t{file_time[1]}\n')

# Initialize the class
# containing the parser
# method. The latter can
# create a generator that
# one by one returns DB-file
# strings with found values.
prs = Prs(args.db_file_path,
          args.idx_prefix)

# Instant search for your values in the
# indexed DB file. There are two search
# strategies: 1. you create an array
# of values and pass it to the parser;
# 2. you create the parser multiple times,
# passing to it a single value each
# time. Values can be of any type: the
# parser itself converts them to str.
trg_file1_path = os.path.join(args.trg_dir_path,
                              'out_1.txt')
with open(trg_file1_path, 'w') as trg_file1_opened:
    for db_zst_found_line in prs.prs(['332',
                                      '1192046386']):
        trg_file1_opened.write(db_zst_found_line)
trg_file2_path = os.path.join(args.trg_dir_path,
                              'out_2.txt')
with open(trg_file2_path, 'w') as trg_file2_opened:
    for rs_id in [332, 1192046386]:
        for db_zst_found_line in prs.prs(rs_id):
            trg_file2_opened.write(db_zst_found_line)
```
