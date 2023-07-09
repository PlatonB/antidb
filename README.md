# antidb
## Quick start
```
from pprint import pprint
from antidb import Idx, Prs, count_exec_time

dbsnp_vcf_path = '/mnt/Storage/databases/dbSNP_platon/GCF_000001405.40.vcf'
dbsnp_idx_prefix = 'all_rsids'

dbsnp_idx = Idx(dbsnp_vcf_path,
                dbsnp_idx_prefix)


@dbsnp_idx.idx
def get_rsid(dbsnp_zst_line):
    return dbsnp_zst_line.split('\t')[2]


get_rsid()

dbsnp_prs = Prs(dbsnp_vcf_path,
                dbsnp_idx_prefix)


@count_exec_time
def search_rsid_lines(dbsnp_prs):
    for dbsnp_zst_line in dbsnp_prs.prs(['rs1009150',
                                         'rs12044852',
                                         'rs4902496']):
        print(dbsnp_zst_line)


pprint(search_rsid_lines(dbsnp_prs))
```

## App example
### Bioinformatic annotator template
```
# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import json
import os
from argparse import ArgumentParser
from datetime import datetime
from antidb import Idx, Prs, count_exec_time

arg_parser = ArgumentParser()
arg_parser.add_argument('-S', '--ann-file-path', required=True, metavar='str', dest='ann_file_path', type=str,
                        help='Path to table with rsIDs column (uncompressed)')
arg_parser.add_argument('-D', '--dbsnp-file-path', required=True, metavar='str', dest='dbsnp_file_path', type=str,
                        help='Path to official dbSNP VCF (uncompressed or compressed via Seekable zstd)')
arg_parser.add_argument('-R', '--rsmerged-file-path', required=True, metavar='str', dest='rsmerged_file_path', type=str,
                        help='Path to official refsnp-merged JSON (uncompressed or compressed via Seekable zstd)')
arg_parser.add_argument('-T', '--trg-dir-path', required=True, metavar='str', dest='trg_dir_path', type=str,
                        help='Path to directory for results')
arg_parser.add_argument('-c', '--rsids-col-num', metavar='1', default=1, dest='rsids_col_num', type=int,
                        help='rsIDs-column number in source table')
args = arg_parser.parse_args()

dbsnp_idx = Idx(args.dbsnp_file_path,
                'rsids__gnomad_cln')


@dbsnp_idx.idx
def parse_dbsnp_line(dbsnp_zst_line):
    if 'GnomAD' in dbsnp_zst_line \
            and 'CLN' in dbsnp_zst_line:
        return dbsnp_zst_line.split('\t')[2]
    return None


rsmerged_idx = Idx(args.rsmerged_file_path,
                   'rsids')


@rsmerged_idx.idx
def parse_rsmerged_line(rsmerged_zst_line):
    rsmerged_zst_obj = json.loads(rsmerged_zst_line)
    rsids = [rsmerged_zst_obj['refsnp_id']] + \
        rsmerged_zst_obj['merged_snapshot_data']['merged_into']
    return rsids


parse_dbsnp_line()
parse_rsmerged_line()

perf = dbsnp_idx.perf + rsmerged_idx.perf

dbsnp_prs = Prs(args.dbsnp_file_path,
                'rsids__gnomad_cln')
rsmerged_prs = Prs(args.rsmerged_file_path,
                   'rsids')


@count_exec_time
def ann(args, res_files_crt_time, dbsnp_prs, rsmerged_prs):
    trg_file_path = os.path.join(args.trg_dir_path,
                                 f'ann_res_{res_files_crt_time}.txt')
    dump_file_path = os.path.join(args.trg_dir_path,
                                  f'ann_dump_{res_files_crt_time}.txt')
    with open(args.ann_file_path) as ann_file_opened:
        with open(trg_file_path, 'w') as trg_file_opened:
            with open(dump_file_path, 'w') as dump_file_opened:
                for ann_file_line in ann_file_opened:
                    if ann_file_line.startswith('#'):
                        continue
                    empty_res = True
                    ann_file_line = ann_file_line.rstrip()
                    ann_rsid = ann_file_line.split('\t')[args.rsids_col_num - 1]
                    for dbsnp_zst_line in dbsnp_prs.prs(ann_rsid):
                        empty_res = False
                        trg_file_opened.write(ann_file_line +
                                              dbsnp_zst_line)
                    if empty_res:
                        for rsmerged_zst_line in rsmerged_prs.prs(ann_rsid):
                            ann_rsid_syns = parse_rsmerged_line.__wrapped__(rsmerged_zst_line)
                            for dbsnp_zst_line in dbsnp_prs.prs(ann_rsid_syns):
                                empty_res = False
                                trg_file_opened.write(ann_file_line +
                                                      dbsnp_zst_line)
                            if not empty_res:
                                break
                        else:
                            dump_file_opened.write(ann_file_line + '\n')


res_files_crt_time = datetime.now()
perf.append(ann(args,
                res_files_crt_time,
                dbsnp_prs,
                rsmerged_prs))

perf_file_path = os.path.join(args.trg_dir_path,
                              f'ann_perf_{res_files_crt_time}.txt')
with open(perf_file_path, 'w') as perf_file_opened:
    for func_name, exec_time in perf:
        perf_file_opened.write(f'{func_name}\t{exec_time}\n')
```
