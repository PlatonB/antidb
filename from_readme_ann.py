# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import json
import os
from argparse import ArgumentParser
from datetime import datetime
from antidb.idx import (Idx,
                        count_exec_time)
from antidb.prs import Prs

if __name__ == '__main__':
    __version__ = 'v1.3.2'


def prs_dbsnp_line(dbsnp_zst_line: str) -> str | None:
    if 'GnomAD' in dbsnp_zst_line \
            and 'CLN' in dbsnp_zst_line:
        return dbsnp_zst_line.split('\t')[2]
    return None


def prs_rsmerged_line(rsmerged_zst_line: str) -> tuple:
    rsmerged_zst_obj = json.loads(rsmerged_zst_line)
    rsids = tuple(map(lambda rsid: f'rs{rsid}',
                      ([rsmerged_zst_obj['refsnp_id']] +
                      rsmerged_zst_obj['merged_snapshot_data']['merged_into'])))
    return rsids


def query_rsid(rsid, dbsnp_prs_obj,
               rsmerged_prs_obj, prs_rsmerged_line) -> str:
    for dbsnp_zst_line in dbsnp_prs_obj.eq(rsid):
        return dbsnp_zst_line
    for rsmerged_zst_line in rsmerged_prs_obj.eq(rsid):
        rsid_syns = prs_rsmerged_line(rsmerged_zst_line)
        for dbsnp_zst_line in dbsnp_prs_obj.eq(*rsid_syns):
            return dbsnp_zst_line
    return None


@count_exec_time
def ann(args, res_files_crt_time,
        dbsnp_prs_obj, rsmerged_prs_obj,
        prs_rsmerged_line) -> None:
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
                    ann_file_line = ann_file_line.rstrip()
                    ann_rsid = ann_file_line.split('\t')[args.rsids_col_num - 1]
                    dbsnp_zst_line = query_rsid(ann_rsid,
                                                dbsnp_prs_obj,
                                                rsmerged_prs_obj,
                                                prs_rsmerged_line)
                    if dbsnp_zst_line:
                        trg_file_opened.write(ann_file_line + '\t' +
                                              dbsnp_zst_line)
                    else:
                        dump_file_opened.write(ann_file_line + '\n')


if __name__ == '__main__':
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
                    'rsids__gnomad_cln',
                    prs_dbsnp_line,
                    lambda rsid: rsid)
    dbsnp_idx.idx()
    rsmerged_idx = Idx(args.rsmerged_file_path,
                       'rsids',
                       prs_rsmerged_line,
                       lambda rsid: rsid)
    rsmerged_idx.idx()
    perf = {'dbsnp_idx': dbsnp_idx.perf,
            'rsmerged_idx': rsmerged_idx.perf}
    dbsnp_prs_obj = Prs(args.dbsnp_file_path,
                        'rsids__gnomad_cln',
                        lambda rsid: rsid)
    rsmerged_prs_obj = Prs(args.rsmerged_file_path,
                           'rsids',
                           lambda rsid: rsid)
    res_files_crt_time = datetime.now()
    perf['ann'] = ann(args,
                      res_files_crt_time,
                      dbsnp_prs_obj,
                      rsmerged_prs_obj,
                      prs_rsmerged_line)[2]
    perf_file_path = os.path.join(args.trg_dir_path,
                                  f'ann_perf_{res_files_crt_time}.json')
    with open(perf_file_path, 'w') as perf_file_opened:
        json.dump(perf, perf_file_opened, indent=4)
