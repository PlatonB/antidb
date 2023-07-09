# antidb
## Quick start
```
from antidb import (Idx,
                    Prs,
                    count_exec_time)

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
def get_rsid_lines(dbsnp_prs):
    for dbsnp_zst_line in dbsnp_prs.prs(['rs1009150',
                                         'rs12044852',
                                         'rs4902496']):
        print(dbsnp_zst_line)


print(get_rsid_lines(dbsnp_prs))
```
```
NC_000022.11    36306254        rs1009150       C       T       .       .       RS=1009150;dbSNPBuildID=86;SSR=0;GENEINFO=MYH9:4627;VC=SNV;PUB;INT;GNO;FREQ=1000Genomes:0.569,0.431|ALSPAC:0.2906,0.7094|Estonian:0.269,0.731|GENOME_DK:0.35,0.65|GnomAD:0.4415,0.5585|GoNL:0.3126,0.6874|HapMap:0.5881,0.4119|KOREAN:0.7334,0.2666|MGP:0.8652,0.1348|NorthernSweden:0.315,0.685|Qatari:0.5463,0.4537|SGDP_PRJ:0.2929,0.7071|Siberian:0.3043,0.6957|TOMMO:0.7117,0.2883|TOPMED:0.4596,0.5404|TWINSUK:0.2869,0.7131|dbGaP_PopFreq:0.3304,0.6696;COMMON;CLNVI=.,;CLNORIGIN=.,1;CLNSIG=.,2;CLNDISDB=.,MedGen:CN517202;CLNDN=.,not_provided;CLNREVSTAT=.,single;CLNACC=.,RCV001695529.1;CLNHGVS=NC_000022.11:g.36306254=,NC_000022.11:g.36306254C>T

NC_000001.11    116545157       rs12044852      C       A       .       .       RS=12044852;dbSNPBuildID=120;SSR=0;GENEINFO=CD58:965|LOC105378925:105378925;VC=SNV;PUB;INT;GNO;FREQ=1000Genomes:0.7473,0.2527|ALSPAC:0.8957,0.1043|Chileans:0.7396,0.2604|Estonian:0.9125,0.0875|GENOME_DK:0.875,0.125|GnomAD:0.8826,0.1174|GoNL:0.9078,0.09218|HapMap:0.787,0.213|KOREAN:0.3945,0.6055|Korea1K:0.3892,0.6108|NorthernSweden:0.895,0.105|PRJEB37584:0.439,0.561|Qatari:0.8704,0.1296|SGDP_PRJ:0.3373,0.6627|Siberian:0.3846,0.6154|TOMMO:0.4146,0.5854|TOPMED:0.8671,0.1329|TWINSUK:0.8972,0.1028|Vietnamese:0.4486,0.5514|dbGaP_PopFreq:0.8864,0.1136;COMMON

NC_000014.9     67588896        rs4902496       C       G,T     .       .       RS=4902496;dbSNPBuildID=111;SSR=0;GENEINFO=PIGH:5283|GPHN:10243|PLEKHH1:57475;VC=SNV;PUB;U3;INT;R3;GNO;FREQ=1000Genomes:0.3357,0.6643,.|ALSPAC:0.2019,0.7981,.|Estonian:0.1518,0.8482,.|GENOME_DK:0.125,0.875,.|GoNL:0.1703,0.8297,.|HapMap:0.3639,0.6361,.|KOREAN:0.3399,0.6601,.|MGP:0.3558,0.6442,.|NorthernSweden:0.1817,0.8183,.|Qatari:0.2176,0.7824,.|SGDP_PRJ:0.189,0.811,.|Siberian:0.1429,0.8571,.|TOMMO:0.2816,0.7184,.|TOPMED:0.285,0.715,.|TWINSUK:0.1888,0.8112,.|Vietnamese:0.4533,0.5467,.|dbGaP_PopFreq:0.2712,0.7288,0;COMMON

('get_rsid_lines', '0:00:00.007858')
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
from antidb import (Idx,
                    Prs,
                    count_exec_time)

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

perf = {'dbsnp_idx': dbsnp_idx.perf,
        'rsmerged_idx': rsmerged_idx.perf}

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

perf['ann'] = ann(args,
                  res_files_crt_time,
                  dbsnp_prs,
                  rsmerged_prs)[1]

perf_file_path = os.path.join(args.trg_dir_path,
                              f'ann_perf_{res_files_crt_time}.json')
with open(perf_file_path, 'w') as perf_file_opened:
    json.dump(perf, perf_file_opened, indent=4)
```

#### Performance measurement results
##### ann_perf_2023-07-09 20:45:36.102376.json
- `dbsnp_idx` - indexing `GnomAD`- and `CLN`-containing lines of dbSNP VCF;
  - `crt_db_zst` - compressing indexable file (output is further called "DB");
  - `crt_full_idx_tmp` - indexing DB (output is further called "temporary full index");
  - `crt_full_idx_tmp_srtd` - sorting temporary full index by indexed DB elements;
  - `crt_full_idx` - compressing sorted temporary full index (output is further called "full index");
  - `crt_mem_idx` - selective indexing of full index;
- `rsmerged_idx` - indexing all lines of rsmerged JSON;
  - <...>
- `ann` - querying 2842 rsIDs by indexed dbSNP VCF and indexed rsmerged JSON.
```
{
    "dbsnp_idx": [
        [
            "crt_db_zst",
            "0:39:02.127938"
        ],
        [
            "crt_full_idx_tmp",
            "1:06:13.698458"
        ],
        [
            "crt_full_idx_tmp_srtd",
            "0:00:00.928633"
        ],
        [
            "crt_full_idx",
            "0:00:00.577710"
        ],
        [
            "crt_mem_idx",
            "0:00:00.280014"
        ]
    ],
    "rsmerged_idx": [
        [
            "crt_db_zst",
            "0:02:44.068920"
        ],
        [
            "crt_full_idx_tmp",
            "0:04:43.153807"
        ],
        [
            "crt_full_idx_tmp_srtd",
            "0:00:30.015826"
        ],
        [
            "crt_full_idx",
            "0:00:17.204649"
        ],
        [
            "crt_mem_idx",
            "0:00:08.811190"
        ]
    ],
    "ann": "0:00:06.995505"
}
```
