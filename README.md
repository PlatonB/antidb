# antidb
## Quick start
```
pip3 install antidb
```
```
from antidb.idx import (Idx,
                        count_exec_time)
from antidb.prs import Prs

__version__ = 'v1.4.0'

dbsnp_vcf_path = '/path/to/GCF_000001405.40[.zst]'
dbsnp_idx_prefix = 'all_rsids'
dbsnp_idx = Idx(db_file_path=dbsnp_vcf_path,
                adb_name_prefix=dbsnp_idx_prefix,
                db_line_prs=(lambda dbsnp_zst_line:
                             dbsnp_zst_line.split('\t')[2]),
                adb_srt_rule=lambda rsid: rsid)
dbsnp_idx.idx()
dbsnp_prs = Prs(db_file_path=dbsnp_vcf_path,
                adb_name_prefix=dbsnp_idx_prefix,
                adb_srt_rule=lambda rsid: rsid)


@count_exec_time
def get_rsid_lines(dbsnp_prs: Prs):
    for dbsnp_vcfzst_line in dbsnp_prs.eq('rs1009150',
                                          'rs12044852',
                                          'rs4902496'):
        print(dbsnp_vcfzst_line)


print(get_rsid_lines(dbsnp_prs))
```
```
NC_000022.11	36306254	rs1009150	C	T	.	.	RS=1009150;dbSNPBuildID=86;SSR=0;GENEINFO=MYH9:4627;VC=SNV;PUB;INT;GNO;FREQ=1000Genomes:0.569,0.431|ALSPAC:0.2906,0.7094|Estonian:0.269,0.731|GENOME_DK:0.35,0.65|GnomAD:0.4415,0.5585|GoNL:0.3126,0.6874|HapMap:0.5881,0.4119|KOREAN:0.7334,0.2666|MGP:0.8652,0.1348|NorthernSweden:0.315,0.685|Qatari:0.5463,0.4537|SGDP_PRJ:0.2929,0.7071|Siberian:0.3043,0.6957|TOMMO:0.7117,0.2883|TOPMED:0.4596,0.5404|TWINSUK:0.2869,0.7131|dbGaP_PopFreq:0.3304,0.6696;COMMON;CLNVI=.,;CLNORIGIN=.,1;CLNSIG=.,2;CLNDISDB=.,MedGen:CN517202;CLNDN=.,not_provided;CLNREVSTAT=.,single;CLNACC=.,RCV001695529.1;CLNHGVS=NC_000022.11:g.36306254=,NC_000022.11:g.36306254C>T

NC_000001.11	116545157	rs12044852	C	A	.	.	RS=12044852;dbSNPBuildID=120;SSR=0;GENEINFO=CD58:965|LOC105378925:105378925;VC=SNV;PUB;INT;GNO;FREQ=1000Genomes:0.7473,0.2527|ALSPAC:0.8957,0.1043|Chileans:0.7396,0.2604|Estonian:0.9125,0.0875|GENOME_DK:0.875,0.125|GnomAD:0.8826,0.1174|GoNL:0.9078,0.09218|HapMap:0.787,0.213|KOREAN:0.3945,0.6055|Korea1K:0.3892,0.6108|NorthernSweden:0.895,0.105|PRJEB37584:0.439,0.561|Qatari:0.8704,0.1296|SGDP_PRJ:0.3373,0.6627|Siberian:0.3846,0.6154|TOMMO:0.4146,0.5854|TOPMED:0.8671,0.1329|TWINSUK:0.8972,0.1028|Vietnamese:0.4486,0.5514|dbGaP_PopFreq:0.8864,0.1136;COMMON

NC_000014.9	67588896	rs4902496	C	G,T	.	.	RS=4902496;dbSNPBuildID=111;SSR=0;GENEINFO=PIGH:5283|GPHN:10243|PLEKHH1:57475;VC=SNV;PUB;U3;INT;R3;GNO;FREQ=1000Genomes:0.3357,0.6643,.|ALSPAC:0.2019,0.7981,.|Estonian:0.1518,0.8482,.|GENOME_DK:0.125,0.875,.|GoNL:0.1703,0.8297,.|HapMap:0.3639,0.6361,.|KOREAN:0.3399,0.6601,.|MGP:0.3558,0.6442,.|NorthernSweden:0.1817,0.8183,.|Qatari:0.2176,0.7824,.|SGDP_PRJ:0.189,0.811,.|Siberian:0.1429,0.8571,.|TOMMO:0.2816,0.7184,.|TOPMED:0.285,0.715,.|TWINSUK:0.1888,0.8112,.|Vietnamese:0.4533,0.5467,.|dbGaP_PopFreq:0.2712,0.7288,0;COMMON

('get_rsid_lines', None, '0:00:00.015202')
```

## Features
- As in classical DBMSs, you spend time indexing once and then run queries in hundredths of a second.
- Designed on a laptop for laptops. It is unlikely to overflow RAM. There is such a risk when indexing, but certainly not when parsing.
- Instead of the typical database hidden in the system directory, you will see an index file neighboring your multiline text file. It's easy to publish them or save to a USB drive.
- Compared to _[tabix](https://www.htslib.org/doc/tabix.html)_, there is no need to sort the data yourself before indexing.
- You write the function for pulling indexable values yourself. This means complete freedom to choose what and how to index. Note that queried values must correspond (e.g., by data type) to the values returned by your indexing function.
- The sort key is also created by you. Just don't forget to consider the sort order when making queries.
- It is possible to create phantom indexes, i.e. when the values themselves are not physically present in the file, but they are present in the index. For example, the index may contain the lengths of reference allele sequences calculated for deletions.
- The _antidb_ syntax is extremely simple and doesn't require bulky API docs. Simply look at the example scripts/tools here.

## Query syntax
It is designed that _antidb_ supports only the simplest queries. A good work scenario is when you reduce the data by simple query to RAM-friendly sizes and post-process it in _pandas_ or something else.

`Prs.eq(*queries)`: creates a generator capable to return lines of indexed file containing element that exactly match your argument. Each argument is a separate query. If nothing matches the query, the generator will not throw an exception, but just not return anything.

`Prs.rng(query_start, query_end)`: creates a generator capable to return lines of indexed file containing elements in the range you specify. Performance note: queries covering a large quantity of lines may run slowly.

## App examples
### Bioinformatic annotator template
It would seem that finding rsIDs by rsIDs is easy. But, unlike genomic coordinates, rsIDs are quite often updated. Therefore, rsIDs should be queried by dbSNP, and in case of failure - by the source of rsID synonyms with further attempt to find a synonym again by dbSNP. This code demonstrates how _antidb_ helps quickly retrieve data from two sources, easily switching between them when needed.

```
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
```

#### Performance measurement results
##### Annotation of 2842 SNPs by dbSNP VCF and refsnp-merged JSON
The calculations were done on a Maibenben P687 laptop with an AMD Ryzen 7 8845HS processor and 32GB RAM. In `*_idx` is the time of indexing steps, in `ann` is the search time.

```
{
    "dbsnp_idx": [
        [
            "presrt_idxs",
            null,
            "0:41:02.303547"
        ],
        [
            "crt_adb",
            null,
            "0:00:03.331788"
        ]
    ],
    "rsmerged_idx": [
        [
            "presrt_idxs",
            null,
            "0:04:37.101617"
        ],
        [
            "crt_adb",
            null,
            "0:02:00.120555"
        ]
    ],
    "ann": "0:00:24.904676"
}
```