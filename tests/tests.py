# autopep8: off
import sys
import unittest
import os
import json
from pathlib import PurePath
from decimal import Decimal
import pyzstd
sys.dont_write_bytecode = True
par_dir_path = PurePath(__file__).parent.parent
sys.path.append(par_dir_path.joinpath('src/antidb').as_posix())
from antidb import (Idx, Prs)
# autopep8: on

__version__ = 'v1.0.0'
__authors__ = ['Platon Bykadorov (platon.work@gmail.com), 2023']


def remove_old_files(idx_obj):
    if os.path.exists(idx_obj.db_zst_path):
        os.remove(idx_obj.db_zst_path)
    if os.path.exists(idx_obj.full_idx_path):
        os.remove(idx_obj.full_idx_path)
    if os.path.exists(idx_obj.full_idx_tmp_path):
        os.remove(idx_obj.full_idx_tmp_path)
    if os.path.exists(idx_obj.full_idx_tmp_srtd_path):
        os.remove(idx_obj.full_idx_tmp_srtd_path)
    if os.path.exists(idx_obj.mem_idx_path):
        os.remove(idx_obj.mem_idx_path)


def remove_new_files(idx_obj):
    os.remove(idx_obj.db_zst_path)
    os.remove(idx_obj.full_idx_path)
    os.remove(idx_obj.mem_idx_path)


class RefsnpChrmtJsonTest(unittest.TestCase):
    mt_jsonbz2_url = 'https://ftp.ncbi.nih.gov/snp/archive/b156/JSON/refsnp-chrMT.json.bz2'
    mt_json_path = PurePath(PurePath(__file__).parent,
                            PurePath(mt_jsonbz2_url).name[:-4]).as_posix()

    def test_mt(self):
        if not os.path.exists(self.mt_json_path):
            os.system(f'''
                      wget -q -O - {self.mt_jsonbz2_url} |
                      bzip2 -d > {self.mt_json_path}''')
        mt_idx = Idx(self.mt_json_path,
                     'rsids')
        remove_old_files(mt_idx)

        @mt_idx.idx
        def parse_mt_line(mt_zst_line):
            return json.loads(mt_zst_line)['refsnp_id']

        parse_mt_line()
        with pyzstd.open(mt_idx.full_idx_path,
                         mode='rt') as full_idx_opened:
            full_idx_lines_cnt = 0
            for full_idx_line in full_idx_opened:
                full_idx_lines_cnt += 1
            self.assertEqual(full_idx_lines_cnt, 7891)
        with pyzstd.open(mt_idx.mem_idx_path,
                         mode='rt') as mem_idx_opened:
            mem_idx_lines_cnt = 0
            for mem_idx_line in mem_idx_opened:
                mem_idx_lines_cnt += 1
            self.assertEqual(mem_idx_lines_cnt, 8)
        mt_prs = Prs(self.mt_json_path,
                     'rsids')
        mt_prs_res = []
        for mt_zst_line in mt_prs.prs(['2124599696',
                                       8936,
                                       Decimal('368463610'),
                                       'nonexistent']):
            mt_prs_res.append(json.loads(mt_zst_line)['refsnp_id'])
        self.assertEqual(mt_prs_res,
                         ['2124599696',
                          '8936',
                          '368463610'])
        remove_new_files(mt_idx)

    def test_mt_clinical(self):
        if not os.path.exists(self.mt_json_path):
            os.system(f'''
                      wget -q -O - {self.mt_jsonbz2_url} |
                      bzip2 -d > {self.mt_json_path}''')
        mt_cln_idx = Idx(self.mt_json_path,
                         'rsids_cln')
        remove_old_files(mt_cln_idx)

        @mt_cln_idx.idx
        def parse_mt_cln_line(mt_zst_line):
            mt_zst_obj = json.loads(mt_zst_line)
            if len(mt_zst_obj['primary_snapshot_data']['allele_annotations'][0]['clinical']) > 0:
                return mt_zst_obj['refsnp_id']
            return None

        parse_mt_cln_line()
        with pyzstd.open(mt_cln_idx.full_idx_path,
                         mode='rt') as full_idx_opened:
            full_idx_lines_cnt = 0
            for full_idx_line in full_idx_opened:
                full_idx_lines_cnt += 1
            self.assertEqual(full_idx_lines_cnt, 2)
        with pyzstd.open(mt_cln_idx.mem_idx_path,
                         mode='rt') as mem_idx_opened:
            mem_idx_lines_cnt = 0
            for mem_idx_line in mem_idx_opened:
                mem_idx_lines_cnt += 1
            self.assertEqual(mem_idx_lines_cnt, 1)
        mt_cln_prs = Prs(self.mt_json_path,
                         'rsids_cln')
        mt_cln_prs_res = []
        for mt_zst_line in mt_cln_prs.prs(['qwerty',
                                           '2001030',
                                           'uiopas',
                                           '1556422499',
                                           'dfghjk']):
            mt_cln_prs_res.append(json.loads(mt_zst_line)['refsnp_id'])
        self.assertEqual(mt_cln_prs_res,
                         ['2001030',
                          '1556422499'])
        remove_new_files(mt_cln_idx)


if __name__ == "__main__":
    unittest.main()
