# autopep8: off
import sys
import unittest
import os
import json
from pathlib import PurePath
from decimal import Decimal
from random import shuffle
import pyzstd
sys.dont_write_bytecode = True
par_dir_path = PurePath(__file__).parent.parent
sys.path.append(par_dir_path.joinpath('src/antidb').as_posix())
from antidb import (Idx,
                    Prs)
from antisrt import (SrtRules,
                     Srt)
# autopep8: on

__version__ = 'v2.2.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023'}]


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


class AntidbTests(unittest.TestCase):
    mt_jsonbz2_url = 'https://ftp.ncbi.nih.gov/snp/archive/b156/JSON/refsnp-chrMT.json.bz2'
    mt_json_path = PurePath(PurePath(__file__).parent,
                            PurePath(mt_jsonbz2_url).name[:-4]).as_posix()
    if not os.path.exists(mt_json_path):
        os.system(f'''
                  wget -q -O - {mt_jsonbz2_url} |
                  bzip2 -d > {mt_json_path}''')
    with open(mt_json_path) as mt_json_opened:
        mt_json_content = mt_json_opened.read()
    os.remove(mt_json_path)

    def test_mt(self):
        if not os.path.exists(self.mt_json_path):
            with open(self.mt_json_path, 'w') as mt_json_opened:
                mt_json_opened.write(self.mt_json_content)
        mt_idx = Idx(self.mt_json_path,
                     'rsids')
        remove_old_files(mt_idx)
        mt_idx.idx(lambda mt_zst_line:
                   json.loads(mt_zst_line)['refsnp_id'])
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
            self.assertEqual(mem_idx_lines_cnt, 9)
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

    def test_mt_by_mt(self):
        if not os.path.exists(self.mt_json_path):
            with open(self.mt_json_path, 'w') as mt_json_opened:
                mt_json_opened.write(self.mt_json_content)
        mt_idx = Idx(self.mt_json_path,
                     'rsids',
                     compr_frame_size=1024,
                     compr_chunk_size=1024,
                     compr_chunk_elems_quan=10,
                     unidx_lines_quan=10)
        remove_old_files(mt_idx)

        def parse_mt_line(mt_zst_line):
            return json.loads(mt_zst_line)['refsnp_id']

        mt_idx.idx(parse_mt_line)
        with pyzstd.open(mt_idx.mem_idx_path,
                         mode='rt') as mem_idx_opened:
            unidx_lines_quan = int(mem_idx_opened.readline().rstrip().split('=')[1])
            self.assertEqual(unidx_lines_quan, 10)
        mt_prs = Prs(self.mt_json_path,
                     'rsids')
        with pyzstd.open(mt_prs.db_zst_path,
                         mode='rt') as mt_zst_opened:
            mt_zst_rsids = [json.loads(mt_zst_line)['refsnp_id']
                            for mt_zst_line in mt_zst_opened]
        mt_prs_res = [parse_mt_line(mt_zst_line)
                      for mt_zst_line in mt_prs.prs(mt_zst_rsids)]
        self.assertEqual(mt_zst_rsids,
                         mt_prs_res)
        mt_zst_rsids_mxd = mt_zst_rsids[:]
        shuffle(mt_zst_rsids_mxd)
        mt_prs_res_mxd = [parse_mt_line(mt_zst_line)
                          for mt_zst_line in mt_prs.prs(mt_zst_rsids_mxd)]
        self.assertNotEqual(mt_zst_rsids_mxd,
                            mt_zst_rsids)
        self.assertNotEqual(mt_prs_res_mxd,
                            mt_prs_res)
        self.assertEqual(len(mt_zst_rsids_mxd),
                         len(mt_zst_rsids))
        self.assertEqual(len(mt_prs_res_mxd),
                         len(mt_prs_res))
        self.assertEqual(len(mt_prs_res_mxd),
                         len(mt_zst_rsids))
        remove_new_files(mt_idx)

    def test_mt_clinical(self):
        if not os.path.exists(self.mt_json_path):
            with open(self.mt_json_path, 'w') as mt_json_opened:
                mt_json_opened.write(self.mt_json_content)
        mt_cln_idx = Idx(self.mt_json_path,
                         'rsids_cln')
        remove_old_files(mt_cln_idx)

        def parse_mt_cln_line(mt_zst_line):
            mt_zst_obj = json.loads(mt_zst_line)
            if len(mt_zst_obj['primary_snapshot_data']['allele_annotations'][0]['clinical']) > 0:
                return mt_zst_obj['refsnp_id']
            return None

        mt_cln_idx.idx(parse_mt_cln_line)
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
            self.assertEqual(mem_idx_lines_cnt, 2)
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


class SrtRulesTests(unittest.TestCase):
    srt_rules = SrtRules()

    def test_natur_srt_rule(self):
        self.assertEqual(self.srt_rules.natur('10'),
                         [[10]])
        self.assertEqual(self.srt_rules.natur('01'),
                         [[1]])
        self.assertEqual(self.srt_rules.natur('rs10'),
                         [[float('+inf'), 'rs', 10]])
        self.assertEqual(self.srt_rules.natur('rs01'),
                         [[float('+inf'), 'rs', 1]])
        self.assertEqual(self.srt_rules.natur('val1.5'),
                         [[float('+inf'), 'val', 1.5]])
        self.assertEqual(self.srt_rules.natur('val1.05'),
                         [[float('+inf'), 'val', 1.05]])
        self.assertEqual(self.srt_rules.natur('val1.05suff'),
                         [[float('+inf'), 'val', 1.05, 'suff']])
        self.assertEqual(self.srt_rules.natur('val1,05'),
                         [[float('+inf'), 'val', 1, ',', 5]])
        self.assertEqual(self.srt_rules.natur('val1,05suff'),
                         [[float('+inf'), 'val', 1, ',', 5, 'suff']])
        self.assertEqual(self.srt_rules.natur('I2a2a1b2a2a2-ZS20'),
                         [[float('+inf'), 'I', 2, 'a', 2, 'a',
                           1, 'b', 2, 'a', 2, 'a', 2, '-ZS', 20]])
        self.assertEqual(self.srt_rules.natur('10\t11'),
                         [[10], [11]])
        self.assertEqual(self.srt_rules.natur('rs10\t11'),
                         [[float('+inf'), 'rs', 10], [11]])
        self.assertEqual(self.srt_rules.natur('10,11.1',
                                              delimiter=','),
                         [[10], [11.1]])
        self.assertEqual(self.srt_rules.natur('10.1,11',
                                              delimiter=',',
                                              src_file_colinds=None),
                         [[10.1], [11]])
        self.assertEqual(self.srt_rules.natur('10,11.1',
                                              delimiter=',',
                                              src_file_colinds=1),
                         [[11.1]])
        self.assertEqual(self.srt_rules.natur('10.1,11',
                                              delimiter=',',
                                              src_file_colinds=[1, 0]),
                         [[11], [10.1]])
        self.assertEqual(self.srt_rules.natur('10,11.1',
                                              delimiter='\t',
                                              src_file_colinds=[0]),
                         [[10, ',', 11.1]])


class SrtTests(unittest.TestCase):
    trf_bedgz_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz'
    trf_bed_path = PurePath(PurePath(__file__).parent,
                            PurePath(trf_bedgz_url).name[:-3]).as_posix()
    if not os.path.exists(trf_bed_path):
        os.system(f'''
                  wget -q -O - {trf_bedgz_url} |
                  gzip -d > {trf_bed_path}''')
    with open(trf_bed_path) as trf_bed_opened:
        trf_bed_content = trf_bed_opened.readlines()
    os.remove(trf_bed_path)

    def test_trf_floats(self):
        if not os.path.exists(self.trf_bed_path):
            with open(self.trf_bed_path, 'w') as trf_bed_opened:
                for trf_bed_line in self.trf_bed_content:
                    trf_bed_opened.write(trf_bed_line)
        trf_bed_floats = []
        trf_bed_lines_quan = 0
        for trf_bed_line in self.trf_bed_content:
            trf_bed_lines_quan += 1
            trf_bed_floats.append(float(trf_bed_line.split('\t')[5]))
        trf_bed_floats.sort()
        self.assertEqual(trf_bed_lines_quan,
                         432604)
        srt = Srt(self.trf_bed_path,
                  SrtRules().natur,
                  delimiter='\t',
                  src_file_colinds=5)
        srt.pre_srt(chunk_elems_quan=float('+inf'))
        nat_presrtd_trf_path = f'{self.trf_bed_path}.1'
        self.assertTrue(os.path.exists(nat_presrtd_trf_path))
        self.assertEqual(srt.presrtd_file_paths,
                         [nat_presrtd_trf_path])
        with open(nat_presrtd_trf_path) as nat_presrtd_trf_opened:
            nat_presrtd_trf_lines_quan = 0
            nat_presrtd_trf_floats = []
            for nat_presrtd_trf_line in nat_presrtd_trf_opened:
                nat_presrtd_trf_lines_quan += 1
                nat_presrtd_trf_row = nat_presrtd_trf_line.split('\t')
                self.assertEqual(len(nat_presrtd_trf_row),
                                 16)
                nat_presrtd_trf_floats.append(float(nat_presrtd_trf_row[5]))
        self.assertEqual(nat_presrtd_trf_lines_quan,
                         432604)
        self.assertEqual(trf_bed_floats,
                         nat_presrtd_trf_floats)
        srt.mrg_srt()
        nat_srtd_trf_path = f'{self.trf_bed_path}.srtd'
        self.assertTrue(os.path.exists(nat_srtd_trf_path))
        self.assertFalse(os.path.exists(nat_presrtd_trf_path))
        self.assertFalse(srt.presrtd_file_paths)
        with open(nat_srtd_trf_path) as nat_srtd_trf_opened:
            nat_srtd_trf_lines_quan = 0
            nat_srtd_trf_floats = []
            for nat_srtd_trf_line in nat_srtd_trf_opened:
                nat_srtd_trf_lines_quan += 1
                nat_srtd_trf_row = nat_srtd_trf_line.split('\t')
                self.assertEqual(len(nat_srtd_trf_row),
                                 16)
                nat_srtd_trf_floats.append(float(nat_srtd_trf_row[5]))
        self.assertEqual(nat_srtd_trf_lines_quan,
                         432604)
        self.assertEqual(trf_bed_floats,
                         nat_srtd_trf_floats)
        os.remove(nat_srtd_trf_path)
        srt.srt_rule = lambda src_file_line: float(src_file_line.split('\t')[5])
        srt.srt_rule_kwargs = {}
        srt.pre_srt(chunk_elems_quan=108151)
        float_presrtd_trf_paths = [f'{self.trf_bed_path}.1',
                                   f'{self.trf_bed_path}.2',
                                   f'{self.trf_bed_path}.3',
                                   f'{self.trf_bed_path}.4']
        for float_presrtd_trf_path in float_presrtd_trf_paths:
            self.assertTrue(os.path.exists(float_presrtd_trf_path))
            with open(float_presrtd_trf_path) as float_presrtd_trf_opened:
                float_presrtd_trf_lines_quan = 0
                for float_presrtd_trf_line in float_presrtd_trf_opened:
                    self.assertEqual(len(float_presrtd_trf_line.split('\t')),
                                     16)
                    float_presrtd_trf_lines_quan += 1
            self.assertEqual(float_presrtd_trf_lines_quan,
                             108151)
        self.assertFalse(os.path.exists(f'{self.trf_bed_path}.0'))
        self.assertFalse(os.path.exists(f'{self.trf_bed_path}.5'))
        self.assertEqual(sorted(srt.presrtd_file_paths),
                         float_presrtd_trf_paths)
        srt.mrg_srt(mrgd_file_suff='sorted',
                    del_presrtd_files=False)
        float_srtd_trf_path = f'{self.trf_bed_path}.sorted'
        self.assertTrue(os.path.exists(float_srtd_trf_path))
        with open(float_srtd_trf_path) as float_srtd_trf_opened:
            float_srtd_trf_lines_quan = 0
            float_srtd_trf_floats = []
            for float_srtd_trf_line in float_srtd_trf_opened:
                float_srtd_trf_lines_quan += 1
                float_srtd_trf_row = float_srtd_trf_line.split('\t')
                self.assertEqual(len(float_srtd_trf_row),
                                 16)
                float_srtd_trf_floats.append(float(float_srtd_trf_row[5]))
        self.assertEqual(float_srtd_trf_lines_quan,
                         432604)
        self.assertEqual(trf_bed_floats,
                         float_srtd_trf_floats)
        for float_presrtd_trf_path in float_presrtd_trf_paths:
            self.assertTrue(os.path.exists(float_presrtd_trf_path))
            os.remove(float_presrtd_trf_path)
        os.remove(float_srtd_trf_path)
        os.remove(self.trf_bed_path)


if __name__ == "__main__":
    unittest.main()
