# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import unittest
import os
import json
import re
from pathlib import PurePath
from decimal import Decimal
from random import shuffle
import pyzstd
from src.antidb.antidb import (Idx,
                               Prs)
from src.antidb.antisrt import (DelimitersMatchError,
                                SrtRules,
                                Srt)

__version__ = 'v3.1.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023-2024'}]


def get_test_content(compr_file_url: str,
                     compr_alg: str):
    file_path = PurePath(PurePath(__file__).parent,
                         PurePath(compr_file_url).name.rsplit('.',
                                                              maxsplit=1)[0]).as_posix()
    if not os.path.exists(file_path):
        os.system(f'''
                  wget -q -O - {compr_file_url} |
                  {compr_alg} -d > {file_path}''')
    with open(file_path) as file_opened:
        file_content = file_opened.read()
    os.remove(file_path)
    return file_path, file_content


def remove_antidb_test_files(idx_obj):
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


class AntidbTests(unittest.TestCase):
    mt_jsonbz2_url = 'https://ftp.ncbi.nih.gov/snp/archive/b156/JSON/refsnp-chrMT.json.bz2'
    mt_json_path, mt_json_content = get_test_content(mt_jsonbz2_url,
                                                     'bzip2')
    mt_vcfgz_url = 'https://ftp.ensembl.org/pub/release-111/variation/vcf/homo_sapiens/homo_sapiens-chrMT.vcf.gz'
    mt_vcf_path, mt_vcf_content = get_test_content(mt_vcfgz_url,
                                                   'gzip')

    def test_mt(self):
        if not os.path.exists(self.mt_json_path):
            with open(self.mt_json_path, 'w') as mt_json_opened:
                mt_json_opened.write(self.mt_json_content)
        mt_idx = Idx(self.mt_json_path,
                     'rsids',
                     lambda mt_zst_line:
                     json.loads(mt_zst_line)['refsnp_id'])
        self.assertEqual(mt_idx.srt_rule.__name__,
                         'natur')
        self.assertFalse(mt_idx.perf)
        remove_antidb_test_files(mt_idx)
        mt_idx.idx()
        self.assertIsInstance(mt_idx.perf,
                              list)
        self.assertEqual(len(mt_idx.perf),
                         5)
        self.assertEqual(mt_idx.perf[0][0:2],
                         ('crt_db_zst',
                          None))
        self.assertEqual(mt_idx.perf[1][0:2],
                         ('crt_full_idx_tmp',
                          None))
        self.assertEqual(mt_idx.perf[2][0:2],
                         ('crt_full_idx_tmp_srtd',
                          None))
        self.assertEqual(mt_idx.perf[3][0:2],
                         ('crt_full_idx',
                          None))
        self.assertEqual(mt_idx.perf[4][0:2],
                         ('crt_mem_idx',
                          None))
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
            self.assertEqual(mem_idx_lines_cnt, 11)
        mt_prs = Prs(self.mt_json_path,
                     'rsids')
        self.assertEqual(mt_prs.idx_srt_rule_name,
                         'natur')
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
        remove_antidb_test_files(mt_idx)

    def test_mt_by_mt(self):
        if not os.path.exists(self.mt_json_path):
            with open(self.mt_json_path, 'w') as mt_json_opened:
                mt_json_opened.write(self.mt_json_content)

        def parse_mt_line(mt_zst_line):
            return json.loads(mt_zst_line)['refsnp_id']

        def alpha_srt_rule(full_idx_line):
            return full_idx_line.split('\t')[0]

        mt_idx = Idx(self.mt_json_path,
                     'rsids',
                     parse_mt_line,
                     compr_frame_size=1024,
                     compr_chunk_size=1024,
                     compr_chunk_elems_quan=10,
                     unidx_lines_quan=10,
                     srt_rule=alpha_srt_rule)
        self.assertEqual(mt_idx.your_line_parser.__name__,
                         'parse_mt_line')
        self.assertEqual(mt_idx.srt_rule.__name__,
                         'alpha_srt_rule')
        self.assertEqual(mt_idx.srt_rule_kwargs,
                         {})
        self.assertEqual(mt_idx.srt_rule_settings,
                         {'col_inds': None,
                          'cols_delimiter': '\t'})
        remove_antidb_test_files(mt_idx)
        mt_idx.idx()
        with pyzstd.open(mt_idx.mem_idx_path,
                         mode='rt') as mem_idx_opened:
            idx_srt_rule_name = mem_idx_opened.readline().rstrip().split('=')[1]
            idx_srt_rule_settings = eval(mem_idx_opened.readline().rstrip().split('=')[1])
            unidx_lines_quan = int(mem_idx_opened.readline().rstrip().split('=')[1])
            self.assertEqual(idx_srt_rule_name,
                             'alpha_srt_rule')
            self.assertEqual(idx_srt_rule_settings,
                             {'col_inds': None,
                              'cols_delimiter': '\t'})
            self.assertEqual(unidx_lines_quan,
                             10)
        mt_prs = Prs(self.mt_json_path,
                     'rsids',
                     srt_rule=alpha_srt_rule)
        self.assertIsNone(mt_prs.your_line_parser)
        self.assertEqual(mt_prs.srt_rule.__name__,
                         'alpha_srt_rule')
        self.assertEqual(mt_prs.srt_rule_kwargs,
                         {})
        self.assertEqual(mt_prs.srt_rule_settings,
                         {'col_inds': None,
                          'cols_delimiter': '\t'})
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
        remove_antidb_test_files(mt_idx)

    def test_mt_clinical(self):
        if not os.path.exists(self.mt_json_path):
            with open(self.mt_json_path, 'w') as mt_json_opened:
                mt_json_opened.write(self.mt_json_content)

        def parse_mt_cln_line(mt_zst_line):
            mt_zst_obj = json.loads(mt_zst_line)
            if len(mt_zst_obj['primary_snapshot_data']['allele_annotations'][0]['clinical']) > 0:
                return mt_zst_obj['refsnp_id']
            return None

        mt_cln_idx = Idx(self.mt_json_path,
                         'rsids_cln',
                         parse_mt_cln_line)
        self.assertEqual(mt_cln_idx.your_line_parser.__name__,
                         'parse_mt_cln_line')
        self.assertEqual(mt_cln_idx.srt_rule.__name__,
                         'natur')
        remove_antidb_test_files(mt_cln_idx)
        mt_cln_idx.idx()
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
            self.assertEqual(mem_idx_lines_cnt, 4)
        mt_cln_prs = Prs(self.mt_json_path,
                         'rsids_cln')
        self.assertEqual(mt_cln_prs.srt_rule.__name__,
                         'natur')
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
        remove_antidb_test_files(mt_cln_idx)

    def test_mt_info(self):
        if not os.path.exists(self.mt_vcf_path):
            with open(self.mt_vcf_path, 'w') as mt_vcf_opened:
                mt_vcf_opened.write(self.mt_vcf_content)

        def parse_mt_vcfzst_line(mt_vcfzst_line: str,
                                 mt_vcfzst_metasymb: str = '',
                                 info_keys: None | str | list = None):
            if mt_vcfzst_line.startswith(mt_vcfzst_metasymb):
                return None
            info = mt_vcfzst_line.rstrip().split('\t')[-1]
            if not info_keys:
                return info
            elif type(info_keys) is str:
                info_keys = [info_keys]
            info_row = info.split(';')
            idx_vals = []
            for info_cell in info_row:
                for info_key in info_keys:
                    if info_key == info_cell:
                        idx_vals.append(info_cell)
                    elif info_key == info_cell.split('=')[0]:
                        idx_vals.append(info_cell.split('=')[1])
            if idx_vals:
                return idx_vals
            else:
                return None

        def assign_idxval_type(idx_line,
                               datatype):
            idx_val = datatype(idx_line.split(',')[0])
            return idx_val

        mt_info_idx = Idx(self.mt_vcf_path,
                          'info',
                          your_line_parser=parse_mt_vcfzst_line,
                          your_line_parser_kwargs={'mt_vcfzst_metasymb': '#',
                                                   'info_keys': ['CLIN_pathogenic',
                                                                 'TSA']},
                          srt_rule=assign_idxval_type,
                          srt_rule_kwargs={'datatype': str})
        self.assertEqual(mt_info_idx.your_line_parser.__name__,
                         'parse_mt_vcfzst_line')
        self.assertEqual(mt_info_idx.your_line_parser_kwargs,
                         {'mt_vcfzst_metasymb': '#',
                          'info_keys': ['CLIN_pathogenic',
                                        'TSA']})
        self.assertEqual(mt_info_idx.srt_rule.__name__,
                         'assign_idxval_type')
        self.assertEqual(mt_info_idx.srt_rule_kwargs,
                         {'datatype': str})
        remove_antidb_test_files(mt_info_idx)
        mt_info_idx.idx()
        self.assertFalse(os.path.exists(mt_info_idx.db_file_path))
        self.assertTrue(os.path.exists(mt_info_idx.full_idx_path))
        self.assertFalse(os.path.exists(mt_info_idx.full_idx_tmp_path))
        self.assertFalse(os.path.exists(mt_info_idx.full_idx_tmp_srtd_path))
        self.assertTrue(os.path.exists(mt_info_idx.mem_idx_path))
        mt_vcf_content = self.mt_vcf_content.split('\n')
        content_idx_vals = []
        for mt_vcf_line in mt_vcf_content:
            if mt_vcf_line.startswith('#'):
                continue
            if 'CLIN_pathogenic' in mt_vcf_line:
                content_idx_vals.append('CLIN_pathogenic')
            tsa = re.search(r'(?<=TSA=)\w+',
                            mt_vcf_line)
            if tsa:
                content_idx_vals.append(tsa.group())
        content_idx_vals.sort()
        full_idx_vals = []
        with pyzstd.open(mt_info_idx.full_idx_path,
                         mode='rt') as full_idx_opened:
            for full_idx_line in full_idx_opened:
                full_idx_vals.append(full_idx_line.split(',')[0])
        for ind in range(len(content_idx_vals)):
            self.assertEqual(content_idx_vals[ind],
                             full_idx_vals[ind])
        mt_info_prs = Prs(self.mt_vcf_path,
                          'info',
                          srt_rule=assign_idxval_type,
                          srt_rule_kwargs={'datatype': str})
        self.assertEqual(str(mt_info_prs.srt_rule_settings),
                         str({'cols_delimiter': '\t',
                              'col_inds': None,
                              'datatype': str}))
        content_sub_ins_lines = sorted([content_line for content_line in mt_vcf_content
                                        if 'substitution' in content_line
                                        or 'insertion' in content_line],
                                       key=lambda content_line:
                                       content_line.rstrip().split('\t')[7].split(';')[1],
                                       reverse=True)
        self.assertEqual(len(content_sub_ins_lines),
                         36)
        prs_sub_ins_lines = [prs_sub_ins_line.rstrip()
                             for prs_sub_ins_line in mt_info_prs.prs(['substitution',
                                                                      'insertion'])]
        self.assertEqual(len(prs_sub_ins_lines),
                         36)
        for ind in range(36):
            self.assertEqual(content_sub_ins_lines[ind],
                             prs_sub_ins_lines[ind])
        self.assertEqual(len(list(mt_info_prs.prs('CLIN_pathogenic'))),
                         262)
        self.assertFalse(list(mt_info_prs.prs('CLIN')))
        remove_antidb_test_files(mt_info_idx)


class SrtRulesTests(unittest.TestCase):
    srt_rules = SrtRules()

    def test_get_cols_srt_rule(self):
        self.assertEqual(self.srt_rules.cols_delimiter,
                         '\t')
        self.assertIsNone(self.srt_rules.col_inds)
        self.assertEqual(self.srt_rules.get_cols('abc\tdef'),
                         ['abc', 'def'])
        self.assertEqual(self.srt_rules.get_cols('123\t456'),
                         ['123', '456'])
        self.srt_rules.cols_delimiter = ','
        self.assertEqual(self.srt_rules.get_cols('chr10\t1\t100\tA\tC,CC'),
                         ['chr10\t1\t100\tA\tC', 'CC'])
        self.srt_rules.cols_delimiter = None
        self.assertEqual(self.srt_rules.get_cols('chr20\t2\t200\tA\tC,CC'),
                         ['chr20\t2\t200\tA\tC,CC'])
        self.srt_rules.cols_delimiter = '\t'
        self.assertEqual(self.srt_rules.get_cols('chr30\t3\t300\tA\tC,CC'),
                         ['chr30', '3', '300', 'A', 'C,CC'])
        self.srt_rules.col_inds = [0, -1]
        self.assertEqual(self.srt_rules.get_cols('chr40\t4\t400\tA\tC,CC'),
                         ['chr40', 'C,CC'])
        self.srt_rules.col_inds = None

    def test_natur_srt_rule(self):
        self.assertEqual(self.srt_rules.cols_delimiter,
                         '\t')
        self.assertIsNone(self.srt_rules.col_inds)
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
        self.srt_rules.cols_delimiter = ','
        self.assertEqual(self.srt_rules.natur('10,11.1'),
                         [[10], [11.1]])
        self.srt_rules.col_inds = None
        self.assertEqual(self.srt_rules.natur('10.1,11'),
                         [[10.1], [11]])
        self.srt_rules.col_inds = 1
        self.assertEqual(self.srt_rules.natur('10,11.1'),
                         [[11.1]])
        self.srt_rules.col_inds = [1, 0]
        self.assertEqual(self.srt_rules.natur('10.1,11'),
                         [[11], [10.1]])
        self.srt_rules.cols_delimiter = '\t'
        self.srt_rules.col_inds = [0]
        self.assertEqual(self.srt_rules.natur('10,11.1'),
                         [[10, ',', 11.1]])
        self.srt_rules.col_inds = None
        self.assertEqual(self.srt_rules.natur('10.1,11',
                                              dec_delimiter=','),
                         [[10, '.', 1.11]])
        self.assertEqual(self.srt_rules.natur('+'),
                         [[float('+inf'), '+']])
        self.assertEqual(self.srt_rules.natur('-'),
                         [[float('+inf'), '-']])
        self.assertEqual(self.srt_rules.natur('1+1'),
                         [[1, '+', 1]])
        self.assertEqual(self.srt_rules.natur('-1-1'),
                         [[-1, -1]])
        self.assertEqual(self.srt_rules.natur('+001-001'),
                         [[float('+inf'), '+', 1, -1]])
        self.assertEqual(self.srt_rules.natur('-1.23e'),
                         [[-1.23, 'e']])
        self.assertEqual(self.srt_rules.natur('123E-3'),
                         [[0.123]])
        self.assertEqual(self.srt_rules.natur('123e-02'),
                         [[1.23]])
        self.assertEqual(self.srt_rules.natur('pref1.23e2suff'),
                         [[float('+inf'), 'pref', 123.0, 'suff']])
        self.assertEqual(self.srt_rules.natur('e1.23e-1e'),
                         [[float('+inf'), 'e', 0.123, 'e']])
        self.assertEqual(self.srt_rules.natur('-E1.23E+01-E'),
                         [[float('+inf'), '-E', 12.3, '-E']])
        self.assertEqual(self.srt_rules.natur('-e-1.23e-1-e'),
                         [[float('+inf'), '-e', -0.123, '-e']])
        self.assertEqual(self.srt_rules.natur('+1230E-1+1'),
                         [[float('+inf'), '+', 123.0, '+', 1]])
        self.assertEqual(self.srt_rules.natur('0.123ee+2'),
                         [[0.123, 'ee+', 2]])
        self.assertEqual(self.srt_rules.natur('-12.3+EE-02'),
                         [[-12.3, '+EE', -2]])
        self.assertEqual(self.srt_rules.natur('-0e10'),
                         [[0]])
        self.assertEqual(self.srt_rules.natur('qwerty\t0.1E2'),
                         [[float('+inf'), 'qwerty'], [10.0]])
        self.assertEqual(self.srt_rules.natur('0,1e+2\tqwerty',
                                              dec_delimiter=','),
                         [[10.0], [float('+inf'), 'qwerty']])
        self.srt_rules.cols_delimiter = ', '
        self.assertEqual(self.srt_rules.natur('1,1, -2,2, str',
                                              dec_delimiter=',',
                                              nums_first=False),
                         [[1.1], [-2.2], [float('-inf'), 'str']])
        self.srt_rules.cols_delimiter = ','
        self.assertRaises(DelimitersMatchError,
                          self.srt_rules.natur,
                          '111,111',
                          dec_delimiter=',')
        self.srt_rules.cols_delimiter = '\t'
        self.srt_rules.col_inds = None

    def test_letts_nums_srt_rule(self):
        self.assertEqual(self.srt_rules.cols_delimiter,
                         '\t')
        self.assertIsNone(self.srt_rules.col_inds)
        self.assertEqual(self.srt_rules.letts_nums('rs1'),
                         [['rs', 1]])
        self.assertEqual(self.srt_rules.letts_nums('rs010'),
                         [['rs', 10]])
        self.srt_rules.col_inds = 0
        self.assertEqual(self.srt_rules.letts_nums('ENSG000\trs000'),
                         [['ENSG', 0]])
        self.srt_rules.col_inds = [1, 0]
        self.assertEqual(self.srt_rules.letts_nums('ENSG000\trs000'),
                         [['rs', 0], ['ENSG', 0]])
        self.srt_rules.col_inds = None
        self.assertRaises(AttributeError,
                          self.srt_rules.letts_nums,
                          'rs')
        self.assertRaises(AttributeError,
                          self.srt_rules.letts_nums,
                          '1dvatri')
        self.assertRaises(AttributeError,
                          self.srt_rules.letts_nums,
                          '123')
        self.assertRaises(AttributeError,
                          self.srt_rules.letts_nums,
                          'id1.1')
        self.assertRaises(AttributeError,
                          self.srt_rules.letts_nums,
                          'id-1')


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
            trf_bed_row = trf_bed_line.split('\t')
            self.assertEqual(len(trf_bed_row),
                             16)
            trf_bed_floats.append(float(trf_bed_row[5]))
        trf_bed_floats.sort()
        self.assertEqual(trf_bed_lines_quan,
                         432604)
        srt = Srt(self.trf_bed_path,
                  SrtRules().natur,
                  cols_delimiter='\t',
                  col_inds=5,
                  nums_first=True)
        self.assertEqual(srt.srt_rule.__name__,
                         'natur')
        self.assertEqual(srt.srt_rule_kwargs,
                         {'nums_first': True})
        self.assertEqual(srt.srt_rule_settings,
                         {'cols_delimiter': '\t',
                          'col_inds': 5,
                          'nums_first': True})
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

    def test_trf_2_cols(self):
        if not os.path.exists(self.trf_bed_path):
            with open(self.trf_bed_path, 'w') as trf_bed_opened:
                for trf_bed_line in self.trf_bed_content:
                    trf_bed_opened.write(trf_bed_line)
        str_float_gnusrtd_path = f'{self.trf_bed_path}.gnusrtd'
        os.system(f"""LC_ALL=C sort -s -t '{chr(9)}' -k16,16 -k15,15n \
                  {self.trf_bed_path} > {str_float_gnusrtd_path}""")
        self.assertTrue(os.path.exists(str_float_gnusrtd_path))
        with open(str_float_gnusrtd_path) as str_float_gnusrtd_opened:
            str_float_gnusrtd_lines_quan = 0
            gnusrtd_strs_floats = []
            for str_float_gnusrtd_line in str_float_gnusrtd_opened:
                str_float_gnusrtd_lines_quan += 1
                str_float_gnusrtd_row = str_float_gnusrtd_line.split('\t')
                self.assertEqual(len(str_float_gnusrtd_row),
                                 16)
                gnusrtd_strs_floats.append([float(str_float_gnusrtd_row[14]),
                                            str_float_gnusrtd_row[15]])
        self.assertEqual(str_float_gnusrtd_lines_quan,
                         432604)
        srt = Srt(self.trf_bed_path,
                  srt_rule=None,
                  cols_delimiter='\t',
                  col_inds=[15, 14])
        self.assertEqual(srt.srt_rule.__name__,
                         'natur')
        self.assertEqual(srt.srt_rule_kwargs,
                         {})
        self.assertEqual(srt.srt_rule_settings,
                         {'cols_delimiter': '\t',
                          'col_inds': [15, 14]})
        srt.pre_srt(chunk_elems_quan=108150)
        self.assertEqual(len(srt.presrtd_file_paths),
                         5)
        for presrtd_file_path in srt.presrtd_file_paths:
            self.assertTrue(os.path.exists(presrtd_file_path))
        str_float_mrgd_path = srt.mrg_srt(mrgd_file_suff='mrgd')
        self.assertFalse(srt.presrtd_file_paths)
        self.assertTrue(str_float_mrgd_path)
        with open(str_float_mrgd_path) as str_float_mrgd_opened:
            str_float_mrgd_lines_quan = 0
            mrgd_strs_floats = []
            for str_float_mrgd_line in str_float_mrgd_opened:
                str_float_mrgd_lines_quan += 1
                str_float_mrgd_row = str_float_mrgd_line.split('\t')
                self.assertEqual(len(str_float_mrgd_row),
                                 16)
                mrgd_strs_floats.append([float(str_float_mrgd_row[14]),
                                         str_float_mrgd_row[15]])
        self.assertEqual(str_float_mrgd_lines_quan,
                         432604)
        self.assertEqual(mrgd_strs_floats,
                         gnusrtd_strs_floats)
        os.remove(self.trf_bed_path)
        os.remove(str_float_gnusrtd_path)
        os.remove(str_float_mrgd_path)


if __name__ == "__main__":
    unittest.main()
