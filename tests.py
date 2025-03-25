# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import unittest
from antidb.srt import *
from antidb.idx import *
from antidb.prs import *

if __name__ == 'main':
    __version__ = 'v5.0.0'
    __authors__ = [{'name': 'Platon Bykadorov',
                    'email': 'platon.work@gmail.com',
                    'years': '2023-2025'}]


def del_files(*file_paths: str) -> None:
    for file_path in file_paths:
        if os.path.exists(file_path):
            os.remove(file_path)


class BedTests(unittest.TestCase):
    src_bed = ['1\t116545156\t116545157\trs12044852\n',
               '1\t241782991\t241782992\trs952084\n',
               '1\t154527612\t154527613\trs4131514\n',
               '1\t201015351\t201015352\trs12122721\n',
               '1\t92515681\t92515682\trs17371561\n',
               '1\t92543755\t92543756\trs11804321\n',
               '1\t92580419\t92580420\trs17380378\n',
               '1\t24977084\t24977085\trs10903122\n',
               '1\t92516860\t92516861\trs11581176\n',
               '1\t92543755\t92543756\trs11804321\n',
               '1\t92516860\t92516861\trs11581176\n',
               '1\t86877126\t86877127\trs581405\n',
               '1\t237115473\t237115474\trs10925318\n',
               '1\t86876786\t86876787\trs479341\n',
               '1\t66265029\t66265030\trs1321172\n',
               '1\t92543755\t92543756\trs11804321\n']
    src_file_path = os.path.join(os.getcwd(),
                                 'bed.bed')
    db_zst_path = os.path.join(os.getcwd(),
                               'bed.bed.zst')

    def test_common(self):
        adb_path = os.path.join(os.getcwd(),
                                'bed.bed.rsids.adb')
        with open(self.src_file_path, 'w') as src_file_opened:
            for src_bed_line in self.src_bed:
                src_file_opened.write(src_bed_line)
        self.assertTrue(os.path.isfile(self.src_file_path))
        del_files(self.db_zst_path,
                  adb_path)
        idx_obj = Idx(db_file_path=self.src_file_path,
                      adb_name_prefix='rsids',
                      db_line_prs=(lambda line:
                                   line.split('\t')[-1]),
                      adb_srt_rule=SrtRules.natur,
                      presrt_chunk_len=8,
                      lstarts_idx_div=2,
                      lstarts_idx_len=5)
        idx_obj.idx()
        self.assertEqual(len(idx_obj.presrtd_idxs_opened), 2)
        self.assertTrue(os.path.isfile(self.src_file_path))
        self.assertTrue(os.path.isfile(self.db_zst_path))
        self.assertTrue(os.path.isfile(adb_path))
        with ZipFile(adb_path) as adb_opened_r:
            adb_content = adb_opened_r.namelist()
            paths_idx_paths = []
            lstarts_idx_paths = []
            for file_path in adb_content:
                if os.path.basename(file_path) == 'paths':
                    paths_idx_paths.append(file_path)
                elif os.path.basename(file_path) == 'lstarts':
                    lstarts_idx_paths.append(file_path)
            self.assertEqual(sorted(paths_idx_paths),
                             ['1/paths',
                              '2/paths',
                              'paths'])
            self.assertEqual(sorted(lstarts_idx_paths),
                             ['1/1/lstarts',
                              '1/2/lstarts',
                              '2/1/lstarts',
                              '2/2/lstarts'])
            self.assertEqual(len(adb_content),
                             13)
            root_paths_idx_path = 'paths'
            with adb_opened_r.open(root_paths_idx_path) as root_paths_idx_opened:
                root_paths_idx_obj = load(root_paths_idx_opened)
                self.assertEqual(root_paths_idx_obj[0][0],
                                 [[inf, 'rs', 479341]])
                self.assertEqual(root_paths_idx_obj[0][1],
                                 [[inf, 'rs', 11581176]])
                self.assertEqual(root_paths_idx_obj[1][0],
                                 '1/paths')
                self.assertEqual(root_paths_idx_obj[1][1],
                                 '2/paths')
                self.assertEqual(len(root_paths_idx_obj),
                                 2)
                self.assertEqual(len(root_paths_idx_obj[0]),
                                 2)
                self.assertEqual(len(root_paths_idx_obj[1]),
                                 2)
            fir_paths_idx_path = '1/paths'
            with adb_opened_r.open(fir_paths_idx_path) as fir_paths_idx_opened:
                fir_paths_idx_obj = load(fir_paths_idx_opened)
                self.assertEqual(fir_paths_idx_obj[0][0],
                                 [[inf, 'rs', 479341]])
                self.assertEqual(fir_paths_idx_obj[0][1],
                                 [[inf, 'rs', 4131514]])
                self.assertEqual(fir_paths_idx_obj[1][0],
                                 '1/1/lstarts')
                self.assertEqual(fir_paths_idx_obj[1][1],
                                 '1/2/lstarts')
                self.assertEqual(len(fir_paths_idx_obj),
                                 2)
                self.assertEqual(len(fir_paths_idx_obj[0]),
                                 2)
                self.assertEqual(len(fir_paths_idx_obj[1]),
                                 2)
            sec_paths_idx_path = '2/paths'
            with adb_opened_r.open(sec_paths_idx_path) as sec_paths_idx_opened:
                sec_paths_idx_obj = load(sec_paths_idx_opened)
                self.assertEqual(sec_paths_idx_obj[0][0],
                                 [[inf, 'rs', 11581176]])
                self.assertEqual(sec_paths_idx_obj[0][1],
                                 [[inf, 'rs', 12044852]])
                self.assertEqual(sec_paths_idx_obj[1][0],
                                 '2/1/lstarts')
                self.assertEqual(sec_paths_idx_obj[1][1],
                                 '2/2/lstarts')
                self.assertEqual(len(sec_paths_idx_obj),
                                 2)
                self.assertEqual(len(sec_paths_idx_obj[0]),
                                 2)
                self.assertEqual(len(sec_paths_idx_obj[1]),
                                 2)
            self.assertEqual(sorted(self.src_bed,
                                    key=lambda line:
                                    int(line.rstrip().split('\t')[-1][2:])),
                             ['1\t86876786\t86876787\trs479341\n',
                              '1\t86877126\t86877127\trs581405\n',
                              '1\t241782991\t241782992\trs952084\n',
                              '1\t66265029\t66265030\trs1321172\n',
                              '1\t154527612\t154527613\trs4131514\n',
                              '1\t24977084\t24977085\trs10903122\n',
                              '1\t237115473\t237115474\trs10925318\n',
                              '1\t92516860\t92516861\trs11581176\n',
                              '1\t92516860\t92516861\trs11581176\n',
                              '1\t92543755\t92543756\trs11804321\n',
                              '1\t92543755\t92543756\trs11804321\n',
                              '1\t92543755\t92543756\trs11804321\n',
                              '1\t116545156\t116545157\trs12044852\n',
                              '1\t201015351\t201015352\trs12122721\n',
                              '1\t92515681\t92515682\trs17371561\n',
                              '1\t92580419\t92580420\trs17380378\n'])
            fir_lstarts_idx_path = '1/1/lstarts'
            with ZstdFile(adb_opened_r.open(fir_lstarts_idx_path)) as fir_lstarts_idx_opened:
                fir_lstarts_idx_obj = load(fir_lstarts_idx_opened)
                self.assertEqual(fir_lstarts_idx_obj[0][0],
                                 [[inf, 'rs', 479341]])
                self.assertEqual(fir_lstarts_idx_obj[0][1],
                                 [[inf, 'rs', 581405]])
                self.assertEqual(fir_lstarts_idx_obj[0][2],
                                 [[inf, 'rs', 952084]])
                self.assertEqual(fir_lstarts_idx_obj[0][3],
                                 [[inf, 'rs', 1321172]])
                self.assertEqual(len(fir_lstarts_idx_obj),
                                 2)
                self.assertEqual(len(fir_lstarts_idx_obj[0]),
                                 4)
                self.assertEqual(len(fir_lstarts_idx_obj[1]),
                                 4)
            sec_lstarts_idx_path = '1/2/lstarts'
            with ZstdFile(adb_opened_r.open(sec_lstarts_idx_path)) as sec_lstarts_idx_opened:
                sec_lstarts_idx_obj = load(sec_lstarts_idx_opened)
                self.assertEqual(sec_lstarts_idx_obj[0][0],
                                 [[inf, 'rs', 4131514]])
                self.assertEqual(sec_lstarts_idx_obj[0][1],
                                 [[inf, 'rs', 10903122]])
                self.assertEqual(sec_lstarts_idx_obj[0][2],
                                 [[inf, 'rs', 10925318]])
                self.assertEqual(sec_lstarts_idx_obj[0][3],
                                 [[inf, 'rs', 11581176]])
                self.assertEqual(len(sec_lstarts_idx_obj),
                                 2)
                self.assertEqual(len(sec_lstarts_idx_obj[0]),
                                 4)
                self.assertEqual(len(sec_lstarts_idx_obj[1]),
                                 4)
            thi_lstarts_idx_path = '2/1/lstarts'
            with ZstdFile(adb_opened_r.open(thi_lstarts_idx_path)) as thi_lstarts_idx_opened:
                thi_lstarts_idx_obj = load(thi_lstarts_idx_opened)
                self.assertEqual(thi_lstarts_idx_obj[0][0],
                                 [[inf, 'rs', 11581176]])
                self.assertEqual(thi_lstarts_idx_obj[0][1],
                                 [[inf, 'rs', 11804321]])
                self.assertEqual(thi_lstarts_idx_obj[0][2],
                                 [[inf, 'rs', 11804321]])
                self.assertEqual(thi_lstarts_idx_obj[0][3],
                                 [[inf, 'rs', 11804321]])
                self.assertEqual(len(thi_lstarts_idx_obj),
                                 2)
                self.assertEqual(len(thi_lstarts_idx_obj[0]),
                                 4)
                self.assertEqual(len(thi_lstarts_idx_obj[1]),
                                 4)
            fou_lstarts_idx_path = '2/2/lstarts'
            with ZstdFile(adb_opened_r.open(fou_lstarts_idx_path)) as fou_lstarts_idx_opened:
                fou_lstarts_idx_obj = load(fou_lstarts_idx_opened)
                self.assertEqual(fou_lstarts_idx_obj[0][0],
                                 [[inf, 'rs', 12044852]])
                self.assertEqual(fou_lstarts_idx_obj[0][1],
                                 [[inf, 'rs', 12122721]])
                self.assertEqual(fou_lstarts_idx_obj[0][2],
                                 [[inf, 'rs', 17371561]])
                self.assertEqual(fou_lstarts_idx_obj[0][3],
                                 [[inf, 'rs', 17380378]])
                self.assertEqual(len(fou_lstarts_idx_obj),
                                 2)
                self.assertEqual(len(fou_lstarts_idx_obj[0]),
                                 4)
                self.assertEqual(len(fou_lstarts_idx_obj[1]),
                                 4)
        prs_obj = Prs(db_file_path=self.src_file_path,
                      adb_name_prefix='rsids',
                      adb_srt_rule=SrtRules.natur)
        self.assertEqual(prs_obj.read_lstarts_idx('1/1/lstarts'),
                         fir_lstarts_idx_obj)
        self.assertEqual(prs_obj.read_lstarts_idx('1/2/lstarts'),
                         sec_lstarts_idx_obj)
        self.assertEqual(prs_obj.read_lstarts_idx('2/1/lstarts'),
                         thi_lstarts_idx_obj)
        self.assertEqual(prs_obj.read_lstarts_idx('2/2/lstarts'),
                         fou_lstarts_idx_obj)
        self.assertEqual(list(prs_obj.eq('rs12044852')),
                         ['1\t116545156\t116545157\trs12044852\n'])
        self.assertEqual(list(prs_obj.eq('rs952084')),
                         ['1\t241782991\t241782992\trs952084\n'])
        self.assertEqual(list(prs_obj.eq('rs4131514')),
                         ['1\t154527612\t154527613\trs4131514\n'])
        self.assertEqual(list(prs_obj.eq('rs12122721')),
                         ['1\t201015351\t201015352\trs12122721\n'])
        self.assertEqual(list(prs_obj.eq('rs17371561')),
                         ['1\t92515681\t92515682\trs17371561\n'])
        self.assertEqual(list(prs_obj.eq('rs11804321')),
                         ['1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n'])
        self.assertEqual(list(prs_obj.eq('rs17380378')),
                         ['1\t92580419\t92580420\trs17380378\n'])
        self.assertEqual(list(prs_obj.eq('rs10903122')),
                         ['1\t24977084\t24977085\trs10903122\n'])
        self.assertEqual(list(prs_obj.eq('rs11581176')),
                         ['1\t92516860\t92516861\trs11581176\n',
                          '1\t92516860\t92516861\trs11581176\n'])
        self.assertEqual(list(prs_obj.eq('rs581405')),
                         ['1\t86877126\t86877127\trs581405\n'])
        self.assertEqual(list(prs_obj.eq('rs10925318')),
                         ['1\t237115473\t237115474\trs10925318\n'])
        self.assertEqual(list(prs_obj.eq('rs479341')),
                         ['1\t86876786\t86876787\trs479341\n'])
        self.assertEqual(list(prs_obj.eq('rs1321172')),
                         ['1\t66265029\t66265030\trs1321172\n'])
        self.assertEqual(list(prs_obj.eq('rs00000')),
                         [])
        self.assertEqual(list(prs_obj.eq('rs11804321',
                                         'rs11581176',
                                         'hz12345')),
                         ['1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t92516860\t92516861\trs11581176\n'])
        self.assertEqual(list(prs_obj.rng('rs479341',
                                          'rs952084')),
                         ['1\t86876786\t86876787\trs479341\n',
                          '1\t86877126\t86877127\trs581405\n',
                          '1\t241782991\t241782992\trs952084\n'])
        self.assertEqual(list(prs_obj.rng('rs4131514',
                                          'rs10903122')),
                         ['1\t154527612\t154527613\trs4131514\n',
                          '1\t24977084\t24977085\trs10903122\n'])
        self.assertEqual(list(prs_obj.rng('rs11804321',
                                          'rs11804321')),
                         ['1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n'])
        self.assertEqual(list(prs_obj.rng('rs12122721',
                                          'rs17380378')),
                         ['1\t201015351\t201015352\trs12122721\n',
                          '1\t92515681\t92515682\trs17371561\n',
                          '1\t92580419\t92580420\trs17380378\n'])
        self.assertEqual(list(prs_obj.rng('rs00000',
                                          'rs480000')),
                         ['1\t86876786\t86876787\trs479341\n'])
        self.assertEqual(list(prs_obj.rng('rs4000000',
                                          'rs11900000')),
                         ['1\t154527612\t154527613\trs4131514\n',
                          '1\t24977084\t24977085\trs10903122\n',
                          '1\t237115473\t237115474\trs10925318\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n'])
        self.assertEqual(list(prs_obj.rng('rs17380000',
                                          'rs99999999')),
                         ['1\t92580419\t92580420\trs17380378\n'])
        self.assertRaises(QueryStartGtEndError,
                          lambda query_start, query_end:
                          list(prs_obj.rng(query_start,
                                           query_end)),
                          'rs11900000',
                          'rs4000000')
        del_files(self.src_file_path,
                  self.db_zst_path,
                  adb_path)

    def test_coords(self):
        adb_path = os.path.join(os.getcwd(),
                                'bed.bed.coords.adb')
        with open(self.src_file_path, 'w') as src_file_opened:
            for src_bed_line in self.src_bed:
                src_file_opened.write(src_bed_line)
        del_files(self.db_zst_path,
                  adb_path)

        def get_coords(src_bed_line: str):
            src_bed_row = src_bed_line.rstrip().split('\t')
            coords = [f'chr{src_bed_row[0]}',
                      int(src_bed_row[1]),
                      int(src_bed_row[2])]
            return coords

        idx_obj = Idx(db_file_path=self.src_file_path,
                      adb_name_prefix='coords',
                      db_line_prs=get_coords,
                      adb_srt_rule=lambda val: val,
                      presrt_chunk_len=12,
                      lstarts_idx_div=2,
                      lstarts_idx_len=8)
        idx_obj.idx()
        self.assertEqual(sorted(self.src_bed,
                                key=get_coords),
                         ['1\t24977084\t24977085\trs10903122\n',
                          '1\t66265029\t66265030\trs1321172\n',
                          '1\t86876786\t86876787\trs479341\n',
                          '1\t86877126\t86877127\trs581405\n',
                          '1\t92515681\t92515682\trs17371561\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92580419\t92580420\trs17380378\n',
                          '1\t116545156\t116545157\trs12044852\n',
                          '1\t154527612\t154527613\trs4131514\n',
                          '1\t201015351\t201015352\trs12122721\n',
                          '1\t237115473\t237115474\trs10925318\n',
                          '1\t241782991\t241782992\trs952084\n'])
        with ZipFile(adb_path) as adb_opened_r:
            fir_lstarts_idx_path = '1/1/lstarts'
            with ZstdFile(adb_opened_r.open(fir_lstarts_idx_path)) as fir_lstarts_idx_opened:
                fir_lstarts_idx_obj = load(fir_lstarts_idx_opened)
                self.assertEqual(fir_lstarts_idx_obj[0][0],
                                 ['chr1', 24977084, 24977085])
                self.assertEqual(fir_lstarts_idx_obj[0][1],
                                 ['chr1', 66265029, 66265030])
                self.assertEqual(fir_lstarts_idx_obj[0][2],
                                 ['chr1', 86876786, 86876787])
                self.assertEqual(fir_lstarts_idx_obj[0][3],
                                 ['chr1', 86877126, 86877127])
                self.assertEqual(fir_lstarts_idx_obj[0][4],
                                 ['chr1', 92515681, 92515682])
                self.assertEqual(fir_lstarts_idx_obj[0][5],
                                 ['chr1', 92516860, 92516861])
                self.assertEqual(len(fir_lstarts_idx_obj[0]), 6)
            sec_lstarts_idx_path = '1/2/lstarts'
            with ZstdFile(adb_opened_r.open(sec_lstarts_idx_path)) as sec_lstarts_idx_opened:
                sec_lstarts_idx_obj = load(sec_lstarts_idx_opened)
                self.assertEqual(sec_lstarts_idx_obj[0][0],
                                 ['chr1', 92516860, 92516861])
                self.assertEqual(sec_lstarts_idx_obj[0][1],
                                 ['chr1', 92543755, 92543756])
                self.assertEqual(sec_lstarts_idx_obj[0][2],
                                 ['chr1', 92543755, 92543756])
                self.assertEqual(sec_lstarts_idx_obj[0][3],
                                 ['chr1', 92543755, 92543756])
                self.assertEqual(sec_lstarts_idx_obj[0][4],
                                 ['chr1', 92580419, 92580420])
                self.assertEqual(sec_lstarts_idx_obj[0][5],
                                 ['chr1', 116545156, 116545157])
                self.assertEqual(len(sec_lstarts_idx_obj[0]),
                                 6)
            thi_lstarts_idx_path = '2/lstarts'
            with ZstdFile(adb_opened_r.open(thi_lstarts_idx_path)) as thi_lstarts_idx_opened:
                thi_lstarts_idx_obj = load(thi_lstarts_idx_opened)
                self.assertEqual(thi_lstarts_idx_obj[0][0],
                                 ['chr1', 154527612, 154527613])
                self.assertEqual(thi_lstarts_idx_obj[0][1],
                                 ['chr1', 201015351, 201015352])
                self.assertEqual(thi_lstarts_idx_obj[0][2],
                                 ['chr1', 237115473, 237115474])
                self.assertEqual(thi_lstarts_idx_obj[0][3],
                                 ['chr1', 241782991, 241782992])
                self.assertEqual(len(thi_lstarts_idx_obj[0]),
                                 4)
        prs_obj = Prs(db_file_path=self.src_file_path,
                      adb_name_prefix='coords',
                      adb_srt_rule=lambda val: val)
        self.assertEqual(list(prs_obj.eq(['chr1', 24977084, 24977085])),
                         ['1\t24977084\t24977085\trs10903122\n'])
        self.assertEqual(list(prs_obj.eq(['chr1', 92543755, 92543756])),
                         ['1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n'])
        self.assertEqual(list(prs_obj.eq(['chr1', 237115473, 237115474],
                                         ['chr1', 241782991, 241782992])),
                         ['1\t237115473\t237115474\trs10925318\n',
                          '1\t241782991\t241782992\trs952084\n'])
        self.assertEqual(list(prs_obj.rng(['chr1', 92516000, 92516001],
                                          ['chr1', 92543800, 92543801])),
                         ['1\t92516860\t92516861\trs11581176\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92543755\t92543756\trs11804321\n'])
        del_files(self.src_file_path,
                  self.db_zst_path,
                  adb_path)

    def test_homogen_idx(self):
        adb_path = os.path.join(os.getcwd(),
                                'bed.bed.chroms.adb')
        with open(self.src_file_path, 'w') as src_file_opened:
            for src_bed_line in self.src_bed:
                src_file_opened.write(src_bed_line)
        del_files(self.db_zst_path,
                  adb_path)
        idx_obj = Idx(db_file_path=self.src_file_path,
                      adb_name_prefix='chroms',
                      db_line_prs=(lambda src_bed_line:
                                   src_bed_line.split('\t')[0]),
                      adb_srt_rule=lambda val: val,
                      presrt_chunk_len=100,
                      lstarts_idx_div=1,
                      lstarts_idx_len=100)
        idx_obj.idx()
        with ZipFile(adb_path) as adb_opened_r:
            lstarts_idx_path = '1/lstarts'
            with ZstdFile(adb_opened_r.open(lstarts_idx_path)) as lstarts_idx_opened:
                lstarts_idx_obj = load(lstarts_idx_opened)
                for lstarts_idx_elem in lstarts_idx_obj[0]:
                    self.assertEqual(lstarts_idx_elem,
                                     '1')
                self.assertEqual(len(lstarts_idx_obj[0]),
                                 16)
        prs_obj = Prs(db_file_path=self.src_file_path,
                      adb_name_prefix='chroms',
                      adb_srt_rule=lambda val: val)
        self.assertEqual(sorted(prs_obj.
                                adb_opened_r.
                                namelist()),
                         ['1/',
                          '1/lstarts',
                          'paths'])
        self.assertEqual(list(prs_obj.eq('1')),
                         self.src_bed)
        self.assertEqual(list(prs_obj.rng('1',
                                          '1')),
                         self.src_bed)
        self.assertEqual(list(prs_obj.rng('-1',
                                          '2')),
                         self.src_bed)
        del_files(self.src_file_path,
                  self.db_zst_path,
                  adb_path)

    def test_rsids_by_len(self):
        adb_path = os.path.join(os.getcwd(),
                                'bed.bed.rsids_len.adb')
        with open(self.src_file_path, 'w') as src_file_opened:
            for src_bed_line in self.src_bed:
                src_file_opened.write(src_bed_line)
        del_files(self.db_zst_path,
                  adb_path)
        idx_obj = Idx(db_file_path=self.src_file_path,
                      adb_name_prefix='rsids_len',
                      db_line_prs=(lambda src_bed_line:
                                   len(src_bed_line.split('\t')[3])),
                      adb_srt_rule=lambda val: val,
                      presrt_chunk_len=8,
                      lstarts_idx_div=3,
                      lstarts_idx_len=3)
        idx_obj.idx()
        self.assertEqual(sorted(self.src_bed,
                                key=lambda line:
                                len(line.rstrip().split('\t')[-1])),
                         ['1\t241782991\t241782992\trs952084\n',
                          '1\t86877126\t86877127\trs581405\n',
                          '1\t86876786\t86876787\trs479341\n',
                          '1\t154527612\t154527613\trs4131514\n',
                          '1\t66265029\t66265030\trs1321172\n',
                          '1\t116545156\t116545157\trs12044852\n',
                          '1\t201015351\t201015352\trs12122721\n',
                          '1\t92515681\t92515682\trs17371561\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92580419\t92580420\trs17380378\n',
                          '1\t24977084\t24977085\trs10903122\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t237115473\t237115474\trs10925318\n',
                          '1\t92543755\t92543756\trs11804321\n'])
        with ZipFile(adb_path) as adb_opened_r:
            fir_lstarts_idx_path = '1/1/lstarts'
            with ZstdFile(adb_opened_r.open(fir_lstarts_idx_path)) as fir_lstarts_idx_opened:
                fir_lstarts_idx_obj = load(fir_lstarts_idx_opened)
                self.assertEqual(fir_lstarts_idx_obj[0],
                                 (8, 8, 8, 9))
            sec_lstarts_idx_path = '1/2/lstarts'
            with ZstdFile(adb_opened_r.open(sec_lstarts_idx_path)) as sec_lstarts_idx_opened:
                sec_lstarts_idx_obj = load(sec_lstarts_idx_opened)
                self.assertEqual(sec_lstarts_idx_obj[0],
                                 (9, 10, 10, 10))
            thi_lstarts_idx_path = '2/1/lstarts'
            with ZstdFile(adb_opened_r.open(thi_lstarts_idx_path)) as thi_lstarts_idx_opened:
                thi_lstarts_idx_obj = load(thi_lstarts_idx_opened)
                self.assertEqual(thi_lstarts_idx_obj[0],
                                 (10, 10, 10, 10))
            fou_lstarts_idx_path = '2/2/lstarts'
            with ZstdFile(adb_opened_r.open(fou_lstarts_idx_path)) as fou_lstarts_idx_opened:
                fou_lstarts_idx_obj = load(fou_lstarts_idx_opened)
                self.assertEqual(fou_lstarts_idx_obj[0],
                                 (10, 10, 10, 10))
        prs_obj = Prs(db_file_path=self.src_file_path,
                      adb_name_prefix='rsids_len',
                      adb_srt_rule=lambda val: val)
        self.assertEqual(list(prs_obj.eq(8)),
                         ['1\t241782991\t241782992\trs952084\n',
                          '1\t86877126\t86877127\trs581405\n',
                          '1\t86876786\t86876787\trs479341\n'])
        self.assertEqual(list(prs_obj.eq(9)),
                         ['1\t154527612\t154527613\trs4131514\n',
                          '1\t66265029\t66265030\trs1321172\n'])
        self.assertEqual(list(prs_obj.eq(10)),
                         ['1\t116545156\t116545157\trs12044852\n',
                          '1\t201015351\t201015352\trs12122721\n',
                          '1\t92515681\t92515682\trs17371561\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92580419\t92580420\trs17380378\n',
                          '1\t24977084\t24977085\trs10903122\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t92543755\t92543756\trs11804321\n',
                          '1\t92516860\t92516861\trs11581176\n',
                          '1\t237115473\t237115474\trs10925318\n',
                          '1\t92543755\t92543756\trs11804321\n'])
        self.assertEqual(list(prs_obj.rng(-9, 9)),
                         ['1\t241782991\t241782992\trs952084\n',
                          '1\t86877126\t86877127\trs581405\n',
                          '1\t86876786\t86876787\trs479341\n',
                          '1\t154527612\t154527613\trs4131514\n',
                          '1\t66265029\t66265030\trs1321172\n'])
        self.assertFalse(list(prs_obj.eq(11)))
        del_files(self.src_file_path,
                  self.db_zst_path,
                  adb_path)


class VcfTests(unittest.TestCase):
    src_vcf = ['##fileformat=VCFv4.1\n',
               '##source=IlluminaPlatinumGenomes, version: 2016-1.0\n',
               '##INFO=<ID=MTD,Number=.,Type=String,Description="Methods that call this site pedigree consistently.">\n',
               '##INFO=<ID=KM,Number=1,Type=String,Description="Minimum normalized k-mer counts">\n',
               '##INFO=<ID=KFP,Number=1,Type=Integer,Description="K-mer GT failures in the main pedigree">\n',
               '##INFO=<ID=KFF,Number=1,Type=Integer,Description="K-mer failures in the founders">\n',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
               '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12877\n',
               'chr1\t126113\t.\tC\tA\t.\tPASS\tMTD=isaac_strelka,bwa_freebayes,bwa_platypus,bwa_gatk;KM=17.27;KFP=0;KFF=0\tGT\t1|1\n',
               'chr1\t567239\t.\tCG\tC\t.\tPASS\tMTD=isaac_strelka,bwa_platypus,bwa_gatk;KM=11.15;KFP=0;KFF=0\tGT\t1|1\n',
               'chr1\t724137\t.\tTAATGG\tTAATGGAATGGAATGGAATGG,TAATGGAATGG\t.\tPASS\tMTD=bwa_gatk;KM=26.29;KFP=0;KFF=0\tGT\t2|1\n',
               'chr1\t725516\t.\tA\tAGAATG\t.\tPASS\tMTD=isaac_strelka;KM=13.37;KFP=0;KFF=0\tGT\t0|1\n',
               'chr1\t753844\t.\tCCT\tC\t.\tPASS\tMTD=isaac_strelka,bwa_platypus;KM=4.21;KFP=0;KFF=0\tGT\t1|0\n',
               'chr1\t761957\t.\tA\tAT\t.\tPASS\tMTD=isaac_strelka,bwa_freebayes,bwa_platypus,bwa_gatk;KM=6.63;KFP=0;KFF=0\tGT\t1|0\n',
               'chr1\t763769\t.\tAT\tA,ATT\t.\tPASS\tMTD=bwa_freebayes,bwa_platypus,bwa_gatk;KM=11.26;KFP=0;KFF=0\tGT\t2|1\n',
               'chr1\t767780\t.\tG\tA\t.\tPASS\tMTD=cgi,bwa_freebayes,bwa_platypus,isaac_strelka,bwa_gatk;KM=11.12;KFP=0;KFF=0\tGT\t1|1\n',
               'chr1\t768116\t.\tAGTTTT\tAGTTTTGTTTT,A\t.\tPASS\tMTD=bwa_freebayes,bwa_platypus,bwa_gatk;KM=18.37;KFP=0;KFF=0\tGT\t1|2\n',
               'chr1\t769138\t.\tCAT\tC\t.\tPASS\tMTD=isaac_strelka,bwa_freebayes,bwa_platypus,bwa_gatk;KM=12.00;KFP=0;KFF=0\tGT\t0|1\n',
               'chr14\t56412076\t.\tG\tGCATACATA\t.\tPASS\tMTD=isaac_strelka,bwa_platypus,bwa_gatk;KM=28.62;KFP=0;KFF=0\tGT\t1|1\n',
               'chr14\t56422021\t.\tAAAAC\tA\t.\tPASS\tMTD=isaac_strelka,bwa_freebayes,bwa_gatk;KM=9.43;KFP=0;KFF=0\tGT\t0|1\n',
               'chr14\t56551760\t.\tC\tCTATCTATA,CTATCTATCTATA\t.\tPASS\tMTD=bwa_gatk;KM=5.77;KFP=0;KFF=0\tGT\t1|2\n',
               'chr14\t56564010\t.\tAACACACACAC\tA,AAC\t.\tPASS\tMTD=bwa_freebayes;KM=38.71;KFP=0;KFF=0\tGT\t2|1\n',
               'chr14\t56664634\t.\tGACACACAC\tGAC\t.\tPASS\tMTD=bwa_freebayes,bwa_platypus,bwa_gatk;KM=6.83;KFP=0;KFF=0\tGT\t0|1\n',
               'chr14\t56669715\t.\tCA\tC\t.\tPASS\tMTD=isaac_strelka,bwa_freebayes,bwa_platypus,bwa_gatk;KM=10.85;KFP=0;KFF=0\tGT\t0|1\n',
               'chr14\t56783534\t.\tT\tC,G\t.\tPASS\tMTD=isaac_strelka,bwa_platypus,bwa_gatk;KM=8.95;KFP=0;KFF=0\tGT\t1|2\n',
               'chr14\t56868236\t.\tTA\tTAA,T\t.\tPASS\tMTD=bwa_freebayes,bwa_platypus,bwa_gatk;KM=4.00;KFP=0;KFF=0\tGT\t1|2\n',
               'chr14\t56898904\t.\tTTTCC\tTTTCCTTCC\t.\tPASS\tMTD=bwa_freebayes,bwa_gatk;KM=21.83;KFP=0;KFF=0\tGT\t0|1\n',
               'chr14\t57002112\t.\tAAGAG\tAAGAGAGAG,A\t.\tPASS\tMTD=bwa_gatk;KM=15.71;KFP=0;KFF=0\tGT\t2|1\n']
    src_file_path = os.path.join(os.getcwd(),
                                 'vcf.vcf')
    db_zst_path = os.path.join(os.getcwd(),
                               'vcf.vcf.zst')

    def test_pseudo_tabix(self):
        adb_path = os.path.join(os.getcwd(),
                                'vcf.vcf.tabix.adb')
        with open(self.src_file_path, 'w') as src_file_opened:
            for src_vcf_line in self.src_vcf:
                src_file_opened.write(src_vcf_line)
        del_files(self.db_zst_path,
                  adb_path)

        def get_intervals(vcf_line: str):
            vcf_row = vcf_line.split('\t')
            intervals = []
            for ref_lett_ind in range(len(vcf_row[3])):
                chrom = vcf_row[0]
                pos = int(vcf_row[1]) + ref_lett_ind
                intervals.append([chrom, pos])
            return tuple(intervals)

        idx_obj = Idx(db_file_path=self.src_file_path,
                      adb_name_prefix='tabix',
                      db_line_prs=get_intervals,
                      adb_srt_rule=SrtRules.natur)
        idx_obj.idx()
        prs_obj = Prs(db_file_path=self.src_file_path,
                      adb_name_prefix='tabix',
                      adb_srt_rule=SrtRules.natur)

        def trun_vcf_line(vcf_line):
            return vcf_line.split('\tMTD')[0]

        self.assertEqual(list(map(trun_vcf_line,
                                  prs_obj.eq(['chr1', 126113],
                                             ['chr1', 567239],
                                             ['chr1', 567240]))),
                         ['chr1\t126113\t.\tC\tA\t.\tPASS',
                          'chr1\t567239\t.\tCG\tC\t.\tPASS',
                          'chr1\t567239\t.\tCG\tC\t.\tPASS'])
        self.assertFalse(list(prs_obj.eq(['chr1', 567241])))
        self.assertEqual(list(map(trun_vcf_line,
                                  prs_obj.eq(['chr1', 768116]))),
                         ['chr1\t768116\t.\tAGTTTT\tAGTTTTGTTTT,A\t.\tPASS'])
        self.assertEqual(list(map(trun_vcf_line,
                                  prs_obj.eq(['chr1', 768121]))),
                         ['chr1\t768116\t.\tAGTTTT\tAGTTTTGTTTT,A\t.\tPASS'])
        self.assertEqual(list(map(trun_vcf_line,
                                  prs_obj.rng(['chr1', 0],
                                              ['chr1', 126113]))),
                         ['chr1\t126113\t.\tC\tA\t.\tPASS'])
        self.assertEqual(list(map(trun_vcf_line,
                                  prs_obj.rng(['chr14', 57002112],
                                              ['chr14', inf]))),
                         ['chr14\t57002112\t.\tAAGAG\tAAGAGAGAG,A\t.\tPASS'
                          for dupl_ind in range(5)])
        self.assertEqual(list(map(trun_vcf_line,
                                  prs_obj.rng(['chr14', 56669710],
                                              ['chr14', 56868240]))),
                         ['chr14\t56669715\t.\tCA\tC\t.\tPASS',
                          'chr14\t56669715\t.\tCA\tC\t.\tPASS',
                          'chr14\t56783534\t.\tT\tC,G\t.\tPASS',
                          'chr14\t56868236\t.\tTA\tTAA,T\t.\tPASS',
                          'chr14\t56868236\t.\tTA\tTAA,T\t.\tPASS'])
        self.assertEqual(list(map(trun_vcf_line,
                                  prs_obj.rng(['chr1', 769138],
                                              ['chr14', 56412076]))),
                         ['chr1\t769138\t.\tCAT\tC\t.\tPASS',
                          'chr1\t769138\t.\tCAT\tC\t.\tPASS',
                          'chr1\t769138\t.\tCAT\tC\t.\tPASS',
                          'chr14\t56412076\t.\tG\tGCATACATA\t.\tPASS'])
        self.assertFalse(list(prs_obj.rng(['chr1', 126110],
                                          ['chr1', 126112])))
        self.assertFalse(list(prs_obj.rng(['chr1', 57002113],
                                          ['chr1', 57002115])))
        del_files(self.src_file_path,
                  self.db_zst_path,
                  adb_path)


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
                         [[10, '\t', 11]])
        self.assertEqual(self.srt_rules.natur('rs10\t11'),
                         [[float('+inf'), 'rs', 10, '\t', 11]])
        self.assertEqual(self.srt_rules.natur('10,11.1'),
                         [[10, ',', 11.1]])
        self.assertEqual(self.srt_rules.natur('10.1,11'),
                         [[10.1, ',', 11]])
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
                         [[float('+inf'), 'qwerty\t', 10.0]])
        self.assertEqual(self.srt_rules.natur('0,1e+2\tqwerty',
                                              dec_delimiter=','),
                         [[10.0, '\tqwerty']])
        self.assertEqual(self.srt_rules.natur('1,1, -2,2, str',
                                              dec_delimiter=','),
                         [[1.1, ', ', -2.2, ', str']])
        self.assertEqual(self.srt_rules.natur(['1,1', '-2,2', 'str'],
                                              dec_delimiter=',',
                                              nums_first=False),
                         [[1.1], [-2.2], [float('-inf'), 'str']])
        self.assertEqual(self.srt_rules.natur(['chr14', 1, 10]),
                         [[float('+inf'), 'chr', 14], [1], [10]])

    def test_letts_nums_srt_rule(self):
        self.assertEqual(self.srt_rules.letts_nums('rs1'),
                         ['rs', 1])
        self.assertEqual(self.srt_rules.letts_nums('rs010'),
                         ['rs', 10])
        self.assertEqual(self.srt_rules.letts_nums('ENSG000'),
                         ['ENSG', 0])
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


if __name__ == "__main__":
    unittest.main()
