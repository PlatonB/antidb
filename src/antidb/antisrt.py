import re
import os
from heapq import merge
from functools import partial

__version__ = 'v1.0.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023'}]


class SrtRules():
    @staticmethod
    def natur(src_file_line,
              delimiter='\t',
              src_file_colinds=None):
        src_file_row = src_file_line.rstrip().split(delimiter)
        if type(src_file_colinds) is int:
            src_file_row = [src_file_row[src_file_colinds]]
        elif type(src_file_colinds) in [list, tuple]:
            src_file_row = [src_file_row[src_file_colind]
                            for src_file_colind in src_file_colinds]
        spl_file_row = []
        for src_file_cell in src_file_row:
            src_file_subcells = list(filter(lambda src_file_subcell:
                                            src_file_subcell,
                                            re.split(r'(\d+)',
                                                     src_file_cell)))
            for src_file_subcell_ind in range(len(src_file_subcells)):
                try:
                    src_file_subcells[src_file_subcell_ind] = int(src_file_subcells[src_file_subcell_ind])
                except ValueError:
                    pass
            if type(src_file_subcells[0]) is str:
                src_file_subcells.insert(0, float('+inf'))
            spl_file_row.append(src_file_subcells)
        return spl_file_row


class Srt():
    def __init__(self,
                 src_file_path,
                 srt_rule,
                 **srt_rule_kwargs):
        self.src_file_path = os.path.normpath(src_file_path)
        self.presrtd_file_paths = []
        self.srt_rule = srt_rule
        if srt_rule_kwargs:
            self.srt_rule_kwargs = srt_rule_kwargs

    @staticmethod
    def iter_file(file_path):
        with open(file_path) as file_opened:
            for file_line in file_opened:
                yield file_line

    def pre_srt(self,
                chunk_elems_quan=10000000):
        with open(self.src_file_path) as src_file_opened:
            while True:
                src_file_lstart = src_file_opened.tell()
                if not src_file_opened.readline().startswith('#'):
                    src_file_opened.seek(src_file_lstart)
                    break
            chunk, chunk_len, chunk_num = [], 0, 0
            for src_file_line in src_file_opened:
                chunk.append(src_file_line)
                chunk_len += 1
                if chunk_len == chunk_elems_quan:
                    chunk_num += 1
                    presrtd_file_path = f'{self.src_file_path}.{chunk_num}'
                    self.presrtd_file_paths.append(presrtd_file_path)
                    chunk.sort(key=partial(self.srt_rule,
                                           **self.srt_rule_kwargs))
                    with open(presrtd_file_path, mode='w') as presrtd_file_opened:
                        for presrtd_file_line in chunk:
                            presrtd_file_opened.write(presrtd_file_line)
                    chunk, chunk_len = [], 0
            if chunk:
                chunk_num += 1
                presrtd_file_path = f'{self.src_file_path}.{chunk_num}'
                self.presrtd_file_paths.append(presrtd_file_path)
                chunk.sort(key=partial(self.srt_rule,
                                       **self.srt_rule_kwargs))
                with open(presrtd_file_path, mode='w') as presrtd_file_opened:
                    for presrtd_file_line in chunk:
                        presrtd_file_opened.write(presrtd_file_line)

    def mrg_srt(self,
                mrgd_file_suff='srtd',
                del_presrtd_files=True):
        if not self.presrtd_file_paths:
            return None
        presrtd_file_common_path = re.sub(r'\.\d+$',
                                          '',
                                          self.presrtd_file_paths[0])
        mrgd_file_path = f'{presrtd_file_common_path}.{mrgd_file_suff}'
        if len(self.presrtd_file_paths) == 1:
            os.rename(self.presrtd_file_paths[0],
                      mrgd_file_path)
            return mrgd_file_path
        with open(mrgd_file_path, mode='w') as mrgd_file_opened:
            for mrgd_file_line in merge(*map(self.iter_file,
                                             self.presrtd_file_paths),
                                        key=partial(self.srt_rule,
                                                    **self.srt_rule_kwargs)):
                mrgd_file_opened.write(mrgd_file_line)
        if del_presrtd_files:
            for presrtd_file_path in self.presrtd_file_paths:
                os.remove(presrtd_file_path)
        return mrgd_file_path
