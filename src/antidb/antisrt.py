import re
import os
from typing import (Callable,
                    Any)
from inspect import stack
from heapq import merge
from functools import partial

__version__ = 'v2.1.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023'}]


class DelimitersMatchError(Exception):
    def __init__(self,
                 cols_delimiter,
                 dec_delimiter):
        err_msg = f'''\nColums delimiter ({cols_delimiter})
matches decimal delimiter ({dec_delimiter})'''
        super().__init__(err_msg)


class NoSrcFilesError(Exception):
    def __init__(self,
                 func_name):
        err_msg = f'''\nThere are no source file(s)
for {func_name} function/method'''
        super().__init__(err_msg)


class SrtRules():
    @staticmethod
    def natur(src_file_line: str,
              cols_delimiter: str = '\t',
              dec_delimiter: str = '.',
              col_inds: None | int | list | tuple = None,
              nums_first: bool = True) -> list:
        if cols_delimiter == dec_delimiter:
            raise DelimitersMatchError(cols_delimiter,
                                       dec_delimiter)
        src_file_row = src_file_line.rstrip().split(cols_delimiter)
        if type(col_inds) is int:
            src_file_row = [src_file_row[col_inds]]
        elif type(col_inds) in [list, tuple]:
            src_file_row = [src_file_row[col_ind]
                            for col_ind in col_inds]
        if dec_delimiter == '.':
            split_cell = r'(-?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)'
        elif dec_delimiter == ',':
            split_cell = r'(-?\d+(?:,\d*)?(?:[Ee][+-]?\d+)?)'
        spl_file_row = []
        for cell in src_file_row:
            subcells = list(filter(lambda subcell:
                                   subcell,
                                   re.split(split_cell,
                                            cell)))
            for subcell_ind in range(len(subcells)):
                try:
                    subcells[subcell_ind] = int(subcells[subcell_ind])
                except ValueError:
                    try:
                        subcells[subcell_ind] = float(subcells[subcell_ind])
                    except ValueError:
                        if dec_delimiter == ',':
                            try:
                                subcells[subcell_ind] = float(subcells[subcell_ind].replace(',', '.'))
                            except ValueError:
                                pass
            if type(subcells[0]) is str:
                if nums_first:
                    subcells.insert(0, float('+inf'))
                else:
                    subcells.insert(0, float('-inf'))
            spl_file_row.append(subcells)
        return spl_file_row


class Srt(SrtRules):
    def __init__(self,
                 unsrtd_file_path: None | str = None,
                 presrtd_file_paths: None | list = None,
                 srt_rule: None | Callable = None,
                 **srt_rule_kwargs: Any):
        if unsrtd_file_path:
            self.unsrtd_file_path = os.path.normpath(unsrtd_file_path)
        else:
            self.unsrtd_file_path = None
        if presrtd_file_paths:
            self.presrtd_file_paths = presrtd_file_paths
        else:
            self.presrtd_file_paths = []
        if srt_rule:
            self.srt_rule = srt_rule
        else:
            self.srt_rule = self.natur
        if srt_rule_kwargs:
            self.srt_rule_kwargs = srt_rule_kwargs
        else:
            self.srt_rule_kwargs = {}

    @staticmethod
    def iter_file(file_path: str) -> str:
        with open(file_path) as file_opened:
            for file_line in file_opened:
                yield file_line

    def pre_srt(self,
                chunk_elems_quan: int = 10000000):
        if not self.unsrtd_file_path:
            raise NoSrcFilesError(stack()[0].function)
        self.presrtd_file_paths = []
        with open(self.unsrtd_file_path) as src_file_opened:
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
                    presrtd_file_path = f'{self.unsrtd_file_path}.{chunk_num}'
                    self.presrtd_file_paths.append(presrtd_file_path)
                    chunk.sort(key=partial(self.srt_rule,
                                           **self.srt_rule_kwargs))
                    with open(presrtd_file_path, mode='w') as presrtd_file_opened:
                        for presrtd_file_line in chunk:
                            presrtd_file_opened.write(presrtd_file_line)
                    chunk, chunk_len = [], 0
            if chunk:
                chunk_num += 1
                presrtd_file_path = f'{self.unsrtd_file_path}.{chunk_num}'
                self.presrtd_file_paths.append(presrtd_file_path)
                chunk.sort(key=partial(self.srt_rule,
                                       **self.srt_rule_kwargs))
                with open(presrtd_file_path, mode='w') as presrtd_file_opened:
                    for presrtd_file_line in chunk:
                        presrtd_file_opened.write(presrtd_file_line)

    def mrg_srt(self,
                mrgd_file_suff: str = 'srtd',
                del_presrtd_files: bool = True) -> None | str:
        if not self.presrtd_file_paths:
            raise NoSrcFilesError(stack()[0].function)
        presrtd_file_common_path = re.sub(r'\.\d+$',
                                          '',
                                          self.presrtd_file_paths[0])
        mrgd_file_path = f'{presrtd_file_common_path}.{mrgd_file_suff}'
        if len(self.presrtd_file_paths) == 1:
            os.rename(self.presrtd_file_paths[0],
                      mrgd_file_path)
            self.presrtd_file_paths = []
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
            self.presrtd_file_paths = []
        return mrgd_file_path
