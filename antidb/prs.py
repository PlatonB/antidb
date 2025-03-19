# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import os
from typing import (Callable,
                    Any,
                    Generator)
from zipfile import (ZipFile,
                     Path)
from pickle import load
from io import TextIOWrapper
from math import inf
from bisect import (bisect_left,
                    bisect_right)
from copy import deepcopy
from .idx import Idx
from .err import QueryStartGtEndError
from pyzstd import (SeekableZstdFile,
                    ZstdFile)

if __name__ == 'main':
    __version__ = 'v5.1.1'
    __authors__ = [{'name': 'Platon Bykadorov',
                    'email': 'platon.work@gmail.com',
                    'years': '2023-2025'}]


class Prs(Idx):
    def __init__(self,
                 db_file_path: str,
                 idx_name_prefix: str,
                 idx_srt_rule: Callable,
                 idx_srt_rule_kwargs: None | dict = None):
        super().__init__(db_file_path=db_file_path,
                         idx_name_prefix=idx_name_prefix,
                         db_line_prs=None,
                         idx_srt_rule=idx_srt_rule,
                         idx_srt_rule_kwargs=idx_srt_rule_kwargs)
        self.adb_path_obj = Path(self.adb_path)
        self.adb_opened_r = ZipFile(self.adb_path)
        self.db_zst_opened_r = TextIOWrapper(SeekableZstdFile(self.db_zst_path))

    def prep_query(self,
                   query_start: Any,
                   query_end: Any = None) -> list[Any,
                                                  Any]:
        if not query_end:
            query_end = query_start
        prepd_query_start = self.idx_srt_rule(query_start,
                                              **self.idx_srt_rule_kwargs)
        prepd_query_end = self.idx_srt_rule(query_end,
                                            **self.idx_srt_rule_kwargs)
        if prepd_query_start > prepd_query_end:
            raise QueryStartGtEndError(prepd_query_start,
                                       prepd_query_end)
        prepd_query_bords = [prepd_query_start,
                             prepd_query_end]
        return prepd_query_bords

    def sel_dir_names(self,
                      prepd_query_bords: list[Any,
                                              Any],
                      cur_dir_path: str = '') -> Generator:
        chi_path_objs = list(self.adb_path_obj.
                             joinpath(cur_dir_path).
                             iterdir())
        if chi_path_objs[0].is_file():
            neces_idx_path = chi_path_objs[0].at
            yield neces_idx_path
        else:
            chi_chunk_begins = sorted(map(lambda chi_path_obj:
                                          eval(chi_path_obj.stem),
                                          chi_path_objs))
            start_dir_ind = bisect_left(chi_chunk_begins,
                                        prepd_query_bords[0]) - 1
            if start_dir_ind < 0:
                start_dir_ind = 0
            end_dir_ind = bisect_right(chi_chunk_begins,
                                       prepd_query_bords[1]) - 1
            if end_dir_ind >= 0:
                neces_chunk_begins = chi_chunk_begins[start_dir_ind:
                                                      end_dir_ind + 1]
                prev_neces_chunk_begin = None
                for neces_chunk_begin in neces_chunk_begins:
                    if neces_chunk_begin != prev_neces_chunk_begin:
                        name_dupl_num = 1
                        prev_neces_chunk_begin = deepcopy(neces_chunk_begin)
                    else:
                        name_dupl_num += 1
                    if type(neces_chunk_begin) is str:
                        neces_dir_name = f"'{neces_chunk_begin}'.{name_dupl_num}/"
                    else:
                        neces_dir_name = f'{neces_chunk_begin}.{name_dupl_num}/'
                    neces_dir_path = os.path.join(cur_dir_path,
                                                  neces_dir_name)
                    for neces_idx_path in self.sel_dir_names(prepd_query_bords,
                                                             neces_dir_path):
                        yield neces_idx_path

    def read_idx(self,
                 idx_path: str) -> list:
        with ZstdFile(self.adb_opened_r.open(idx_path)) as idx_opened:
            idx = load(idx_opened)
            return idx

    def eq(self,
           *queries: Any) -> Generator:
        for query in queries:
            prepd_query_bords = self.prep_query(query)
            for neces_idx_path in self.sel_dir_names(prepd_query_bords):
                neces_idx = self.read_idx(neces_idx_path)
                neces_lstarts = self.read_idx(f'{neces_idx_path[:-4]}.b')
                start_idxval_ind = bisect_left(neces_idx,
                                               prepd_query_bords[0])
                if start_idxval_ind == len(neces_idx) \
                        or prepd_query_bords[0] != neces_idx[start_idxval_ind]:
                    continue
                end_idxval_ind = bisect_right(neces_idx,
                                              prepd_query_bords[1]) - 1
                if prepd_query_bords[1] != neces_idx[end_idxval_ind]:
                    continue
                for line_start_ind in range(start_idxval_ind,
                                            end_idxval_ind + 1):
                    self.db_zst_opened_r.seek(neces_lstarts[line_start_ind])
                    found_line = self.db_zst_opened_r.readline()
                    yield found_line

    def rng(self,
            query_start: Any,
            query_end: Any) -> Generator:
        prepd_query_bords = self.prep_query(query_start,
                                            query_end)
        for neces_idx_path in self.sel_dir_names(prepd_query_bords):
            neces_idx = self.read_idx(neces_idx_path)
            neces_line_starts = self.read_idx(f'{neces_idx_path[:-4]}.b')
            start_idxval_ind = bisect_left(neces_idx,
                                           prepd_query_bords[0])
            if start_idxval_ind == len(neces_idx):
                continue
            end_idxval_ind = bisect_right(neces_idx,
                                          prepd_query_bords[1]) - 1
            for line_start_ind in range(start_idxval_ind,
                                        end_idxval_ind + 1):
                self.db_zst_opened_r.seek(neces_line_starts[line_start_ind])
                found_line = self.db_zst_opened_r.readline()
                yield found_line
