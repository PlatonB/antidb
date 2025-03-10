# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
from typing import (Callable,
                    Any,
                    Generator)
from zipfile import ZipFile
from pickle import load
from io import TextIOWrapper
from math import inf
from bisect import (bisect_left,
                    bisect_right)
from idx import Idx
from err import (NoIdxsError,
                 QueryStartGtEndError)
from pyzstd import (SeekableZstdFile,
                    ZstdFile)

if __name__ == 'main':
    __version__ = 'v3.5.0'
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
        self.adb_opened_r = ZipFile(self.adb_path)
        self.db_zst_opened_r = TextIOWrapper(SeekableZstdFile(self.db_zst_path))
        self.idx_names = list(filter(lambda name:
                                     name.endswith('.idx'),
                                     self.adb_opened_r.namelist()))
        if not self.idx_names:
            raise NoIdxsError()
        self.idx_begins = sorted(map(lambda idx_name:
                                     eval(idx_name[:-4].split('.')[0]),
                                     self.idx_names))

    def read_idx(self,
                 idx_name: str) -> list:
        with ZstdFile(self.adb_opened_r.open(idx_name)) as idx_opened:
            idx = load(idx_opened)
            return idx

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

    def select_idx_names(self,
                         prepd_query_bords: list[Any,
                                                 Any]) -> list:
        start_idx_ind = bisect_left(self.idx_begins,
                                    prepd_query_bords[0]) - 1
        if start_idx_ind < 0:
            start_idx_ind = 0
        end_idx_ind = bisect_right(self.idx_begins,
                                   prepd_query_bords[1]) - 1
        if end_idx_ind < 0:
            neces_idx_names = []
        else:
            neces_idx_begins = self.idx_begins[start_idx_ind:end_idx_ind + 1]
            neces_idx_names = []
            for neces_idx_begin in neces_idx_begins:
                name_dupl_num = 0
                while True:
                    if type(neces_idx_begin) is str:
                        neces_idx_name = f"'{neces_idx_begin}'.{name_dupl_num}.idx"
                    else:
                        neces_idx_name = f'{neces_idx_begin}.{name_dupl_num}.idx'
                    if neces_idx_name in self.idx_names:
                        neces_idx_names.append(neces_idx_name)
                        name_dupl_num += 1
                    else:
                        break
        return neces_idx_names

    def eq(self,
           *queries: Any) -> Generator:
        for query in queries:
            prepd_query_bords = self.prep_query(query)
            neces_idx_names = self.select_idx_names(prepd_query_bords)
            for neces_idx_name in neces_idx_names:
                neces_idx = self.read_idx(neces_idx_name)
                neces_line_starts = self.read_idx(f'{neces_idx_name[:-4]}.b')
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
                    self.db_zst_opened_r.seek(neces_line_starts[line_start_ind])
                    found_line = self.db_zst_opened_r.readline()
                    yield found_line

    def rng(self,
            query_start: Any,
            query_end: Any) -> Generator:
        prepd_query_bords = self.prep_query(query_start,
                                            query_end)
        neces_idx_names = self.select_idx_names(prepd_query_bords)
        for neces_idx_name in neces_idx_names:
            neces_idx = self.read_idx(neces_idx_name)
            neces_line_starts = self.read_idx(f'{neces_idx_name[:-4]}.b')
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
