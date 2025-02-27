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
from .idx import Idx
from .err import (NoIdxsError,
                  QueryStartGtEndError)
from pyzstd import (SeekableZstdFile,
                    ZstdFile)

__version__ = 'v3.0.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023-2025'}]


class Prs(Idx):
    def __init__(self,
                 db_file_path: str,
                 idx_prefix: str,
                 srt_rule: None | Callable = None,
                 srt_rule_kwargs: None | dict = None,
                 srt_rule_cols_delimiter: str | None = '\t',
                 srt_rule_col_inds: None | int | list | tuple = None):
        super().__init__(db_file_path,
                         idx_prefix,
                         your_line_prs=None,
                         srt_rule=srt_rule,
                         srt_rule_kwargs=srt_rule_kwargs,
                         srt_rule_cols_delimiter=srt_rule_cols_delimiter,
                         srt_rule_col_inds=srt_rule_col_inds)
        self.adb_opened_r = ZipFile(self.adb_path)
        self.db_zst_opened_r = TextIOWrapper(SeekableZstdFile(self.db_zst_path))
        self.idx_names = list(filter(lambda name:
                                     name.endswith('.idx'),
                                     self.adb_opened_r.namelist()))
        if not self.idx_names:
            raise NoIdxsError()
        self.idx_begins = sorted(map(lambda idx_name:
                                     eval(idx_name[:-4]),
                                     self.idx_names))
        self.meta = self.read_meta()

    def read_meta(self):
        with TextIOWrapper(self.adb_opened_r.open('meta.txt')) as meta_opened:
            meta = {}
            for key_val in meta_opened.readlines():
                key, val = key_val[:-1].split('=')
                try:
                    meta[key] = eval(val)
                except (NameError,
                        SyntaxError):
                    meta[key] = val
        return meta

    def prep_query(self,
                   query_start: Any,
                   query_end: Any = None) -> tuple[list[Any,
                                                        Any],
                                                   list]:
        if not query_end:
            query_end = query_start
        prepd_query_start = self.srt_rule(query_start,
                                          **self.srt_rule_kwargs)
        prepd_query_end = self.srt_rule(query_end,
                                        **self.srt_rule_kwargs)
        if prepd_query_start > prepd_query_end:
            raise QueryStartGtEndError(prepd_query_start,
                                       prepd_query_end)
        prepd_query_bords = [prepd_query_start,
                             prepd_query_end]
        start_idx_ind = bisect_left(self.idx_begins,
                                    prepd_query_start) - 1
        if start_idx_ind < 0:
            start_idx_ind = 0
        end_idx_ind = bisect_right(self.idx_begins,
                                   prepd_query_end) - 1
        if end_idx_ind < 0:
            neces_idx_names = []
        else:
            neces_idx_begins = self.idx_begins[start_idx_ind:end_idx_ind + 1]
            neces_idx_names = [f'{neces_idx_begin}.idx'
                               for neces_idx_begin in neces_idx_begins]
        return (prepd_query_bords,
                neces_idx_names)

    def eq(self,
           *queries: Any) -> Generator:
        for query in queries:
            prepd_query_bords, neces_idx_names = self.prep_query(query)
            for neces_idx_name in neces_idx_names:
                with ZstdFile(self.adb_opened_r.open(neces_idx_name)) as neces_idx_opened:
                    neces_idx = load(neces_idx_opened)
                with ZstdFile(self.adb_opened_r.open(f'{neces_idx_name[:-4]}.b')) as neces_b_opened:
                    neces_line_starts = load(neces_b_opened)
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
        prepd_query_bords, neces_idx_names = self.prep_query(query_start,
                                                             query_end)
        for neces_idx_name in neces_idx_names:
            with ZstdFile(self.adb_opened_r.open(neces_idx_name)) as neces_idx_opened:
                neces_idx = load(neces_idx_opened)
            with ZstdFile(self.adb_opened_r.open(f'{neces_idx_name[:-4]}.b')) as neces_b_opened:
                neces_line_starts = load(neces_b_opened)
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
