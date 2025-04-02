# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import os
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
from .err import QueryStartGtEndError
from pyzstd import (SeekableZstdFile,
                    ZstdFile)

if __name__ == 'main':
    __version__ = 'v6.1.0'
    __authors__ = [{'name': 'Platon Bykadorov',
                    'email': 'platon.work@gmail.com',
                    'years': '2023-2025'}]


class Prs(Idx):
    def __init__(self,
                 db_file_path: str,
                 adb_name_prefix: str,
                 adb_srt_rule: Callable,
                 adb_srt_rule_kwargs: None | dict = None):
        super().__init__(db_file_path=db_file_path,
                         adb_name_prefix=adb_name_prefix,
                         db_line_prs=None,
                         adb_srt_rule=adb_srt_rule,
                         adb_srt_rule_kwargs=adb_srt_rule_kwargs)
        self.adb_opened_r = ZipFile(self.adb_path)
        self.db_zst_opened_r = TextIOWrapper(SeekableZstdFile(self.db_zst_path))

    def prep_query(self,
                   query_start: Any,
                   query_end: Any = None) -> list[Any,
                                                  Any]:
        if not query_end:
            query_end = query_start
        prepd_query_start = self.adb_srt_rule(query_start,
                                              **self.adb_srt_rule_kwargs)
        prepd_query_end = self.adb_srt_rule(query_end,
                                            **self.adb_srt_rule_kwargs)
        if prepd_query_start > prepd_query_end:
            raise QueryStartGtEndError(prepd_query_start,
                                       prepd_query_end)
        prepd_query_bords = [prepd_query_start,
                             prepd_query_end]
        return prepd_query_bords

    def walk_dir_tree(self,
                      prepd_query_bords: list[Any,
                                              Any],
                      any_idx_path: str = 'paths') -> Generator:
        if os.path.basename(any_idx_path) == 'lstarts':
            yield any_idx_path
        else:
            with self.adb_opened_r.open(any_idx_path) as paths_idx_opened:
                paths_idx_obj = load(paths_idx_opened)
            start_gchi_any_idx_ind = bisect_left(paths_idx_obj[0],
                                                 prepd_query_bords[0]) - 1
            if start_gchi_any_idx_ind < 0:
                start_gchi_any_idx_ind = 0
            end_gchi_any_idx_ind = bisect_right(paths_idx_obj[0],
                                                prepd_query_bords[1]) - 1
            if end_gchi_any_idx_ind >= 0:
                for neces_gchi_any_idx_path in paths_idx_obj[1][start_gchi_any_idx_ind:
                                                                end_gchi_any_idx_ind + 1]:
                    for neces_lstarts_idx_path in self.walk_dir_tree(prepd_query_bords,
                                                                     neces_gchi_any_idx_path):
                        yield neces_lstarts_idx_path

    def read_lstarts_idx(self,
                         lstarts_idx_path: str) -> list:
        with ZstdFile(self.adb_opened_r.
                      open(lstarts_idx_path)) as lstarts_idx_opened:
            lstarts_idx = load(lstarts_idx_opened)
            return lstarts_idx

    def eq(self,
           *queries: Any) -> Generator:
        for query in queries:
            prepd_query_bords = self.prep_query(query)
            for neces_lstarts_idx_path in self.walk_dir_tree(prepd_query_bords):
                neces_lstarts_idx_obj = self.read_lstarts_idx(neces_lstarts_idx_path)
                start_lstart_ind = bisect_left(neces_lstarts_idx_obj[0],
                                               prepd_query_bords[0])
                if start_lstart_ind == len(neces_lstarts_idx_obj[0]) \
                        or prepd_query_bords[0] != neces_lstarts_idx_obj[0][start_lstart_ind]:
                    continue
                end_lstart_ind = bisect_right(neces_lstarts_idx_obj[0],
                                              prepd_query_bords[1]) - 1
                if prepd_query_bords[1] != neces_lstarts_idx_obj[0][end_lstart_ind]:
                    continue
                for lstart_ind in range(start_lstart_ind,
                                        end_lstart_ind + 1):
                    self.db_zst_opened_r.seek(neces_lstarts_idx_obj[1][lstart_ind])
                    found_line = self.db_zst_opened_r.readline()
                    yield found_line

    def rng(self,
            query_start: Any,
            query_end: Any) -> Generator:
        prepd_query_bords = self.prep_query(query_start,
                                            query_end)
        for neces_lstarts_idx_path in self.walk_dir_tree(prepd_query_bords):
            neces_lstarts_idx_obj = self.read_lstarts_idx(neces_lstarts_idx_path)
            if prepd_query_bords[0] <= neces_lstarts_idx_obj[0][0]:
                start_lstart_ind = 0
            else:
                start_lstart_ind = bisect_left(neces_lstarts_idx_obj[0],
                                               prepd_query_bords[0])
            neces_lstarts_quan = len(neces_lstarts_idx_obj[0])
            if start_lstart_ind == neces_lstarts_quan:
                continue
            if neces_lstarts_idx_obj[0][-1] <= prepd_query_bords[1]:
                end_lstart_ind = neces_lstarts_quan - 1
            else:
                end_lstart_ind = bisect_right(neces_lstarts_idx_obj[0],
                                              prepd_query_bords[1]) - 1
            for lstart_ind in range(start_lstart_ind,
                                    end_lstart_ind + 1):
                self.db_zst_opened_r.seek(neces_lstarts_idx_obj[1][lstart_ind])
                found_line = self.db_zst_opened_r.readline()
                yield found_line
