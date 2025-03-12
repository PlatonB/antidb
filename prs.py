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
from idx import Idx
from err import (NoIdxsError,
                 QueryStartGtEndError)
from pyzstd import (SeekableZstdFile,
                    ZstdFile)

if __name__ == 'main':
    __version__ = 'v4.0.0'
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
        self.top_dir_names = [path_obj.at for path_obj in Path(self.adb_path).iterdir()]
        if not self.top_dir_names:
            raise NoIdxsError()
        self.top_chunk_begins = sorted(map(lambda dir_name:
                                           eval(dir_name.rsplit('.', maxsplit=1)[0]),
                                           self.top_dir_names))

    def read_idx(self,
                 idx_path: str) -> list:
        with ZstdFile(self.adb_opened_r.open(idx_path)) as idx_opened:
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

    def sel_top_dir_names(self,
                          prepd_query_bords: list[Any,
                                                  Any]) -> Generator:
        start_top_dir_ind = bisect_left(self.top_chunk_begins,
                                        prepd_query_bords[0]) - 1
        if start_top_dir_ind < 0:
            start_top_dir_ind = 0
        end_top_dir_ind = bisect_right(self.top_chunk_begins,
                                       prepd_query_bords[1]) - 1
        if end_top_dir_ind >= 0:
            neces_top_chunk_begins = self.top_chunk_begins[start_top_dir_ind:
                                                           end_top_dir_ind + 1]
            for neces_top_chunk_begin in neces_top_chunk_begins:
                name_dupl_num = 1
                while True:
                    if type(neces_top_chunk_begin) is str:
                        neces_top_dir_name = f"'{neces_top_chunk_begin}'.{name_dupl_num}/"
                    else:
                        neces_top_dir_name = f'{neces_top_chunk_begin}.{name_dupl_num}/'
                    if neces_top_dir_name in self.top_dir_names:
                        yield neces_top_dir_name
                        name_dupl_num += 1
                    else:
                        break

    def sel_idx_paths(self,
                      neces_dir_path: str,
                      prepd_query_bords: list[Any,
                                              Any]) -> Generator:
        child_paths = [path_obj.at
                       for path_obj
                       in Path(self.adb_path,
                               neces_dir_path).iterdir()]
        if Path(self.adb_path,
                child_paths[0]).is_file():
            neces_idx_path = sorted(child_paths)[1]
            yield neces_idx_path
        else:
            child_dir_path_1 = child_paths[0]
            child_dir_path_2 = child_paths[1]
            child_dir_name_1 = os.path.basename(child_dir_path_1[:-1])
            child_dir_name_2 = os.path.basename(child_dir_path_2[:-1])
            chunk_begin_1 = eval(child_dir_name_1.rsplit('.', maxsplit=1)[0])
            chunk_begin_2 = eval(child_dir_name_2.rsplit('.', maxsplit=1)[0])
            if prepd_query_bords[0] <= chunk_begin_1 <= prepd_query_bords[1]:
                for neces_idx_path in self.sel_idx_paths(child_dir_path_1,
                                                         prepd_query_bords):
                    yield neces_idx_path
            if prepd_query_bords[0] <= chunk_begin_2 <= prepd_query_bords[1]:
                for neces_idx_path in self.sel_idx_paths(child_dir_path_2,
                                                         prepd_query_bords):
                    yield neces_idx_path

    def eq(self,
           *queries: Any) -> Generator:
        for query in queries:
            prepd_query_bords = self.prep_query(query)
            for neces_top_dir_name in self.sel_top_dir_names(prepd_query_bords):
                for neces_idx_path in self.sel_idx_paths(neces_top_dir_name,
                                                         prepd_query_bords):
                    neces_idx = self.read_idx(neces_idx_path)
                    neces_line_starts = self.read_idx(f'{neces_idx_path[:-4]}.b')
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
        for neces_top_dir_name in self.sel_top_dir_names(prepd_query_bords):
            for neces_idx_path in self.sel_idx_paths(neces_top_dir_name,
                                                     prepd_query_bords):
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
