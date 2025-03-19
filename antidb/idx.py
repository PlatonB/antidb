# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import os
from typing import (Callable,
                    Any,
                    Generator)
from datetime import datetime
from copy import deepcopy
from zipfile import ZipFile
from tempfile import TemporaryFile
from pickle import (dump,
                    load,
                    HIGHEST_PROTOCOL)
from heapq import merge
from io import TextIOWrapper
from .srt import SrtRules
from pyzstd import (CParameter,
                    SeekableZstdFile,
                    ZstdFile)

if __name__ == 'main':
    __version__ = 'v5.0.2'
    __authors__ = [{'name': 'Platon Bykadorov',
                    'email': 'platon.work@gmail.com',
                    'years': '2023-2025'}]


def count_exec_time(any_func: Callable) -> Callable:
    def wrapper(*args: Any, **kwargs: Any):
        exec_time_start = datetime.now()
        any_func_res = any_func(*args, **kwargs)
        return (any_func.__name__,
                any_func_res,
                str(datetime.now() -
                    exec_time_start))
    return wrapper


class Idx(SrtRules):
    def __init__(self,
                 db_file_path: str,
                 idx_name_prefix: str,
                 db_line_prs: Callable,
                 idx_srt_rule: Callable,
                 db_line_prs_kwargs: None | dict = None,
                 idx_srt_rule_kwargs: None | dict = None,
                 compr_level: int = 3,
                 compr_frame_size: int = 1024 * 1024,
                 compr_chunk_size: int = 1024 * 1024 * 1024,
                 presrt_chunk_elems_quan: int = 10000000,
                 idx_chunk_div: int = 100,
                 idx_max_chunk_elems_quan: int = 100000):
        super().__init__()
        self.db_file_path = os.path.normpath(db_file_path)
        if self.db_file_path.endswith('.zst'):
            self.db_zst_path = deepcopy(self.db_file_path)
            self.adb_path = f'{self.db_file_path[:-4]}.{idx_name_prefix}.adb'
        else:
            self.db_zst_path = self.db_file_path + '.zst'
            self.adb_path = f'{self.db_file_path}.{idx_name_prefix}.adb'
        self.temp_dir_path = os.path.dirname(self.db_file_path)
        self.db_line_prs = db_line_prs
        self.idx_srt_rule = idx_srt_rule
        if db_line_prs_kwargs:
            self.db_line_prs_kwargs = db_line_prs_kwargs
        else:
            self.db_line_prs_kwargs = {}
        if idx_srt_rule_kwargs:
            self.idx_srt_rule_kwargs = idx_srt_rule_kwargs
        else:
            self.idx_srt_rule_kwargs = {}
        self.presrtd_idxs_opened = []
        self.compr_settings = {CParameter.compressionLevel:
                               compr_level}
        self.compr_frame_size = compr_frame_size
        self.compr_chunk_size = compr_chunk_size
        self.presrt_chunk_elems_quan = presrt_chunk_elems_quan
        self.idx_chunk_div = idx_chunk_div
        if self.idx_chunk_div < 2:
            self.idx_chunk_div = 2
        self.idx_max_chunk_elems_quan = idx_max_chunk_elems_quan
        self.perf = []

    def idx(self):
        if not os.path.exists(self.db_zst_path):
            self.perf.append(self.crt_db_zst())
        if not os.path.exists(self.adb_path):
            self.perf.append(self.presrt_idxs())
            self.perf.append(self.crt_adb())
        for presrtd_idx_opened in self.presrtd_idxs_opened:
            presrtd_idx_opened.close()

    @count_exec_time
    def crt_db_zst(self) -> None:
        with open(self.db_file_path) as db_file_opened:
            with TextIOWrapper(SeekableZstdFile(self.db_zst_path,
                                                mode='w',
                                                level_or_option=self.compr_settings,
                                                max_frame_content_size=self.compr_frame_size)) as db_zst_opened:
                while True:
                    db_file_chunk = db_file_opened.read(self.compr_chunk_size)
                    if not db_file_chunk:
                        break
                    db_zst_opened.write(db_file_chunk)

    def presrt_idx(self,
                   chunk: list,
                   line_starts: list) -> None:
        presrtd_data = sorted(zip(chunk,
                                  line_starts))
        presrtd_idx_opened = TemporaryFile(dir=self.temp_dir_path)
        self.presrtd_idxs_opened.append(presrtd_idx_opened)
        dump(len(presrtd_data),
             presrtd_idx_opened)
        for obj in presrtd_data:
            dump(obj,
                 presrtd_idx_opened,
                 HIGHEST_PROTOCOL)
        presrtd_idx_opened.seek(0)

    @count_exec_time
    def presrt_idxs(self) -> None:
        with TextIOWrapper(SeekableZstdFile(self.db_zst_path)) as db_zst_opened:
            while True:
                db_zst_lstart = db_zst_opened.tell()
                if not db_zst_opened.readline().startswith('#'):
                    db_zst_opened.seek(db_zst_lstart)
                    break
            self.presrtd_idxs_opened.clear()
            chunk, line_starts = [], []
            while True:
                db_zst_lstart = db_zst_opened.tell()
                db_zst_line = db_zst_opened.readline().rstrip()
                if not db_zst_line:
                    if chunk:
                        self.presrt_idx(chunk,
                                        line_starts)
                    break
                db_line_prs_out = self.db_line_prs(db_zst_line,
                                                   **self.db_line_prs_kwargs)
                if not db_line_prs_out:
                    continue
                elif type(db_line_prs_out) is tuple:
                    for db_line_prs_out_elem in db_line_prs_out:
                        chunk.append(self.idx_srt_rule(db_line_prs_out_elem,
                                                       **self.idx_srt_rule_kwargs))
                        line_starts.append(db_zst_lstart)
                else:
                    chunk.append(self.idx_srt_rule(db_line_prs_out,
                                                   **self.idx_srt_rule_kwargs))
                    line_starts.append(db_zst_lstart)
                if len(chunk) == self.presrt_chunk_elems_quan:
                    self.presrt_idx(chunk,
                                    line_starts)
                    chunk.clear()
                    line_starts.clear()

    @staticmethod
    def read_presrtd_idx(presrtd_idx_opened: TemporaryFile) -> Generator:
        for obj_ind in range(load(presrtd_idx_opened)):
            obj = load(presrtd_idx_opened)
            yield obj

    def crt_idx(self,
                srtd_chunk: list,
                srtd_lstarts: list,
                low_dir_path: str,
                adb_opened_w: ZipFile) -> None:
        fir_chunk_elem = srtd_chunk[0]
        if type(fir_chunk_elem) is str:
            fir_chunk_elem = f"'{fir_chunk_elem}'"
        idx_path = os.path.join(low_dir_path,
                                f'{fir_chunk_elem}.idx')
        b_path = os.path.join(low_dir_path,
                              f'{fir_chunk_elem}.b')
        with ZstdFile(adb_opened_w.open(idx_path,
                                        mode='w'),
                      mode='w',
                      level_or_option=self.compr_settings) as idx_opened:
            dump(srtd_chunk,
                 idx_opened,
                 HIGHEST_PROTOCOL)
        with ZstdFile(adb_opened_w.open(b_path,
                                        mode='w'),
                      mode='w',
                      level_or_option=self.compr_settings) as b_opened:
            dump(srtd_lstarts,
                 b_opened,
                 HIGHEST_PROTOCOL)

    @staticmethod
    def crt_chi_dir(srtd_chunk_start: Any,
                    prev_srtd_chunk_start: Any,
                    name_dupl_num: int,
                    adb_opened_w: ZipFile,
                    cur_dir_path: str = '') -> tuple[int,
                                                     str]:
        if srtd_chunk_start != prev_srtd_chunk_start:
            name_dupl_num = 1
        else:
            name_dupl_num += 1
        if type(srtd_chunk_start) is str:
            srtd_chunk_start = f"'{srtd_chunk_start}'"
        chi_dir_name = f'{srtd_chunk_start}.{name_dupl_num}'
        chi_dir_path = os.path.join(cur_dir_path,
                                    chi_dir_name)
        adb_opened_w.mkdir(chi_dir_path)
        return name_dupl_num, chi_dir_path

    def crt_dir_tree(self,
                     cur_dir_path: str,
                     cur_srtd_chunk: list,
                     srtd_lstarts: list,
                     adb_opened_w: ZipFile,
                     min_srtd_chunk_flag: bool = False) -> None:
        cur_srtd_chunk_len = len(cur_srtd_chunk)
        if (cur_srtd_chunk_len <= self.idx_max_chunk_elems_quan
                or min_srtd_chunk_flag):
            self.crt_idx(cur_srtd_chunk,
                         srtd_lstarts,
                         cur_dir_path,
                         adb_opened_w)
            return None
        chi_srtd_chunk_len = cur_srtd_chunk_len // self.idx_chunk_div
        if chi_srtd_chunk_len < self.idx_max_chunk_elems_quan:
            idx_chunk_div = cur_srtd_chunk_len // self.idx_max_chunk_elems_quan
            if idx_chunk_div > 1:
                chi_srtd_chunk_len = cur_srtd_chunk_len // idx_chunk_div
            else:
                chi_srtd_chunk_len = cur_srtd_chunk_len // 2
            min_srtd_chunk_flag = True
        bord_inds = list(range(0, cur_srtd_chunk_len,
                               chi_srtd_chunk_len))
        chi_srtd_chunks = [cur_srtd_chunk[bord_ind:
                                          (bord_ind +
                                           chi_srtd_chunk_len)]
                           for bord_ind in bord_inds]
        chi_srtd_lstarts = [srtd_lstarts[bord_ind:
                                         (bord_ind +
                                          chi_srtd_chunk_len)]
                            for bord_ind in bord_inds]
        prev_chi_srtd_chunk_start = None
        name_dupl_num = 1
        for chi_srtd_chunk_ind in range(len(chi_srtd_chunks)):
            chi_srtd_chunk_start = chi_srtd_chunks[chi_srtd_chunk_ind][0]
            name_dupl_num, chi_dir_path = self.crt_chi_dir(chi_srtd_chunk_start,
                                                           prev_chi_srtd_chunk_start,
                                                           name_dupl_num,
                                                           adb_opened_w,
                                                           cur_dir_path)
            self.crt_dir_tree(chi_dir_path,
                              chi_srtd_chunks[chi_srtd_chunk_ind],
                              chi_srtd_lstarts[chi_srtd_chunk_ind],
                              adb_opened_w,
                              min_srtd_chunk_flag)
            prev_chi_srtd_chunk_start = deepcopy(chi_srtd_chunk_start)

    @count_exec_time
    def crt_adb(self) -> None:
        with ZipFile(self.adb_path,
                     mode='w') as adb_opened_w:
            srtd_chunk, srtd_lstarts = [], []
            prev_srtd_chunk_start = None
            name_dupl_num = 1
            for obj in merge(*map(self.read_presrtd_idx,
                                  self.presrtd_idxs_opened)):
                srtd_chunk.append(obj[0])
                srtd_lstarts.append(obj[1])
                if len(srtd_chunk) == self.presrt_chunk_elems_quan:
                    name_dupl_num, chi_dir_name = self.crt_chi_dir(srtd_chunk[0],
                                                                   prev_srtd_chunk_start,
                                                                   name_dupl_num,
                                                                   adb_opened_w)
                    self.crt_dir_tree(chi_dir_name,
                                      srtd_chunk,
                                      srtd_lstarts,
                                      adb_opened_w)
                    prev_srtd_chunk_start = srtd_chunk[0]
                    srtd_chunk.clear()
                    srtd_lstarts.clear()
            if srtd_chunk:
                name_dupl_num, chi_dir_name = self.crt_chi_dir(srtd_chunk[0],
                                                               prev_srtd_chunk_start,
                                                               name_dupl_num,
                                                               adb_opened_w)
                self.crt_dir_tree(chi_dir_name,
                                  srtd_chunk,
                                  srtd_lstarts,
                                  adb_opened_w)
