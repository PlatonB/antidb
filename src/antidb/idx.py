# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import os
from typing import (Callable,
                    Any,
                    Generator)
from datetime import datetime
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

__version__ = 'v3.1.0'
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
                 idx_prefix: str,
                 your_line_prs: Callable,
                 your_line_prs_kwargs: None | dict = None,
                 compr_level: int = 6,
                 compr_frame_size: int = 1024 * 1024,
                 compr_chunk_size: int = 1024 * 1024 * 1024,
                 presrt_chunk_elems_quan: int = 10000000,
                 idx_chunk_elems_quan: int = 1000000,
                 srt_rule: None | Callable = None,
                 srt_rule_kwargs: None | dict = None,
                 srt_rule_cols_delimiter: str | None = '\t',
                 srt_rule_col_inds: None | int | list | tuple = None):
        self.db_file_path = os.path.normpath(db_file_path)
        if self.db_file_path.endswith('.zst'):
            self.db_zst_path = self.db_file_path[:]
            self.adb_path = f'{self.db_file_path[:-4]}.{idx_prefix}.adb'
        else:
            self.db_zst_path = self.db_file_path + '.zst'
            self.adb_path = f'{self.db_file_path}.{idx_prefix}.adb'
        self.temp_dir_path = os.path.dirname(self.db_file_path)
        self.your_line_prs = your_line_prs
        if your_line_prs_kwargs:
            self.your_line_prs_kwargs = your_line_prs_kwargs
        else:
            self.your_line_prs_kwargs = {}
        self.compr_settings = {CParameter.compressionLevel: compr_level}
        self.compr_frame_size = compr_frame_size
        self.compr_chunk_size = compr_chunk_size
        self.presrt_chunk_elems_quan = presrt_chunk_elems_quan
        self.idx_chunk_elems_quan = idx_chunk_elems_quan
        super().__init__(srt_rule_cols_delimiter,
                         srt_rule_col_inds)
        if srt_rule:
            self.srt_rule = srt_rule
        else:
            self.srt_rule = self.natur
        if srt_rule_kwargs:
            self.srt_rule_kwargs = srt_rule_kwargs
        else:
            self.srt_rule_kwargs = {}
        self.presrtd_idxs_opened = []
        self.perf = []

    def idx(self):
        if not os.path.exists(self.db_zst_path):
            self.perf.append(self.crt_db_zst())
        self.perf.append(self.presrt_idxs())
        if not os.path.exists(self.adb_path):
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

    @staticmethod
    def read_presrtd_idx(presrtd_idx_opened: TemporaryFile) -> Generator:
        for obj_ind in range(load(presrtd_idx_opened)):
            obj = load(presrtd_idx_opened)
            yield obj

    def crt_idx(self,
                srtd_chunk: list,
                srtd_line_starts: list,
                adb_opened_w: ZipFile) -> None:
        with ZstdFile(adb_opened_w.open(f'{srtd_chunk[0]}.idx',
                                        mode='w'),
                      mode='w',
                      level_or_option=self.compr_settings) as idx_opened:
            dump(srtd_chunk,
                 idx_opened,
                 HIGHEST_PROTOCOL)
        with ZstdFile(adb_opened_w.open(f'{srtd_chunk[0]}.b',
                                        mode='w'),
                      mode='w',
                      level_or_option=self.compr_settings) as b_opened:
            dump(srtd_line_starts,
                 b_opened,
                 HIGHEST_PROTOCOL)

    def crt_idxs(self,
                 adb_opened_w: ZipFile) -> None:
        srtd_chunk, srtd_line_starts = [], []
        for obj in merge(*map(self.read_presrtd_idx,
                              self.presrtd_idxs_opened)):
            srtd_chunk.append(obj[0])
            srtd_line_starts.append(obj[1])
            if len(srtd_chunk) == self.idx_chunk_elems_quan:
                self.crt_idx(srtd_chunk,
                             srtd_line_starts,
                             adb_opened_w)
                srtd_chunk.clear()
                srtd_line_starts.clear()
        if srtd_chunk:
            self.crt_idx(srtd_chunk,
                         srtd_line_starts,
                         adb_opened_w)

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
                your_line_prs_out = self.your_line_prs(db_zst_line,
                                                       **self.your_line_prs_kwargs)
                if not your_line_prs_out:
                    continue
                elif type(your_line_prs_out) is tuple:
                    for your_line_prs_out_elem in your_line_prs_out:
                        chunk.append(self.srt_rule(your_line_prs_out_elem,
                                                   **self.srt_rule_kwargs))
                        line_starts.append(db_zst_lstart)
                else:
                    chunk.append(self.srt_rule(your_line_prs_out,
                                               **self.srt_rule_kwargs))
                    line_starts.append(db_zst_lstart)
                if len(chunk) == self.presrt_chunk_elems_quan:
                    self.presrt_idx(chunk,
                                    line_starts)
                    chunk.clear()
                    line_starts.clear()

    @count_exec_time
    def crt_adb(self) -> None:
        with ZipFile(self.adb_path,
                     mode='w') as adb_opened_w:
            self.crt_idxs(adb_opened_w)
            with TextIOWrapper(adb_opened_w.open('meta.txt',
                                                 mode='w')) as meta_opened:
                meta_opened.write(f'your_line_prs_name={self.your_line_prs.__name__}\n')
                meta_opened.write(f'your_line_prs_kwargs={self.your_line_prs_kwargs}\n')
                meta_opened.write(f'srt_rule_name={self.srt_rule.__name__}\n')
                meta_opened.write(f'srt_rule_cols_delimiter={self.cols_delimiter}\n')
                meta_opened.write(f'srt_rule_col_inds={self.col_inds}\n')
                meta_opened.write(f'srt_rule_kwargs={self.srt_rule_kwargs}\n')
