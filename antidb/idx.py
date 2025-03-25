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
    __version__ = 'v6.1.0'
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
                 adb_name_prefix: str,
                 db_line_prs: Callable,
                 adb_srt_rule: Callable,
                 db_line_prs_kwargs: None | dict = None,
                 adb_srt_rule_kwargs: None | dict = None,
                 compr_level: int = 3,
                 compr_frame_size: int = 1024 * 1024,
                 compr_chunk_size: int = 1024 * 1024 * 1024,
                 presrt_chunk_len: int = 40000000,
                 lstarts_idx_div: int = 1000,
                 lstarts_idx_len: int = 40000):
        super().__init__()
        self.db_file_path = os.path.normpath(db_file_path)
        if self.db_file_path.endswith('.zst'):
            self.db_zst_path = deepcopy(self.db_file_path)
            self.adb_path = f'{self.db_file_path[:-4]}.{adb_name_prefix}.adb'
        else:
            self.db_zst_path = self.db_file_path + '.zst'
            self.adb_path = f'{self.db_file_path}.{adb_name_prefix}.adb'
        self.temp_dir_path = os.path.dirname(self.db_file_path)
        self.db_line_prs = db_line_prs
        self.adb_srt_rule = adb_srt_rule
        if db_line_prs_kwargs:
            self.db_line_prs_kwargs = db_line_prs_kwargs
        else:
            self.db_line_prs_kwargs = {}
        if adb_srt_rule_kwargs:
            self.adb_srt_rule_kwargs = adb_srt_rule_kwargs
        else:
            self.adb_srt_rule_kwargs = {}
        self.presrtd_idxs_opened = []
        self.compr_settings = {CParameter.compressionLevel:
                               compr_level}
        self.compr_frame_size = compr_frame_size
        self.compr_chunk_size = compr_chunk_size
        self.presrt_chunk_len = presrt_chunk_len
        self.lstarts_idx_div = lstarts_idx_div
        if self.lstarts_idx_div < 2:
            self.lstarts_idx_div = 2
        self.lstarts_idx_len = lstarts_idx_len
        self.perf = []

    def idx(self) -> None:
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
                   vals_n_lstarts: list) -> None:
        vals_n_lstarts.sort()
        presrtd_idx_opened = TemporaryFile(dir=self.temp_dir_path)
        self.presrtd_idxs_opened.append(presrtd_idx_opened)
        dump(len(vals_n_lstarts),
             presrtd_idx_opened)
        for val_n_start in vals_n_lstarts:
            dump(val_n_start,
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
            vals_n_lstarts = []
            while True:
                db_zst_lstart = db_zst_opened.tell()
                db_zst_line = db_zst_opened.readline().rstrip()
                if not db_zst_line:
                    if vals_n_lstarts:
                        self.presrt_idx(vals_n_lstarts)
                    break
                db_line_prs_out = self.db_line_prs(db_zst_line,
                                                   **self.db_line_prs_kwargs)
                if not db_line_prs_out:
                    continue
                elif type(db_line_prs_out) is tuple:
                    for db_line_prs_out_elem in db_line_prs_out:
                        vals_n_lstarts.append([self.adb_srt_rule(db_line_prs_out_elem,
                                                                 **self.adb_srt_rule_kwargs),
                                               db_zst_lstart])
                else:
                    vals_n_lstarts.append([self.adb_srt_rule(db_line_prs_out,
                                                             **self.adb_srt_rule_kwargs),
                                           db_zst_lstart])
                if len(vals_n_lstarts) == self.presrt_chunk_len:
                    self.presrt_idx(vals_n_lstarts)
                    vals_n_lstarts.clear()

    @staticmethod
    def read_presrtd_idx(presrtd_idx_opened: TemporaryFile) -> Generator:
        for obj_ind in range(load(presrtd_idx_opened)):
            obj = load(presrtd_idx_opened)
            yield obj

    def crt_lstarts_idx(self,
                        vals_n_lstarts: list,
                        low_dir_path: str,
                        adb_opened_w: ZipFile) -> str:
        lstarts_idx_path = os.path.join(low_dir_path,
                                        'lstarts')
        with ZstdFile(adb_opened_w.open(lstarts_idx_path,
                                        mode='w'),
                      mode='w',
                      level_or_option=self.compr_settings) as lstarts_idx_opened:
            dump(list(zip(*vals_n_lstarts)),
                 lstarts_idx_opened,
                 HIGHEST_PROTOCOL)
        return lstarts_idx_path

    def crt_paths_idx(self,
                      adb_opened_w: ZipFile,
                      paths_idx_obj: list,
                      dir_path: str = '') -> str:
        paths_idx_path = os.path.join(dir_path,
                                      'paths')
        with adb_opened_w.open(paths_idx_path,
                               mode='w') as paths_idx_opened:
            dump(paths_idx_obj,
                 paths_idx_opened,
                 HIGHEST_PROTOCOL)
        return paths_idx_path

    def crt_dir_tree(self,
                     cur_dir_path: str,
                     cur_vals_n_lstarts: list,
                     adb_opened_w: ZipFile,
                     min_vals_n_lstarts_flag: bool = False) -> str | None:
        cur_vals_n_lstarts_len = len(cur_vals_n_lstarts)
        if (cur_vals_n_lstarts_len <= self.lstarts_idx_len
                or min_vals_n_lstarts_flag):
            lstarts_idx_path = self.crt_lstarts_idx(cur_vals_n_lstarts,
                                                    cur_dir_path,
                                                    adb_opened_w)
            return lstarts_idx_path
        chi_vals_n_lstarts_len = cur_vals_n_lstarts_len // self.lstarts_idx_div
        if chi_vals_n_lstarts_len < self.lstarts_idx_len:
            lstarts_idx_div = cur_vals_n_lstarts_len // self.lstarts_idx_len
            if lstarts_idx_div > 1:
                chi_vals_n_lstarts_len = cur_vals_n_lstarts_len // lstarts_idx_div
            else:
                chi_vals_n_lstarts_len = cur_vals_n_lstarts_len // 2
            min_vals_n_lstarts_flag = True
        bord_inds = list(range(0, cur_vals_n_lstarts_len,
                               chi_vals_n_lstarts_len))
        chi_vals_n_lstarts = [cur_vals_n_lstarts[bord_ind:
                                                 (bord_ind +
                                                  chi_vals_n_lstarts_len)]
                              for bord_ind in bord_inds]
        chi_dir_num = 1
        paths_idx_obj = [[], []]
        for ind in range(len(chi_vals_n_lstarts)):
            chi_dir_path = os.path.join(cur_dir_path,
                                        str(chi_dir_num))
            adb_opened_w.mkdir(chi_dir_path)
            chi_dir_num += 1
            gchi_any_idx_path = self.crt_dir_tree(chi_dir_path,
                                                  chi_vals_n_lstarts[ind],
                                                  adb_opened_w,
                                                  min_vals_n_lstarts_flag)
            paths_idx_obj[0].append(chi_vals_n_lstarts[ind][0][0])
            paths_idx_obj[1].append(gchi_any_idx_path)
        paths_idx_path = self.crt_paths_idx(adb_opened_w,
                                            paths_idx_obj,
                                            cur_dir_path)
        return paths_idx_path

    @count_exec_time
    def crt_adb(self) -> None:
        with ZipFile(self.adb_path,
                     mode='w') as adb_opened_w:
            vals_n_lstarts = []
            chi_dir_num = 1
            paths_idx_obj = [[], []]
            for val_n_lstart in merge(*map(self.read_presrtd_idx,
                                           self.presrtd_idxs_opened)):
                vals_n_lstarts.append(val_n_lstart)
                if len(vals_n_lstarts) == self.presrt_chunk_len:
                    chi_dir_name = str(chi_dir_num)
                    adb_opened_w.mkdir(chi_dir_name)
                    chi_dir_num += 1
                    gchi_any_idx_path = self.crt_dir_tree(chi_dir_name,
                                                          vals_n_lstarts,
                                                          adb_opened_w)
                    paths_idx_obj[0].append(vals_n_lstarts[0][0])
                    paths_idx_obj[1].append(gchi_any_idx_path)
                    vals_n_lstarts.clear()
            if vals_n_lstarts:
                chi_dir_name = str(chi_dir_num)
                adb_opened_w.mkdir(chi_dir_name)
                gchi_any_idx_path = self.crt_dir_tree(chi_dir_name,
                                                      vals_n_lstarts,
                                                      adb_opened_w)
                paths_idx_obj[0].append(vals_n_lstarts[0][0])
                paths_idx_obj[1].append(gchi_any_idx_path)
            self.crt_paths_idx(adb_opened_w,
                               paths_idx_obj)
