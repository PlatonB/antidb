# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import os
from typing import (Callable,
                    Any)
from datetime import datetime
from io import TextIOWrapper
from decimal import Decimal
from warnings import warn
from bisect import bisect
from functools import partial
from .antisrt import (SrtRules,
                      Srt)
from pyzstd import (CParameter,
                    SeekableZstdFile,
                    ZstdFile)

__version__ = 'v2.6.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023'}]


def count_exec_time(any_func: Callable) -> Callable:
    def wrapper(*args: Any, **kwargs: Any):
        exec_time_start = datetime.now()
        any_func_res = any_func(*args, **kwargs)
        return (any_func.__name__,
                any_func_res,
                str(datetime.now() -
                    exec_time_start))
    return wrapper


class FileNotFoundError(Exception):
    def __init__(self, file_path):
        err_msg = f'\n{file_path} is missing'
        super().__init__(err_msg)


class Idx(Srt):
    def __init__(self,
                 db_file_path: str,
                 idx_prefix: str,
                 your_line_parser: Callable,
                 compr_level: int = 6,
                 compr_frame_size: int = 1024 * 1024,
                 compr_chunk_size: int = 1024 * 1024 * 1024,
                 compr_chunk_elems_quan: int = 10000000,
                 unidx_lines_quan: int = 1000,
                 srt_rule: None | Callable = None,
                 **srt_rule_kwargs: Any):
        self.db_file_path = os.path.normpath(db_file_path)
        if os.path.basename(db_file_path).endswith('.zst'):
            self.db_zst_path = self.db_file_path[:]
        else:
            self.db_zst_path = self.db_file_path + '.zst'
        self.idx_prefix = idx_prefix
        self.your_line_parser = your_line_parser
        self.full_idx_path = f'{self.db_zst_path}.{self.idx_prefix}.full'
        self.full_idx_tmp_path = self.full_idx_path + '.tmp'
        self.full_idx_tmp_srtd_path = self.full_idx_tmp_path + '.srtd'
        self.mem_idx_path = f'{self.db_zst_path}.{self.idx_prefix}.mem'
        self.compr_settings = {CParameter.compressionLevel: compr_level}
        self.compr_frame_size = compr_frame_size
        self.compr_chunk_size = compr_chunk_size
        self.compr_chunk_elems_quan = compr_chunk_elems_quan
        self.unidx_lines_quan = unidx_lines_quan
        super().__init__(unsrtd_file_path=self.full_idx_tmp_path,
                         srt_rule=srt_rule,
                         **srt_rule_kwargs)
        self.perf = []

    def idx(self):
        if not os.path.exists(self.db_zst_path):
            self.perf.append(self.crt_db_zst())
            os.remove(self.db_file_path)
        if not os.path.exists(self.full_idx_tmp_path) \
                and not os.path.exists(self.full_idx_path):
            self.perf.append(self.crt_full_idx_tmp())
        if not os.path.exists(self.full_idx_tmp_srtd_path) \
                and not os.path.exists(self.full_idx_path):
            self.perf.append(self.crt_full_idx_tmp_srtd())
            os.remove(self.full_idx_tmp_path)
        if not os.path.exists(self.full_idx_path):
            self.perf.append(self.crt_full_idx())
            os.remove(self.full_idx_tmp_srtd_path)
        if not os.path.exists(self.mem_idx_path):
            self.perf.append(self.crt_mem_idx())

    @staticmethod
    def compr_text_file(src_text_path,
                        trg_zst_path,
                        compr_settings,
                        compr_frame_size,
                        compr_chunk_size):
        with open(src_text_path) as src_text_opened:
            with TextIOWrapper(SeekableZstdFile(trg_zst_path,
                                                mode='w',
                                                level_or_option=compr_settings,
                                                max_frame_content_size=compr_frame_size)) as trg_zst_opened:
                while True:
                    src_text_chunk = src_text_opened.read(compr_chunk_size)
                    if not src_text_chunk:
                        break
                    trg_zst_opened.write(src_text_chunk)

    @count_exec_time
    def crt_db_zst(self):
        self.compr_text_file(self.db_file_path,
                             self.db_zst_path,
                             self.compr_settings,
                             self.compr_frame_size,
                             self.compr_chunk_size)

    @count_exec_time
    def crt_full_idx_tmp(self):
        with TextIOWrapper(SeekableZstdFile(self.db_zst_path,
                                            mode='r')) as db_zst_opened:
            with open(self.full_idx_tmp_path,
                      mode='w') as full_idx_tmp_opened:
                while True:
                    db_zst_lstart = db_zst_opened.tell()
                    if not db_zst_opened.readline().startswith('#'):
                        db_zst_opened.seek(db_zst_lstart)
                        break
                chunk = []
                chunk_len = 0
                while True:
                    db_zst_lstart = db_zst_opened.tell()
                    db_zst_line = db_zst_opened.readline()
                    if not db_zst_line:
                        if chunk:
                            full_idx_tmp_opened.write('\n'.join(chunk) + '\n')
                        break
                    your_line_parser_out = self.your_line_parser(db_zst_line)
                    if not your_line_parser_out:
                        continue
                    elif type(your_line_parser_out) in [str, int, float, Decimal]:
                        chunk.append(f'{your_line_parser_out},{db_zst_lstart}')
                    elif type(your_line_parser_out) in [list, tuple, set]:
                        for your_val in your_line_parser_out:
                            chunk.append(f'{your_val},{db_zst_lstart}')
                    chunk_len += 1
                    if chunk_len == self.compr_chunk_elems_quan:
                        full_idx_tmp_opened.write('\n'.join(chunk) + '\n')
                        chunk = []
                        chunk_len = 0

    @count_exec_time
    def crt_full_idx_tmp_srtd(self):
        self.pre_srt(chunk_elems_quan=self.compr_chunk_elems_quan)
        self.mrg_srt()

    @count_exec_time
    def crt_full_idx(self):
        self.compr_text_file(self.full_idx_tmp_srtd_path,
                             self.full_idx_path,
                             self.compr_settings,
                             self.compr_frame_size,
                             self.compr_chunk_size)

    @count_exec_time
    def crt_mem_idx(self):
        with TextIOWrapper(SeekableZstdFile(self.full_idx_path,
                                            mode='r')) as full_idx_opened:
            with TextIOWrapper(ZstdFile(self.mem_idx_path,
                                        mode='w',
                                        level_or_option=self.compr_settings)) as mem_idx_opened:
                mem_idx_opened.write(f'idx_srt_rule_name={self.srt_rule.__name__}\n')
                mem_idx_opened.write(f'idx_srt_rule_kwargs={self.srt_rule_kwargs}\n')
                mem_idx_opened.write(f'unidx_lines_quan={self.unidx_lines_quan}\n')
                while True:
                    full_idx_lstart = full_idx_opened.tell()
                    full_idx_line = full_idx_opened.readline()
                    if not full_idx_line:
                        break
                    full_idx_your_val = full_idx_line.split(',')[0]
                    mem_idx_opened.write(f'{full_idx_your_val},{full_idx_lstart}\n')
                    for full_idx_line_ind in range(self.unidx_lines_quan):
                        if not full_idx_opened.readline():
                            break


class Prs(Idx):
    def __init__(self,
                 db_file_path: str,
                 idx_prefix: str,
                 srt_rule: None | Callable = None,
                 **srt_rule_kwargs: Any):
        super().__init__(db_file_path,
                         idx_prefix,
                         your_line_parser=None,
                         srt_rule=srt_rule,
                         **srt_rule_kwargs)
        if not os.path.exists(self.db_zst_path):
            raise FileNotFoundError(self.db_zst_path)
        else:
            self.db_zst_opened = TextIOWrapper(SeekableZstdFile(self.db_zst_path,
                                                                mode='r'))
        if not os.path.exists(self.full_idx_path):
            raise FileNotFoundError(self.full_idx_path)
        else:
            self.full_idx_opened = TextIOWrapper(SeekableZstdFile(self.full_idx_path,
                                                                  mode='r'))
        if not os.path.exists(self.mem_idx_path):
            raise FileNotFoundError(self.mem_idx_path)
        else:
            self.mem_idx_opened = TextIOWrapper(ZstdFile(self.mem_idx_path,
                                                         mode='r'))
        mem_idx = self.read_mem_idx()
        self.idx_srt_rule_name, self.idx_srt_rule_kwargs, self.unidx_lines_quan = mem_idx[:3]
        self.mem_idx_your_vals, self.full_idx_lstarts = mem_idx[3:]
        if self.idx_srt_rule_name != self.srt_rule.__name__:
            warn(f"""Your sort key name ({self.srt_rule.__name__}) doesn't
                 match the index sort key name ({self.idx_srt_rule_name})""")
        if self.idx_srt_rule_kwargs != self.srt_rule_kwargs:
            warn(f"""Your sort key arguments ({self.srt_rule_kwargs}) don't
                 match the index sort key arguments ({self.idx_srt_rule_kwargs})""")

    def read_mem_idx(self):
        idx_srt_rule_name = self.mem_idx_opened.readline().rstrip().split('=')[1]
        idx_srt_rule_kwargs = eval(self.mem_idx_opened.readline().rstrip().split('=')[1])
        unidx_lines_quan = int(self.mem_idx_opened.readline().rstrip().split('=')[1])
        mem_idx_your_vals, full_idx_lstarts = [], []
        for mem_idx_line in self.mem_idx_opened:
            mem_idx_row = mem_idx_line.rstrip().split(',')
            mem_idx_your_vals.append(mem_idx_row[0])
            full_idx_lstarts.append(int(mem_idx_row[1]))
        return (idx_srt_rule_name, idx_srt_rule_kwargs, unidx_lines_quan,
                mem_idx_your_vals, full_idx_lstarts)

    def prs(self,
            your_vals: list | tuple | set | dict | str | int | float | Decimal) -> str:
        if type(your_vals) in [str,
                               int,
                               float,
                               Decimal]:
            your_vals = [your_vals]
        for your_val in your_vals:
            your_val = str(your_val)
            mem_idx_left_val_ind = bisect(self.mem_idx_your_vals,
                                          self.srt_rule(your_val,
                                                        **self.srt_rule_kwargs),
                                          key=partial(self.srt_rule,
                                                      **self.srt_rule_kwargs)) - 1
            full_idx_lstart = self.full_idx_lstarts[mem_idx_left_val_ind]
            self.full_idx_opened.seek(full_idx_lstart)
            for line_idx in range(self.unidx_lines_quan + 1):
                full_idx_line = self.full_idx_opened.readline()
                if not full_idx_line:
                    break
                full_idx_row = full_idx_line.rstrip().split(',')
                if your_val == full_idx_row[0]:
                    self.db_zst_opened.seek(int(full_idx_row[1]))
                    yield self.db_zst_opened.readline()
                    for full_idx_line in self.full_idx_opened:
                        full_idx_row = full_idx_line.rstrip().split(',')
                        if your_val != full_idx_row[0]:
                            break
                        self.db_zst_opened.seek(int(full_idx_row[1]))
                        yield self.db_zst_opened.readline()
                    break
