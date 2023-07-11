# autopep8: off
import sys; sys.dont_write_bytecode = True
# autopep8: on
import os
from datetime import datetime
from io import TextIOWrapper
from functools import wraps
from decimal import Decimal
from bisect import bisect
from pyzstd import (CParameter,
                    SeekableZstdFile,
                    ZstdFile)

__version__ = 'v1.3.3'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023'}]


def count_exec_time(any_func):
    def wrapper(*args):
        exec_time_start = datetime.now()
        any_func(*args)
        return (any_func.__name__,
                str(datetime.now() -
                    exec_time_start))
    return wrapper


class FileNotFoundError(Exception):
    def __init__(self, file_path):
        err_msg = f'\n{file_path} is missing'
        super().__init__(err_msg)


class Idx():
    def __init__(self,
                 db_file_path,
                 idx_prefix,
                 compr_level=6,
                 compr_frame_size=1024 * 1024,
                 compr_chunk_size=1024 * 1024 * 1024,
                 compr_chunk_elems_quan=1024 * 1024 * 1024 // 72,
                 unidx_lines_quan=1000):
        self.db_file_path = db_file_path
        if os.path.basename(db_file_path).endswith('.zst'):
            self.db_zst_path = self.db_file_path[:]
        else:
            self.db_zst_path = self.db_file_path + '.zst'
        self.idx_prefix = idx_prefix
        self.full_idx_path = f'{self.db_zst_path}.{self.idx_prefix}.full'
        self.full_idx_tmp_path = self.full_idx_path + '.tmp'
        self.full_idx_tmp_srtd_path = self.full_idx_tmp_path + '.srtd'
        self.mem_idx_path = f'{self.db_zst_path}.{self.idx_prefix}.mem'
        self.compr_settings = {CParameter.compressionLevel: compr_level}
        self.compr_frame_size = compr_frame_size
        self.compr_chunk_size = compr_chunk_size
        self.compr_chunk_elems_quan = compr_chunk_elems_quan
        self.unidx_lines_quan = unidx_lines_quan
        self.perf = []

    def idx(self, your_line_parser):
        @wraps(your_line_parser)
        def mng():
            if not os.path.exists(self.db_zst_path):
                self.perf.append(self.crt_db_zst())
                os.remove(self.db_file_path)
            if not os.path.exists(self.full_idx_tmp_path) \
                    and not os.path.exists(self.full_idx_path):
                self.perf.append(self.crt_full_idx_tmp(your_line_parser))
            if not os.path.exists(self.full_idx_tmp_srtd_path) \
                    and not os.path.exists(self.full_idx_path):
                self.perf.append(self.crt_full_idx_tmp_srtd())
                os.remove(self.full_idx_tmp_path)
            if not os.path.exists(self.full_idx_path):
                self.perf.append(self.crt_full_idx())
                os.remove(self.full_idx_tmp_srtd_path)
            if not os.path.exists(self.mem_idx_path):
                self.perf.append(self.crt_mem_idx())
        return mng

    @count_exec_time
    def crt_db_zst(self):
        with open(self.db_file_path) as src_txt_opened:
            with TextIOWrapper(SeekableZstdFile(self.db_zst_path,
                                                mode='w',
                                                level_or_option=self.compr_settings,
                                                max_frame_content_size=self.compr_frame_size)) as db_zst_opened:
                while True:
                    src_txt_chunk = src_txt_opened.read(self.compr_chunk_size)
                    if not src_txt_chunk:
                        break
                    db_zst_opened.write(src_txt_chunk)

    @count_exec_time
    def crt_full_idx_tmp(self, your_line_parser):
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
                    your_line_parser_out = your_line_parser(db_zst_line)
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
        os.system(f"sort -t ',' {self.full_idx_tmp_path} > {self.full_idx_tmp_srtd_path}")

    @count_exec_time
    def crt_full_idx(self):
        with open(self.full_idx_tmp_srtd_path) as full_idx_tmp_srtd_opened:
            with TextIOWrapper(SeekableZstdFile(self.full_idx_path,
                                                mode='w',
                                                level_or_option=self.compr_settings,
                                                max_frame_content_size=self.compr_frame_size)) as full_idx_opened:
                while True:
                    full_idx_tmp_srtd_chunk = full_idx_tmp_srtd_opened.read(self.compr_chunk_size)
                    if not full_idx_tmp_srtd_chunk:
                        break
                    full_idx_opened.write(full_idx_tmp_srtd_chunk)

    @count_exec_time
    def crt_mem_idx(self):
        with TextIOWrapper(SeekableZstdFile(self.full_idx_path,
                                            mode='r')) as full_idx_opened:
            with TextIOWrapper(ZstdFile(self.mem_idx_path,
                                        mode='w',
                                        level_or_option=self.compr_settings)) as mem_idx_opened:
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
    def __init__(self, db_file_path, idx_prefix):
        super().__init__(db_file_path, idx_prefix)
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
        self.mem_idx_your_vals, self.full_idx_lstarts = self.read_mem_idx()

    def read_mem_idx(self):
        mem_idx_your_vals, full_idx_lstarts = [], []
        for mem_idx_line in self.mem_idx_opened:
            mem_idx_row = mem_idx_line.rstrip().split(',')
            mem_idx_your_vals.append(mem_idx_row[0])
            full_idx_lstarts.append(int(mem_idx_row[1]))
        return mem_idx_your_vals, full_idx_lstarts

    def prs(self, your_vals):
        if type(your_vals) in [str,
                               int,
                               float,
                               Decimal]:
            your_vals = [your_vals]
        for your_val in your_vals:
            your_val = str(your_val)
            mem_idx_left_val_ind = bisect(self.mem_idx_your_vals,
                                          your_val) - 1
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
