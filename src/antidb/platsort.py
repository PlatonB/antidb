import re
import os
from multiprocessing import Pool
from heapq import merge
from functools import partial

__version__ = 'v0.4.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023'}]


def natur_sort_rule(src_file_line,
                    delimiter='\t',
                    src_file_colinds=None):
    src_file_row = src_file_line.rstrip().split(delimiter)
    if type(src_file_colinds) is int:
        src_file_row = [src_file_row[src_file_colinds]]
    elif type(src_file_colinds) in [list, tuple]:
        src_file_row = [src_file_row[src_file_colind]
                        for src_file_colind in src_file_colinds]
    spl_file_row = []
    for src_file_cell in src_file_row:
        src_file_subcells = list(filter(lambda src_file_subcell:
                                        src_file_subcell,
                                        re.split(r'(\d+)',
                                                 src_file_cell)))
        for src_file_subcell_ind in range(len(src_file_subcells)):
            try:
                src_file_subcells[src_file_subcell_ind] = int(src_file_subcells[src_file_subcell_ind])
            except ValueError:
                pass
        if type(src_file_subcells[0]) is str:
            src_file_subcells.insert(0, float('+inf'))
        spl_file_row.append(src_file_subcells)
    return spl_file_row


def iter_file(src_file_path):
    with open(src_file_path) as src_file_opened:
        for src_file_line in src_file_opened:
            yield src_file_line


def pre_srt(src_file_path,
            chunk_elems_quan=10000000,
            delimiter='\t',
            src_file_colinds=None):
    with open(src_file_path) as src_file_opened:
        while True:
            src_file_lstart = src_file_opened.tell()
            if not src_file_opened.readline().startswith('#'):
                src_file_opened.seek(src_file_lstart)
                break
        presrtd_file_paths = []
        chunk, chunk_len, chunk_num = [], 0, 0
        for src_file_line in src_file_opened:
            chunk.append(src_file_line)
            chunk_len += 1
            if chunk_len == chunk_elems_quan:
                chunk_num += 1
                presrtd_file_path = f'{src_file_path}.{chunk_num}'
                presrtd_file_paths.append(presrtd_file_path)
                chunk.sort(key=partial(natur_sort_rule,
                                       delimiter=delimiter,
                                       src_file_colinds=src_file_colinds))
                with open(presrtd_file_path, mode='w') as presrtd_file_opened:
                    for presrtd_file_line in chunk:
                        presrtd_file_opened.write(presrtd_file_line)
                chunk, chunk_len = [], 0
        if chunk:
            chunk_num += 1
            presrtd_file_path = f'{src_file_path}.{chunk_num}'
            presrtd_file_paths.append(presrtd_file_path)
            chunk.sort(key=partial(natur_sort_rule,
                                   delimiter=delimiter,
                                   src_file_colinds=src_file_colinds))
            with open(presrtd_file_path, mode='w') as presrtd_file_opened:
                for presrtd_file_line in chunk:
                    presrtd_file_opened.write(presrtd_file_line)
    return presrtd_file_paths


def mrg_srt(presrtd_file_paths,
            mrg_file_suff='srtd',
            delimiter='\t',
            src_file_colinds=None):
    if not presrtd_file_paths:
        return None
    presrtd_file_common_path = re.sub(r'\.\d+$',
                                      '',
                                      presrtd_file_paths[0])
    mrg_file_path = f'{presrtd_file_common_path}.{mrg_file_suff}'
    if len(presrtd_file_paths) == 1:
        os.rename(presrtd_file_paths[0],
                  mrg_file_path)
        return mrg_file_path
    with open(mrg_file_path, mode='w') as mrg_file_opened:
        for mrg_file_line in merge(*map(iter_file,
                                        presrtd_file_paths),
                                   key=partial(natur_sort_rule,
                                               delimiter=delimiter,
                                               src_file_colinds=src_file_colinds)):
            mrg_file_opened.write(mrg_file_line)
    return mrg_file_path


def srt(src_file_path,
        chunk_elems_quan=10000000,
        delimiter='\t',
        src_file_colinds=None,
        mrg_file_suff='srtd'):
    presrtd_file_paths = pre_srt(src_file_path,
                                 chunk_elems_quan,
                                 delimiter,
                                 src_file_colinds)
    mrg_file_path = mrg_srt(presrtd_file_paths,
                            mrg_file_suff,
                            delimiter,
                            src_file_colinds)
    for presrtd_file_path in presrtd_file_paths:
        os.remove(presrtd_file_path)
    return mrg_file_path
