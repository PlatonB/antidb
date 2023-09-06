import re
import os
from multiprocessing import Pool
from functools import partial
from contextlib import ExitStack
from collections import deque

__version__ = 'v0.2.0'
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


def pre_sort(src_file_path,
             chunk_elems_quan=10000000,
             delimiter='\t',
             src_file_colinds=None):
    with open(src_file_path) as src_file_opened:
        while True:
            src_file_lstart = src_file_opened.tell()
            if not src_file_opened.readline().startswith('#'):
                src_file_opened.seek(src_file_lstart)
                break
        chunk, chunk_len, chunk_num = [], 0, 0
        for src_file_line in src_file_opened:
            chunk.append(src_file_line)
            chunk_len += 1
            if chunk_len == chunk_elems_quan:
                chunk_num += 1
                presrtd_file_path = f'{src_file_path}.{chunk_num}'
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
            chunk.sort(key=partial(natur_sort_rule,
                                   delimiter=delimiter,
                                   src_file_colinds=src_file_colinds))
            with open(presrtd_file_path, mode='w') as presrtd_file_opened:
                for presrtd_file_line in chunk:
                    presrtd_file_opened.write(presrtd_file_line)


def mrg_sort(presrtd_file_paths,
             mrg_file_suff,
             chunk_elems_quan=10000000,
             delimiter='\t',
             src_file_colinds=None):
    if not presrtd_file_paths:
        return None
    mrg_file_path = f'{presrtd_file_paths[0]}.{mrg_file_suff}'
    if len(presrtd_file_paths) == 1:
        os.rename(presrtd_file_paths[0],
                  mrg_file_path)
        return mrg_file_path
    chunk_elems_quan //= len(presrtd_file_paths)
    with ExitStack() as stack:
        presrtd_files_opened = [stack.enter_context(open(presrtd_file_path))
                                for presrtd_file_path in presrtd_file_paths]
        presrtd_files_opened_quan = len(presrtd_files_opened)
        chunks = []
        for presrtd_file_opened in presrtd_files_opened:
            chunks.append(deque([]))
            for chunk_vert_ind in range(chunk_elems_quan):
                presrtd_file_line = presrtd_file_opened.readline()
                if presrtd_file_line:
                    chunks[-1].append(presrtd_file_line)
                else:
                    break
        with open(mrg_file_path, mode='w') as mrg_file_opened:
            while True:
                if presrtd_files_opened_quan == 1:
                    for chunk_line in chunks[0]:
                        mrg_file_opened.write(chunk_line)
                    for presrtd_file_line in presrtd_files_opened[0]:
                        mrg_file_opened.write(presrtd_file_line)
                    break
                chunk_fir_lines = [chunks[chunk_horiz_ind][0]
                                   for chunk_horiz_ind in range(presrtd_files_opened_quan)]
                min_chunk_fir_line = min(chunk_fir_lines,
                                         key=partial(natur_sort_rule,
                                                     delimiter=delimiter,
                                                     src_file_colinds=src_file_colinds))
                for chunk_horiz_ind in range(presrtd_files_opened_quan):
                    if chunk_fir_lines[chunk_horiz_ind] == min_chunk_fir_line:
                        mrg_file_opened.write(chunks[chunk_horiz_ind].popleft())
                        if not chunks[chunk_horiz_ind]:
                            for chunk_vert_ind in range(chunk_elems_quan):
                                presrtd_file_line = presrtd_files_opened[chunk_horiz_ind].readline()
                                if presrtd_file_line:
                                    chunks[chunk_horiz_ind].append(presrtd_file_line)
                                else:
                                    break
                            if not chunks[chunk_horiz_ind]:
                                del presrtd_files_opened[chunk_horiz_ind]
                                presrtd_files_opened_quan -= 1
                                del chunks[chunk_horiz_ind]
                                break
    return mrg_file_path


def pair_sort(src_file_path,
              max_proc_quan=4,
              delimiter='\t',
              src_file_colinds=None):
    src_dir_path, src_file_name = os.path.split(src_file_path)
    presrtd_file_paths = [os.path.join(src_dir_path, any_file_name)
                          for any_file_name in os.listdir(src_dir_path)
                          if re.search(f'{src_file_name}.\d+', any_file_name)]
    presrtd_files_quan = len(presrtd_file_paths)
    if presrtd_files_quan == 1:
        mrg_file_path = mrg_sort(presrtd_file_paths,
                                 'srtd',
                                 delimiter,
                                 src_file_colinds)
        os.rename(mrg_file_path,
                  re.sub(r'(?:\.\d+)+(?=\.srtd)',
                         '', mrg_file_path))
        return None
    elif presrtd_files_quan % 2 == 1:
        del presrtd_file_paths[-1]
        presrtd_files_quan -= 1
    presrtd_file_pairs = [presrtd_file_paths[presrtd_file_ind:presrtd_file_ind + 2]
                          for presrtd_file_ind in range(0, presrtd_files_quan, 2)]
    presrtd_file_pairs_quan = len(presrtd_file_pairs)
    presrtd_files_pair_inds = list(range(1, presrtd_file_pairs_quan + 1))
    proc_quan = min(max_proc_quan,
                    presrtd_file_pairs_quan,
                    os.cpu_count() // 2)
    with Pool(proc_quan) as pool_obj:
        pool_obj.starmap(partial(mrg_sort,
                                 delimiter=delimiter,
                                 src_file_colinds=src_file_colinds),
                         zip(presrtd_file_pairs,
                             presrtd_files_pair_inds))
    for presrtd_file_path in presrtd_file_paths:
        os.remove(presrtd_file_path)
    pair_sort(src_file_path,
              max_proc_quan,
              delimiter,
              src_file_colinds)


def plat_sort(src_file_path,
              chunk_elems_quan=10000000,
              max_proc_quan=8,
              delimiter='\t',
              src_file_colinds=0):
    pre_sort(src_file_path,
             chunk_elems_quan,
             delimiter,
             src_file_colinds)
    pair_sort(src_file_path,
              max_proc_quan,
              delimiter,
              src_file_colinds)
