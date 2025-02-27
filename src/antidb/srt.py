import re
from collections.abc import Iterable
from .err import DelimitersMatchError

__version__ = 'v4.0.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2023-2025'}]


class SrtRules():
    def __init__(self,
                 cols_delimiter: str | None = '\t',
                 col_inds: None | int | list | tuple = None):
        self.cols_delimiter = cols_delimiter
        self.col_inds = col_inds

    def get_cols(self,
                 src_line: str):
        if self.cols_delimiter:
            src_row = src_line.rstrip().split(self.cols_delimiter)
        else:
            src_row = [src_line.rstrip()]
        if type(self.col_inds) is int:
            src_row = [src_row[self.col_inds]]
        elif type(self.col_inds) in [list, tuple]:
            src_row = [src_row[col_ind]
                       for col_ind in self.col_inds]
        return src_row

    def natur(self,
              src_line_or_row: str | list | tuple,
              dec_delimiter: str = '.',
              nums_first: bool = True) -> list:
        if self.cols_delimiter == dec_delimiter:
            raise DelimitersMatchError(self.cols_delimiter,
                                       dec_delimiter)
        if type(src_line_or_row) is str:
            src_row = self.get_cols(src_line_or_row)
        elif isinstance(src_line_or_row,
                        Iterable):
            src_row = src_line_or_row[:]
        if dec_delimiter == '.':
            natur_split_cell = r'(-?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)'
        elif dec_delimiter == ',':
            natur_split_cell = r'(-?\d+(?:,\d*)?(?:[Ee][+-]?\d+)?)'
        spl_row = []
        for cell in src_row:
            subcells = list(filter(lambda subcell:
                                   subcell,
                                   re.split(natur_split_cell,
                                            cell)))
            for subcell_ind in range(len(subcells)):
                try:
                    subcells[subcell_ind] = int(subcells[subcell_ind])
                except ValueError:
                    try:
                        subcells[subcell_ind] = float(subcells[subcell_ind])
                    except ValueError:
                        if dec_delimiter == ',':
                            try:
                                subcells[subcell_ind] = float(subcells[subcell_ind].replace(',', '.'))
                            except ValueError:
                                pass
            if type(subcells[0]) is str:
                if nums_first:
                    subcells.insert(0, float('+inf'))
                else:
                    subcells.insert(0, float('-inf'))
            spl_row.append(subcells)
        return spl_row

    def letts_nums(self,
                   src_line_or_row: str | list | tuple) -> list:
        if type(src_line_or_row) is str:
            src_row = self.get_cols(src_line_or_row)
        else:
            src_row = src_line_or_row[:]
        spl_row = []
        for cell in src_row:
            letts = re.search(r'^[a-zA-Z]+',
                              cell).group()
            nums = int(re.search(f'(?<=^{letts})\\d+$',
                                 cell).group())
            spl_row.append([letts,
                            nums])
        return spl_row
