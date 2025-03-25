import re
from collections.abc import Iterable

if __name__ == 'main':
    __version__ = 'v5.0.2'
    __authors__ = [{'name': 'Platon Bykadorov',
                    'email': 'platon.work@gmail.com',
                    'years': '2023-2025'}]


class SrtRules():
    @staticmethod
    def natur(src_str_or_row: str | Iterable,
              dec_delimiter: str = '.',
              nums_first: bool = True) -> list:
        if type(src_str_or_row) is str:
            src_row = [src_str_or_row]
        elif isinstance(src_str_or_row,
                        Iterable):
            src_row = list(map(str,
                               src_str_or_row))
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

    @staticmethod
    def letts_nums(src_srt: str) -> list:
        letts = re.search(r'^[a-zA-Z]+',
                          src_srt).group()
        nums = int(re.search(f'(?<=^{letts})\\d+$',
                             src_srt).group())
        spl_row = [letts, nums]
        return spl_row
