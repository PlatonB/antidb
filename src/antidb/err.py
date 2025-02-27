__version__ = 'v1.0.0'
__authors__ = [{'name': 'Platon Bykadorov',
                'email': 'platon.work@gmail.com',
                'years': '2025'}]


class DelimitersMatchError(Exception):
    def __init__(self,
                 cols_delimiter,
                 dec_delimiter):
        err_msg = f'''\nColums delimiter ({cols_delimiter})
matches decimal delimiter ({dec_delimiter})'''
        super().__init__(err_msg)


class NoIdxsError(Exception):
    def __init__(self):
        err_msg = f'\nThere are no indexes in DB'
        super().__init__(err_msg)


class QueryStartGtEndError(Exception):
    def __init__(self,
                 query_start,
                 query_end):
        err_msg = f'''\nQuery start ({query_start})
more then query end ({query_end})'''
        super().__init__(err_msg)
