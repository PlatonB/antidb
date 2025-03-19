if __name__ == 'main':
    __version__ = 'v1.2.0'
    __authors__ = [{'name': 'Platon Bykadorov',
                    'email': 'platon.work@gmail.com',
                    'years': '2025'}]


class QueryStartGtEndError(Exception):
    def __init__(self,
                 query_start,
                 query_end):
        err_msg = f'''\nQuery start ({query_start})
more then query end ({query_end})'''
        super().__init__(err_msg)
