from idx import (Idx,
                 count_exec_time)
from prs import Prs

__version__ = 'v1.2.1'

dbsnp_vcf_path = '/path/to/GCF_000001405.40[.zst]'
dbsnp_idx_prefix = 'all_rsids'
dbsnp_idx = Idx(dbsnp_vcf_path,
                dbsnp_idx_prefix,
                lambda dbsnp_zst_line:
                dbsnp_zst_line.split('\t')[2],
                lambda rsid: rsid)
dbsnp_idx.idx()
dbsnp_prs = Prs(dbsnp_vcf_path,
                dbsnp_idx_prefix,
                lambda rsid: rsid)


@count_exec_time
def get_rsid_lines(dbsnp_prs: Prs):
    for dbsnp_vcfzst_line in dbsnp_prs.eq('rs1009150',
                                          'rs12044852',
                                          'rs4902496'):
        print(dbsnp_vcfzst_line)


print(get_rsid_lines(dbsnp_prs))
