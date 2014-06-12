'''
Created on Aug 8, 2013

@author: smirarab
'''
import sys
from pasta.alignment import Alignment

def compact_to_fasta(inf, out):
    a = Alignment()
    a.read_filepath(inf,"COMPACT3")
    a.write(out, 'FASTA')

if __name__ == '__main__':
    compact_to_fasta(sys.argv[1], '%s.FASTA' %sys.argv[1] if len(sys.argv) == 1 else sys.argv[2])
    