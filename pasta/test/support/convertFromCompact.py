'''
Created on Aug 8, 2013

@author: smirarab
'''
import sys
from pasta.alignment import Alignment

def compact_to_fasta(inf):
    a = Alignment()
    a.read_filepath(inf,"COMPACT3")
    a.write('%s.FASTA' %inf, 'FASTA')

if __name__ == '__main__':
    compact_to_fasta(sys.argv[1])
    