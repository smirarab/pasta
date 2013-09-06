'''
Created on Aug 8, 2013

@author: smirarab
'''
import sys
from sate.alignment import Alignment

def compact_to_fasta(inf):
    a = Alignment()
    a.read_filepath(inf,"COMPACT3")
    a.write('%s.compatc3' %inf, 'COMPACT3')

if __name__ == '__main__':
    compact_to_fasta(sys.argv[1])
    