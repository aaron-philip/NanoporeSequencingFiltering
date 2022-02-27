#writing a .rawsig file for cwDTW_nano use

import math
from typing import Optional
import os

path_to_fast5 = '/Users/aaronphilip/ScienceFair/projects/NanoporeSequencingFiltering/DeepSimulator/artificial_human_chr22_DeepSimu/fast5'
fname = 'signal_0_30d8bd50-3aae-4eab-95da-2083e9df0631'
os.system('h5ls -d %s/%s.fast5/Raw/Reads/Read_981/Signal >  %s/%s.data' % (path_to_fast5, fname, path_to_fast5, fname))
raw = []

#temporary data file to write to
with open('%s/%s.data' % (path_to_fast5, fname), 'r') as data:
    data.readline()
    data.readline()
    line = data.readline()
    while line != '':
        splitarr = line.split()
        splitarr.pop(0)
        last = splitarr[-1]
        for current in range(len(splitarr)):
            splitarr[current] = splitarr[current][:-1]
        raw.extend(splitarr)
        line = data.readline()
    raw[-1] = last

with open('%s/%s.rawsig' % (path_to_fast5, fname), 'w') as rawsig:
    for term in range(len(raw)):
        rawsig.write(raw[term])
        rawsig.write('\n')

os.system('rm %s/%s.data' % (path_to_fast5, fname))
