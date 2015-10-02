#! /usr/bin/env python
import sys
import re
from collections import Counter

lines = [l.strip('\n').split('\t') for l in sys.stdin]
outh = sys.stdout

fams = []
print >>outh, 'locus\tfamily\tlocation'
for l in lines:
    if 7000 < int(l[4])-int(l[3]) < 10000:
        attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
        fam = attr['family']
        if fam == 'HML2' or fam == 'HML5':
            fam = 'HERVK'
        if fam == 'HERVFRD':
            fam = 'HERVW'
        fams.append(fam)
        print >>outh, '%s\t%s\t%s:%s-%s' % (attr['locus'],fam,l[0],l[3],l[4])

print >>sys.stderr, 'Number of annotations by family:'
print >>sys.stderr, '\n'.join('%s - %s' % t for t in Counter(fams).most_common())
