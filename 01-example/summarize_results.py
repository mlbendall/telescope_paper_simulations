#! /bin/bash
from subprocess import call, Popen, PIPE
from collections import defaultdict, Counter
import sys

### Coordinates distinguishing the LTR and internal regions
references_bed = '''HML2_chr3_3q12.3	0	969	HML2_chr3_3q12.3_LTR5Hs
HML2_chr3_3q12.3	969	8154	HML2_chr3_3q12.3_INT
HML2_chr3_3q12.3	8154	9123	HML2_chr3_3q12.3_LTR5Hs
HML2_chr1_1q22	0	968	HML2_chr1_1q22_LTR5Hs
HML2_chr1_1q22	968	8212	HML2_chr1_1q22_INT
HML2_chr1_1q22	8212	9180	HML2_chr1_1q22_LTR5Hs
HML2_chr1_1q23.3	0	969	HML2_chr1_1q23.3_LTR5
HML2_chr1_1q23.3	969	8263	HML2_chr1_1q23.3_INT
HML2_chr1_1q23.3	8263	9232	HML2_chr1_1q23.3_LTR5
HML2_chr6_6q25.1	0	1857	HML2_chr6_6q25.1_INT
HML2_chr6_6q25.1	1857	2832	HML2_chr6_6q25.1_LTR5B
'''
with open('references.bed','w') as outh:
    print >>outh, references_bed

samp = 'sample_01'

### Create a file with the coordinates for each read
lines = [l.strip('\n').strip('>') for l in open('%s/%s_1.fasta' % (samp,samp), 'rU') if l.startswith('>')]
with open('%s/%s.bed' % (samp,samp),'w') as outh:
    for l in lines:
        f1,f2,f3 = l.split(';')
        r,chrom = f1.split('/')
        s1,e1 = f2.split(':')[1].split('-')
        s2,e2 = f3.split(':')[1].split('-')
        print >>outh, '%s\t%s\t%s\t%s_1' % (chrom,s1,e1,r)
        print >>outh, '%s\t%s\t%s\t%s_2' % (chrom,s2,e2,r)

##########################################################################################

### Use bedtools to intersect the files

call('module load bedtools && bedtools intersect -a %s/%s.bed -b references.bed -wao > intersect_reads.bed' % (samp,samp), shell=True)

### Calculate the expected number of reads for each feature
##########################################################################################

lines = [l.strip('\n').split('\t') for l in open('intersect_reads.bed','rU')]

# Resolve reads overlapping more than one region
rec = defaultdict(list)
for l in lines:
    rec[l[3]].append(l)

goodlines = []
for k,v in rec.iteritems():
    if len(v) > 1:
        ilen = [int(l[-1]) for l in v]
        # Take the largest overlap
        goodlines.append(v[ilen.index(max(ilen))])
    else:
        goodlines.append(v[0])

c1 = Counter([l[7] for l in goodlines])
c2 = Counter([l[7].split('_')[-1] for l in goodlines])

##########################################################################################
# Make table 1
table = {}
lookup = {'INT':'HERVK-int', 'LTR5Hs':'LTR5_Hs', 'LTR5B': 'LTR5B', 'LTR5': 'LTR5'}
for n,e in c2.most_common():
    name = lookup[n]
    if name not in table:
        table[name] = {'expected':0, 'count':0, 'family':None}
    table[name]['expected'] = int(e)

lines = [l.strip('\n').split('\t') for l in open('%s/RepEnrich/rmsk_ALL.hg19_fraction_counts.txt' % samp,'rU')]
nonzero = [l for l in lines if int(l[-1])>0]
for l in nonzero:
  if l[0] not in table:
      table[l[0]] = {'expected':0, 'count':0, 'family':None}
  table[l[0]]['family'] = l[2]
  table[l[0]]['count'] = int(l[3])

other_row = {'expected':0, 'count':0, 'family':'-'}
allnames = table.keys()
for k in allnames:
    if table[k]['family'] != 'ERVK':
      temprow = table.pop(k)
      other_row['count'] += temprow['count']

table['Other'] = other_row

sortedtable = sorted(table.items(),key=lambda x: x[1]['count'], reverse=True)
tot_exp = sum([v[1]['expected'] for v in sortedtable])
tot_count = sum([v[1]['count'] for v in sortedtable])
if tot_exp - tot_count:
    sortedtable.append(('Unassigned',{'family':'-','count':tot_exp - tot_count, 'expected':0}))

cols = ['subfamily','family','count','prop','expected','prop']

print >>sys.stdout, '### Table 1. RepEnrich estimates of TE expression for simulated HERVK sample.'
print >>sys.stdout, '\t'.join(_.title() for _ in cols)
for row in sortedtable:
    print >>sys.stdout, '%s\t%s\t%d\t%.1f%%\t%d\t%.1f%%' % (row[0], row[1]['family'],
                                                        row[1]['count'], (float(row[1]['count'])/tot_count)*100 ,
                                                        row[1]['expected'], (float(row[1]['expected'])/tot_exp)*100 )

##########################################################################################


master = ['HML2_1q22', 'HML2_1q23.3', 'HML2_3q12.3', 'HML2_6q25.1',]
table2 = {}
for k in master:
  table2[k] = {'unique': 0, 'best': 0, 'frac': 0, 'expected':100}

lines = [l.strip('\n').split('\t') for l in open('%s/tsout/knownHERV-telescope_report.tsv' % samp,'rU')]
for l in lines[2:]:
    if int(l[7])>0 or float(l[9])>0:
        shortname = l[0].split('_')[0] + '_' + '_'.join(l[0].split('_')[2:])
        if shortname not in table2:
            table2[shortname] = {'unique': 0, 'best': 0, 'frac': 0, 'expected':0}
        table2[shortname]['unique'] = int(l[7])
        table2[shortname]['frac'] = float(l[9])


cmd = """module load samtools && samtools view -f 0x40 -F 0x100  sample_01/align_bt2.tagged.bam | sed 's/^.*\sXF:Z:\([^\s]*\).*$/\\1/' """
p1 = Popen(cmd, shell=True, stdout=PIPE)
out,err = p1.communicate()
besthits = out.strip('\n').split('\n')

for k,n in Counter(besthits).most_common():
    shortname = k.split('_')[0] + '_' + '_'.join(k.split('_')[2:])
    if shortname not in table2:
        table2[shortname] = {'unique': 0, 'best': 0, 'frac': 0, 'expected':0}
    table2[shortname]['best'] = n

other_row = {'unique': 0, 'best': 0, 'frac': 0, 'expected':0}
allnames = table2.keys()
for k in allnames:
    if k not in master:
      temprow = table2.pop(k)
      other_row['unique'] += temprow['unique']
      other_row['best'] += temprow['best']
      other_row['frac'] += temprow['frac']

sortedtable = [(k,table2[k]) for k in master]

if other_row['unique'] > 0 or other_row['best'] > 0 or other_row['frac'] > 0:
    sortedtable.append( ('Other', other_row) )

tot_exp = sum([v[1]['expected'] for v in sortedtable])
tot_unique = sum([v[1]['unique'] for v in sortedtable])
tot_best = sum([v[1]['best'] for v in sortedtable])
tot_frac = sum([v[1]['frac'] for v in sortedtable])
# if tot_exp - tot_count:
sortedtable.append(('Unassigned', {'unique': tot_exp - tot_unique , 'best':tot_exp - tot_best, 'frac': int(round(tot_exp - tot_frac))}))

cols = ['locus','unique','best','fractional']

print >>sys.stdout, '### Table 2. Simple estimates of HERV expression.'
print >>sys.stdout, '\t'.join(_.title() for _ in cols)
for row in sortedtable:
    print >>sys.stdout, '%s\t%d\t%d\t%.2f' % (row[0], row[1]['unique'],row[1]['best'], row[1]['frac'])


### Telescope table ######################################################################

lines = [l.strip('\n').split('\t') for l in open('%s/tsout/knownHERV-telescope_report.tsv' % samp,'rU')]
master = ['HML2_1q22','HML2_1q23.3','HML2_3q12.3','HML2_6q25.1',]
table3 = {}
for k in master:
  table3[k] = {'count':0, 'expected':100, 'pi': 0.0}

for l in lines[2:]:
    if float(l[1])>0:
        shortname = l[0].split('_')[0] + '_' + '_'.join(l[0].split('_')[2:])
        if shortname not in table3:
            table3[shortname] = {'count':0, 'expected':0, 'pi': 0.0}
        table3[shortname]['count'] = int(round(float(l[1])))
        table3[shortname]['pi'] = float(l[3])

other_row = {'expected':0, 'count':0, 'pi':0.0}
allnames = table3.keys()
for k in allnames:
    if k not in master:
      temprow = table3.pop(k)
      other_row['count'] += temprow['count']
      other_row['pi'] += temprow['pi']

sortedtable = [(k,table3[k]) for k in master]# sorted(table3.items(), key=lambda x: x[1]['count'], reverse=True)

if other_row['count'] > 0:
    sortedtable.append( ('Other', other_row) )

tot_exp = sum([v[1]['expected'] for v in sortedtable])
tot_count = sum([v[1]['count'] for v in sortedtable])
# if tot_exp - tot_count:
sortedtable.append(('Unassigned', {'pi':'-', 'count':tot_exp - tot_count, 'expected':0}))

cols = ['locus','count', 'pi', 'expected']

print >>sys.stdout, '### Table 3. Telescope estimates of TE expression for simulated HERVK sample.'
print >>sys.stdout, '\t'.join(_.title() for _ in cols)
for row in sortedtable:
    if type(row[1]['pi']) is str:
        print >>sys.stdout, '%s\t%s\t%s\t%d' % (row[0], row[1]['count'], row[1]['pi'], row[1]['expected'])
    else:
        print >>sys.stdout, '%s\t%s\t%.2f%%\t%d' % (row[0], row[1]['count'], row[1]['pi']*100, row[1]['expected'])



