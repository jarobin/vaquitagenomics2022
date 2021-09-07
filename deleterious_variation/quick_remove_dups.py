### quick_remove_dups.py
# Clare Marsden 2015

## a quick method to remove dups.
## note a duplicate is based on the first two columns i.e. chrom and pos
## if the fold site is not the same, then this script will not catch that. it just will print out the first time it sees it.

import sys


print ' \n \n # # # # # # # # # # # # # # #'
print '\nBEWARE This script removes duplicates based on chrom and position columns only\n'

infilename=sys.argv[1]
outfilename='nodups_'+infilename

pos_seen = set() # holds lines already seen

outfile = open(outfilename, "w")
outfile.write('## Duplicates (based on chrom and pos) have been removed\n')
 

with open(infilename, 'rU') as infile:
	for line in open(infilename, "r"):
		cols=line.strip().split('\t')
		if line.startswith('#'):
			outfile.write(line)
		elif len(cols) != 3:
			print line
			sys.exit('## error exiting -  your line doesnt have three columns')
		else:
			pos=cols[0] + cols[1]
#			print 'pos', pos, 'poseen', pos_seen
			if pos not in pos_seen: # not a duplicate
#				print 'not dup'
				outfile.write(line)
				pos_seen.add(pos)
			else:
#				print 'dup'
				pass
outfile.close()

print '\n\n # # # # # DONE # # # # # # #\n'

