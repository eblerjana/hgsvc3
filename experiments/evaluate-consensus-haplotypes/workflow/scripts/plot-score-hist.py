import sys
import matplotlib.pyplot as plt
import statistics


phred_scores = []
ref_ht = None

for line in sys.stdin:
	if line.startswith('ref_HT'):
		continue
	fields = line.strip().split()
	if fields[16] != 'inf':
		phred_scores.append(float(fields[16]))
	ref_ht = '_'.join(fields[0].split('_')[:-1])


plt.figure()
plt.title(ref_ht)
plt.hist(phred_scores, bins = range(int(min(phred_scores)) - 1, int(max(phred_scores)) + 2, 1))
plt.axvline(statistics.median(phred_scores), color='k', linestyle='dashed', linewidth=1)
plt.xlabel('Phred scores')
plt.ylabel('Count')
plt.savefig(sys.argv[1])

print('min: ' + str(min(phred_scores)))
print('max: ' + str(max(phred_scores)))
print('mean: ' + str(statistics.mean(phred_scores)))
print('median: ' + str(statistics.median(phred_scores)))
