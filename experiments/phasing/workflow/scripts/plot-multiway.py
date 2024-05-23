import sys
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(prog='plot-multiway.py', description=__doc__)
parser.add_argument('-tsv', metavar='TSV', nargs='+', required=True, help='Output of whatshap compare multiway comparison.')
parser.add_argument('-trios', metavar='TRIOS', required=True, help="ped file with trio information.")
parser.add_argument('-outname', metavar='OUTNAME', required=True, help="name of the output file.")
args = parser.parse_args()


sample_to_count = defaultdict(lambda: 0)
samples = set([])
sets = set([])

for tsvfile in args.tsv:
	for line in open(tsvfile, 'r'):
		fields = line.strip().split()
		if line.startswith('#'):
			continue
		sample_to_count[(fields[0], fields[2]+fields[3])] += int(fields[4])
		samples.add(fields[0])
		sets.add(fields[2]+fields[3])
trios = {}
for line in open(args.trios, 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	if fields[1] == 'NA' or fields[2] == 'NA':
		continue
	if (fields[0] in samples) and (fields[1] in samples) and (fields[2] in samples):
		trios[fields[0]] = [fields[1], fields[2]]
print('sets', sets)
print('samples', samples)
print('sample_to_count', sample_to_count)
print('trios', trios)

samples_to_plot = []
x_values = []


i = 0
n = 0
for k,v in trios.items():
	samples_to_plot.extend([k, v[0], v[1]])
	x_values.extend([i + n, i+1+n, i+2+n])
	n += 1
	i += 3

print(x_values)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05)  # adjust space between axes

# plot the same data on both axes
for s in sets:
	values = []
	for sample in samples_to_plot:
		print(sample)
		values.append(sample_to_count[(sample, s)])
	ax1.plot(x_values, values, label=s)
	ax2.plot(x_values, values, label=s)


# zoom-in / limit the view to different portions of the data
ax1.set_ylim(1500000, 3000000)  # outliers only
ax2.set_ylim(0, 60000)  # most of the data

ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

plt.xticks(x_values, samples_to_plot, rotation=90)
ax1.set_ylabel('Count')
ax2.set_ylabel('Count')
ax1.legend()
plt.tight_layout()
plt.savefig(args.outname)
