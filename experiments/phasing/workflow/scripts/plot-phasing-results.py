import argparse
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='plot-phasing-results.py', description=__doc__)
parser.add_argument('-tsvfiles', metavar='TSVFILES', nargs='+', help='List of TSV files from whatshap compare.')
parser.add_argument('-truthsetname', metavar='TRUTHSETNAME', help='Name of the truth set for figure title.' )
parser.add_argument('-outname', metavar='OUTNAME', help='Name of the output PDF.')
args = parser.parse_args()


samples = []
switch_error_rates = []

for tsv_file in args.tsvfiles:
	df = pd.read_csv(tsv_file, sep='\t')
	sample = df['#sample'][0]
	switch_errors = sum(df['all_switches'])
	variants = sum(df['all_assessed_pairs'])
	samples.append(sample)
	switch_error_rates.append((switch_errors / variants) * 100.0)


print(samples)
print(switch_error_rates)

# plot the results
plt.figure()

x_values = [x for x in range(1, len(samples) + 1)]

plt.plot(x_values, switch_error_rates, marker='o')
plt.xticks(x_values, samples, rotation=90)
plt.ylabel('Switch Error Rate [%]')
plt.title('Comparison against ' + args.truthsetname + ' phasing')
plt.tight_layout()
plt.savefig(args.outname)
