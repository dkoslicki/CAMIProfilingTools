# This is an example of how to use the script
import ProfilingTools as PF
import numpy as np
import argparse
import sys
import os

__author__ = 'David Koslicki (dmkoslicki@gmail.com, david.koslicki@math.oregonstate.edu)'
__version__ = '1.0.0'
__date__ = '29 Mar 2017'


def read_params(args):
	parser = argparse.ArgumentParser(description='')
	arg = parser.add_argument
	arg('--input', metavar='files_file', type=str, required=True,
					default=None, help="File of CAMI profile files to use")
	arg('--output', metavar='output_file', required=True, default=None, type=str,
					help="Output file (you should have this end in .csv as it is a matrix)")
	arg('--threshold', metavar='threshold', required=True, default=None, type=float,
					help="Value to threshold profiles to before computing EMDUnifrac. "
					"NOTE THIS VALUE IS IN PERCENTAGES so if you want 1% use 1")
	return vars(parser.parse_args())

if __name__ == '__main__':
	par = read_params(sys.argv)
	files_file = par['input']
	output_file = par['output']
	threshold = par['threshold']

	# Get all the profile file names
	files = []
	fid = open(files_file, 'r')
	for line in fid.readlines():
		files.append(line.strip())
	fid.close()

	# Import all the profiles
	profiles = []
	for file_name in files:
		profiles.append(PF.Profile(input_file_name=file_name))

	# Threshold all the profiles
	for profile in profiles:
		profile.threshold(threshold=threshold)

	# Compute EMDUnifrac
	D = np.zeros((len(profiles), len(profiles)))
	for i in xrange(len(profiles)):
		for j in xrange(i+1, len(profiles)):
			val = profiles[i].unifrac(profiles[j])
			D[i, j] = val
			D[j, i] = val

	np.savetxt(output_file, D, delimiter=',', newline='\n')
