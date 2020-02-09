from argparse import ArgumentParser
import logging
import sys
import pysam
import csv
from bx.intervals.intersection import IntervalTree
import numpy as np
import matplotlib.pyplot as plt
import os
from math import log



# The maximum number of characters to take from the prefix of an input BAM
# file to use in the legend of the graph.
MAX_BAM_PREFIX = 7
L = 100
S = L
# default name of log file if none is specified on the command line
DEFAULT_LOG_FILE = "bamcover.log"


def parse_args():
    'Parse the command line arguments for the program.'
    parser = ArgumentParser(
        description="Generate coverage information for BAM files")
    parser.add_argument(
        '--coords', required=True, type=str,
        help='TSV coordinates file (1-based) for region of interest')
    parser.add_argument(
        'bams', nargs='+', type=str,
        help='bam files containing mapped reads')
    parser.add_argument(
	'--sample', required=True, type=str,
        help='Sample Name')
    parser.add_argument(
        '--log', metavar='FILE', type=str, default=DEFAULT_LOG_FILE,
        help='Log progress in FILENAME.')
    return parser.parse_args()


def slices(array, l=L, step=S):
    '''
    rolling slices
    '''

    for i in range(0, len(array), step):
        slice = array[i:i+l]
        yield slice

def get_coords(coords_filename):
    '''Read the input coordinates file in TSV format and check
    that the start and end coordinates are integers. Log any
    coordinates which are skipped.'''
    result = []
    with open(coords_filename) as coords_file:
        reader = csv.reader(coords_file, delimiter='\t')
        for row in reader:
            if len(row) >= 3:
                chrom, start, end,gene,a,strand = row
                if start.isdigit() and end.isdigit():
                    result.append((chrom, int(start), int(end),gene))
            else:
                logging.warn("Skipping invalid coordinate: {}".format(row))
    return result


def bam_name_legend(bam_filename):
    '''Take the first MAX_BAM_PREFIX characters from the prefix of the input BAM
    file name to be used in the legend of the graph.'''
    basename = os.path.basename(bam_filename)
    return basename[:MAX_BAM_PREFIX]

s = ['154063861','154088486']
e = ['-1','-1']


def plot_coverage(coords, bams, sample):
	'''Given the name of a DNA coordinates firl and a list of bam file names,
	plot the read aligment coverage for each bam file for each coordinate.
	One graph per coordinate will be generated. The coverage for each
 	BAM file for a given coordinate will be plotted on the same graph.
	The coordinates file should be in TSV format.'''
	coords = get_coords(coords)
	for chrom,start,end,gene in coords:
		logging.info("processing coord {} {} {}".format(chrom, start, end))
		# Start plotting the graph and generate a name for the output file
		graph_filename = start_graph(chrom, start, end,sample,gene)
		coords_range = range(start, end+1)
		#print chrom + "\t" + str(start) + "\t" + str(end) + "\t"+ str(min(coords_range)) + "\t" + str(max(coords_range))
		for bam_filename in bams:
            		# interval tree tracks the start and end mapped coordinates
            		# of each read in the bam file that lies within our region
            		# of interest.
			interval_tree = IntervalTree()
			with pysam.Samfile(bam_filename, "rb") as bam:
				logging.info("processing bam file {}".format(bam_filename))
				# Collect all the reads from the BAM file which lie in
				# the region of interest.
				# fetch uses 0-based indexing. Our input coordinates are
				# in 1-based coordinates.
				reads = bam.fetch(chrom, start-1, end-1)
				# Insert the start and end of each aligned read into the
				# interval tree.
				for read in reads:
					#print read.positions
					if len(read.positions) > 0:
						# Add 1 to convert from 0-based to 1-based coordinates
						first_pos = read.positions[0] + 1
						last_pos = read.positions[-1] + 1
						interval_tree.add(first_pos, last_pos, None)
						#For each base position in our region of interest,
						#count the number of reads which overlap this position.
						#This computes the coverage for each position in the region.
			counts = [len(interval_tree.find(pos, pos)) for pos in coords_range]
			x_min = []
			x_max = []
			y_val = []    	
			for x, slice in zip(slices(coords_range),slices(counts)):
				y = sum(slice)/float(len(slice))
				y_val.append(y)
				x_min.append(min(x))
				x_max.append(max(x))
				print('{} {} {} {}'.format(graph_filename,min(x),max(x),y))
				
			# Plot the coverage information for this bam file
			legend_text = bam_name_legend(bam_filename)
			try:
				y_val1 = [log(y,10) for y in y_val]
				plot_graph(y_val1, x_min, sample)
			except ValueError:
#				plt.plot(s, e ,'ro ', label='zeros')
				plot_graph1(y_val, x_min, sample)
			#	plt.plot(s, e ,'ro ', label='zeros')
			#Close the drawing of the graph for this set of coordinates
		end_graph(graph_filename,sample)


def start_graph(chrom, start, end, sample,gene):
    '''Begin plotting of a graph for a given coordinate.'''
   # plt.ylabel("Coverage in log scale")
   # plt.xlabel("Position")
    plt.title("Sample:{} Gene:{}\nAmplicon Position:\n{}:{}-{}".format(sample,gene,chrom,start,end))
    return "{}.{}.{}.{}.{}.png".format(sample,gene,chrom, start, end)


def plot_graph(counts, coords_range, sample):
    '''Plot the coverage data for a given bam file.'''
    plt.plot(coords_range, counts, label=sample)
    plt.ylabel("Coverage in log scale")
    plt.xlabel("Position")

def plot_graph1(counts, coords_range, sample):
    '''Plot the coverage data for a given bam file.'''
    plt.plot(coords_range, counts, label=sample)
    plt.ylabel("Coverage")
    plt.xlabel("Position")


#    print str(len(coords_range)) + "\t" + str(len(counts))
   # plt.xticks(np.arange(min(coords_range), max(coords_range)+1, 1000.0))


def end_graph(filename,sample):
    '''Begin plotting of a graph for a given coordinate.'''
    lgd = plt.legend(title="sample",loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(filename,bbox_extra_artists=(lgd,),dpi=300,bbox_inches='tight')
    plt.close()


def start_log(log):
    '''Initiate program logging. If no log file is specified then
    log output goes to DEFAULT_LOGFILE.'''
    logging.basicConfig(
        filename=log,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    # Log the command line that was used to run the program
    logging.info('command line: {0}'.format(' '.join(sys.argv)))


def main():
    '''Program entry point.'''
    args = parse_args()
    start_log(args.log)
    plot_coverage(args.coords, args.bams,args.sample)


if __name__ == '__main__':
    main()
