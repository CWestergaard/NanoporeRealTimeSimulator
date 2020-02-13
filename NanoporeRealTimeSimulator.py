"""
Program: CoverageTimesort
Description: Simulate real-time data generation by sorting input fastq.gz file 
             into smaller files. Files can be sorted based on average genome coverage or sequencing time.
             Reads are sorted based on their timestamp (Nanopore).
Version: 1.0
Author: Casper Westergaard

Example string:
python3 CoverageTimesort.py -i /srv/data/AS/CASW/data/q8/CPO20160077/CPO20160077.chop.q8.fastq.gz
 -o /srv/data/AS/CASW/data/q8/CPO20160077/realtime/ -gs 5m -cl 1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60
"""

#Import libraries
import datetime
import gzip
import sys
import os
import argparse
import shutil

###########################################################################
# GET INPUT
###########################################################################

# Input from commandline
parser = argparse.ArgumentParser(
        description='Simulate real-time data generation by sorting input fastq.gz file' 
                     'into smaller files. Files can be sorted based on average genome coverage or sequencing time.'
                     'Reads are sorted based on their timestamp (Nanopore).')
parser.add_argument('-i', type=str, dest='input_filename', 
                    help='Input fastq.gz file containing Nanopore reads', required=True)
parser.add_argument('-o', type=str, dest='output_folder', 
                    help='Path to output folder', required=True)
parser.add_argument('-gs', type=str, dest='genome_size_string', 
                    help='Size of genome', required=True)
parser.add_argument('-alb', type=str, dest='albacore', 
                    help='If albacore is used for basecalling, use -alb yes', required=False)
parser.add_argument('-cl', type=str, dest='coverage_list', 
                    help='List of coverages, a new file will be created for each element,'
                    ' coverages should be comma-delimitered and ordered from lowest to highest', required=False)
parser.add_argument('-tl', type=str, dest='time_list', 
                    help='List of times (in minutes), a new file will be created for each element,'
                    ' times should be comma-delimitered and ordered from lowest to highest', required=False)
args = parser.parse_args()


# Check that input arguments are valid
coverage_flag = 0
time_flag = 0

# Convert input genome size to int
if args.genome_size_string[-1].upper() == 'K':
    genome_size = int(float(args.genome_size_string[0:-1])) * 1000
elif args.genome_size_string[-1].upper() == 'M':
    genome_size = int(float(args.genome_size_string[0:-1])) * 1000000
elif args.genome_size_string[-1].upper() == 'G':
    genome_size = int(float(args.genome_size_string[0:-1])) * 1000000000
else:
    genome_size = int(args.genome_size_string)
    
if args.time_list:
    time_flag = 1
    # Check thal all times in list are ints and ordered
    try:
        args.time_list = [int(x) for x in args.time_list.split(',')]
    except ValueError as err:
        sys.exit('Times must be integer-values.')
    for i, cov in enumerate(args.time_list):
        if i == len(args.time_list)-1:
            break
        else:
            if cov > args.time_list[i + 1]:
                sys.exit('Times must be ordered from lowest to highest.')  
                
if args.coverage_list:  
    coverage_flag = 1      
    # Check thal all coverages in list are ints and ordered
    try:
        args.coverage_list = [int(x) for x in args.coverage_list.split(',')]
    except ValueError as err:
        sys.exit('Coverages must be integer-values.')
    for i, cov in enumerate(args.coverage_list):
        if i == len(args.coverage_list)-1:
            break
        else:
            if cov > args.coverage_list[i + 1]:
                sys.exit('Coverages must be ordered from lowest to highest.')
        
if time_flag == 0 and coverage_flag == 0:
    sys.exit('Arguments must only include either a coverage list (-cl) and genome size (-gs), or a time list (-tl).')
if time_flag == 1 and coverage_flag == 1:
    sys.exit('Arguments should only include a coverage list (-cl) or a time list (-tl), not both.')    
               
# Check that input file exists
if not os.path.exists(args.input_filename):
    sys.exit('Input file:\n{}\ndoes not exist.'.format(args.input_filename))
 
# Create output directory
args.output_folder += '/'
os.makedirs(args.output_folder, exist_ok=True)
###########################################################################
# Load all reads into memory and sort them based on their timestamp
###########################################################################  

# Initialize variables
reads = list()
sequence_flag = 0
total_bases = 0

# Difference in header between Albacore and Guppy reads. Guppy assumed as default.
if args.albacore:
    if args.albacore.lower() == 'yes':
        time_index = 4
else:
    time_index = 5
    
# Go through each read, save the timestamp, number of bases and the total read
with gzip.open(args.input_filename, 'r') as infile:           
    for line in infile:
        # Sequence - Count bases
        if sequence_flag == 1:
            bases = len(line)-1
            total_bases += bases
            reads[-1][1] = bases
            sequence_flag = 0
        # Header - Save timestamp
        if line.startswith(b'@') and b'start_time' in line:      
            time_string = line.split()[time_index].split(b'=')[1].decode('ascii')   # Get time from timestamp in header
            time = datetime.datetime.strptime(time_string, '%Y-%m-%dT%H:%M:%SZ')
            reads.append([time,None,b'']) # Time, number of bases in read, total read
            sequence_flag = 1
            line_count = 0
        # Save reads for output later
        reads[-1][2] += line

  
# Sort reads based on timestamp        
reads.sort()

###########################################################################
# Figure out the subsets of reads to include in each output file
###########################################################################  

# Coverage
if coverage_flag == 1:
    total_bp = 0
    coverage_reads = list()
    coverage_counter = 0  
    # Go through each read     
    for i in range(len(reads)): 
        read_bp = reads[i][1]
        total_bp += read_bp
    
        # When reads containing enough basepairs to reach a given coverage is reached, continue to next element in the list
        if total_bp > args.coverage_list[coverage_counter]*genome_size:
            coverage_reads.append([args.coverage_list[coverage_counter],i])
            coverage_counter += 1
        
        # Reads for the given coverage intervals are found, no need to continue
        if total_bp > args.coverage_list[-1]*genome_size:
            break 
   
# Time
if time_flag == 1:
    start_time = reads[0][0]
    time_reads = list()
    time_counter = 0  
    # Go through each read     
    for i in range(len(reads)): 
        read_time = reads[i][0]
        total_time = read_time-start_time
        if total_time.days == 0:
            total_time_minutes = ((read_time-start_time).seconds)/60
        else:
            total_time_minutes = (total_time.days*24*60)+((read_time-start_time).seconds)/60
    
        # When time since start reaches a time in the time interval list, continue to next element in the list
        if total_time_minutes > args.time_list[time_counter]:
            time_reads.append([args.time_list[time_counter],i])
            time_counter += 1
        
        # Reads for the given time intervals are found, no need to continue
        if total_time_minutes > args.time_list[-1]:
            break 
      
###########################################################################
# Create output files
###########################################################################

# Create filenames
isolate_filename = os.path.basename(args.input_filename)
isolate = '.'.join(isolate_filename.split('.')[0:-2]) # To not include .fastq.gz
output_filenames = list()

# Coverage
if coverage_flag == 1:
    for i in range(len(coverage_reads)):
        if coverage_reads[i][0] < 10:
            output_filenames.append('{}{}.cov_00{}.fastq.gz'.format(args.output_folder,isolate,coverage_reads[i][0]))        
        elif coverage_reads[i][0] < 100:
            output_filenames.append('{}{}.cov_0{}.fastq.gz'.format(args.output_folder,isolate,coverage_reads[i][0]))  
        else:
            output_filenames.append('{}{}.cov_{}.fastq.gz'.format(args.output_folder,isolate,coverage_reads[i][0]))
            
    # Write to first file
    with gzip.open(output_filenames[0], 'w') as outfile:
        first_read_index = 0
        last_read_index = coverage_reads[0][1]+1
        for i in range(first_read_index,last_read_index):
            outfile.write(reads[i][2])
            
    # Copy last file and append to it            
    for i in range(1,len(coverage_reads)):
        shutil.copyfile(output_filenames[i-1], output_filenames[i])
        first_read_index = last_read_index
        last_read_index = coverage_reads[i][1]+1
        with gzip.open(output_filenames[i], 'a') as outfile:
            for j in range(first_read_index,last_read_index):
                outfile.write(reads[j][2])
                
    # Write coverage information file
    outfilename = ('{}{}.cov_stats.txt'.format(args.output_folder, isolate))
    with open(outfilename, 'w') as outfile:
        
        # Output general information
        outfile.write('Input file: ' + args.input_filename)
        first_read_time = reads[0][0]
        last_read_time = reads[-1][0]
        total_time = (datetime.datetime.min + (last_read_time-first_read_time)).time()
        outfile.write('\nTotal sequencing time: {}'.format(total_time))
        outfile.write('\nGenome size: {}'.format(args.genome_size_string))
        outfile.write('\nTotal bases in file: {}'.format(total_bases))
        outfile.write('\nTotal coverage: {}'.format(total_bases/genome_size))
        outfile.write('\nTime to reach coverage (Coverage: Hours:Minutes:Seconds):\n')
        
        # Output time for to reach coverages
        for i in range(len(coverage_reads)):
            if args.coverage_list[i] < 100:
                outfile.write('0')
            if args.coverage_list[i] < 10:
                outfile.write('0')
            coverage_last_read_time = reads[coverage_reads[i][1]][0]     
            time_since_start = (datetime.datetime.min + (coverage_last_read_time-first_read_time)).time()
            outfile.write('{}: {} \n'.format(args.coverage_list[i],time_since_start))
    
# Time
if time_flag == 1:
    for i in range(len(time_reads)):
        if time_reads[i][0] < 10:
            output_filenames.append('{}{}.time_000{}.fastq.gz'.format(args.output_folder,isolate,time_reads[i][0]))        
        elif time_reads[i][0] < 100:
            output_filenames.append('{}{}.time_00{}.fastq.gz'.format(args.output_folder,isolate,time_reads[i][0]))  
        elif time_reads[i][0] < 1000:
            output_filenames.append('{}{}.time_0{}.fastq.gz'.format(args.output_folder,isolate,time_reads[i][0]))
        else:
            output_filenames.append('{}{}.time_{}.fastq.gz'.format(args.output_folder,isolate,time_reads[i][0]))
            
    # Write to first file
    with gzip.open(output_filenames[0], 'w') as outfile:
        first_read_index = 0
        last_read_index = time_reads[0][1]
        for i in range(first_read_index,last_read_index):
            outfile.write(reads[i][2])
            
    # Copy last file and append to it            
    for i in range(1,len(time_reads)):
        shutil.copyfile(output_filenames[i-1], output_filenames[i])
        first_read_index = last_read_index
        last_read_index = time_reads[i][1]
        with gzip.open(output_filenames[i], 'a') as outfile:
            for j in range(first_read_index,last_read_index):
                outfile.write(reads[j][2])

    # Write time information file
    outfilename = ('{}{}.time_stats.txt'.format(args.output_folder, isolate))
    with open(outfilename, 'w') as outfile:
        
        # Output general information
        outfile.write('Input file: ' + args.input_filename)
        first_read_time = reads[0][0]
        last_read_time = reads[-1][0]
        total_time = (datetime.datetime.min + (last_read_time-first_read_time)).time()
        outfile.write('\nTotal sequencing time: {}'.format(total_time))
        outfile.write('\nGenome size: {}'.format(args.genome_size_string))
        outfile.write('\nTotal bases in file: {}'.format(total_bases))
        outfile.write('\nTotal coverage: {}'.format(total_bases/genome_size))
        outfile.write('\nCoverage after sequencing time (Minutes: Coverage):\n')
        
        # Output time for to reach coverages
        time_basecount = 0    
        for i in range(len(time_reads)):
            if args.time_list[i] < 1000:
                outfile.write('0')
            if args.time_list[i] < 100:
                outfile.write('0')
            if args.time_list[i] < 10:
                outfile.write('0')
               
            # Count and output bases for each time in the time list
            if i > 0:
                time_first_read_index = (time_reads[i-1][1])-1
            else:
                time_first_read_index = 0
            time_last_read_index = time_reads[i][1]
            bases = sum([reads[x][1] for x in range(time_first_read_index,time_last_read_index)])   # Bases in reads added between previous time and now
            time_basecount += bases
            time_coverage = time_basecount/genome_size
            outfile.write('{}: {} \n'.format(args.time_list[i],time_coverage))
            