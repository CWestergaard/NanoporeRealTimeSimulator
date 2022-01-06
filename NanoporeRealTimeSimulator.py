"""
Program: NanoporeRealTimeSimulator
Description: 
Simulate the ability of Nanopore to generate sequencing reads for analysis in 
real-time, by sorting reads in the input fastq.gz file into smaller files, 
based on their timestamp.
Files can be sorted based on:
    - Coverage (Requires genome size)
    - Sequencing time (Minutes)
    - Size (bp, Kbp, Mbs or Gbp)

Version: 2.0
Author: Casper Westergaard

Example of arguments for the 3 modes:
Coverage:
    - Requires genome size (gs) argument in bp, Kbp (K), Mbp (M) or Gbp (G).
    - Requires list of coverages for the output files (cl).
    
python3 NanoporeRealTimeSimulator.py -i SampleReads.fastq.gz
 -o SampleReads/realtime/ -gs 5M -cl 1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60

Time:
    - Requires list of sequencing times (tl) for the output files, in minutes.
    
python3 NanoporeRealTimeSimulator.py -i SampleReads.fastq.gz
 -o SampleReads/realtime/ -tl 10,20,30,40,50,60

Size:
    - Requires list of sizes (sl) for the output files, in bp, Kbp (K), Mbp (M) or Gbp (G).
    
python3 NanoporeRealTimeSimulator.py -i SampleReads.fastq.gz
 -o SampleReads/realtime/ -sl 10M,20M,30M,40M,50M,100M,200M
"""
###########################################################################
# IMPORT LIBRARIES
###########################################################################

#Import libraries
import datetime
import gzip
import sys
import os
import argparse
import shutil

###########################################################################
# FUNCTIONS
###########################################################################

def baseSize_convert(baseSize_string):
    """
        Takes a basesize as string, with potential suffixes, and converts it to an int.
        
        parameters:
            baseSize_string = Size of genome as bp, Kbp, Mbp or Gbp.
        returns:
            baseSize = Size of base as int.
    """    
    # Convert input genome size to int
    if baseSize_string[-1].upper() == 'K':
        baseSize = float(baseSize_string[0:-1]) * 1000
    elif baseSize_string[-1].upper() == 'M':
        baseSize = float(baseSize_string[0:-1]) * 1000000
    elif baseSize_string[-1].upper() == 'G':
        baseSize = float(baseSize_string[0:-1]) * 1000000000
    else:
        baseSize = float(baseSize)
    
    return int(baseSize)

def stringConvert_intList(inputString):
    """
        Takes a string of comma-delimitered inputs and converts them to a
        list of integers. The elements are then sorted based on size.
        
        parameters:
            inputString = String of comma-delimitered elements to convert
        returns:
            intList = List of elements sorted by size
    """    
    
    intList = [int(x) for x in inputString.split(',')]
    intList.sort()

    return intList

###########################################################################
# GET INPUT
###########################################################################

# Input from commandline
parser = argparse.ArgumentParser(
        description='Simulate real-time data generation by sorting input fastq.gz file' 
                     'into smaller files. Files can be sorted based on average genome coverage or sequencing time.'
                     'Reads are sorted based on their timestamp (Nanopore).')

# Required arguments
parser.add_argument('-i', type=str, dest='input_filename', 
                    help='Input fastq.gz file containing Nanopore reads', required=True)
parser.add_argument('-o', type=str, dest='output_folder', 
                    help='Path to output folder', required=True)

# Only required with -cl
parser.add_argument('-gs', type=str, dest='genome_size_string', 
                    help='Basesize of genome, basesize can be listed as bp, Kbp (K), Mbp (M) or Gbp (G). '
                    '-gs is required with -cl and -tl, but not with -sizes.',
                    required='-cl' in sys.argv)

# Mutually exclusive arguments, one of them is required
group = parser.add_mutually_exclusive_group(required=True)

group.add_argument('-cl', type=str, dest='coverage_list', 
                    help='List of coverages, a new file will be created for each element. '
                    'Coverages should be comma-delimitered integers.')

group.add_argument('-tl', type=str, dest='time_list', 
                    help='List of times (in minutes), a new file will be created for each element.'
                    'Times should be comma-delimitered integers.')

group.add_argument('-sl', type=str, dest='size_list', 
                    help='List of basesizes, a new file will be created for each element. '
                    'Basesizes should be comma-delimitered numbers. '
                    'Basesizes can be listed as bp, Kbp (K), Mbp (M) or Gbp (G).')

args = parser.parse_args()

###########################################################################
# Check input arguments
###########################################################################

# Check if user chose to base the analysis on time, coverage or size
# and check for valid input

# Time
if args.time_list:
    mode = "time"
    # Convert input argument string to sorted list of integers
    try:
        args.time_list = stringConvert_intList(args.time_list)
        input_list = args.time_list
    except ValueError:
        sys.exit('Times in -tl must be integer-values.')  
        
        
# Coverage
if args.coverage_list:
    mode = "coverage"
    # Convert input argument string to sorted list of integers
    try:
        args.coverage_list = stringConvert_intList(args.coverage_list)
        input_list = args.coverage_list
    except ValueError:
        sys.exit('Coverages in -cl must be integer-values.')
  
        
# Basesize        
if args.size_list:
    mode = "size"
    # Convert input argument string to sorted list of integers
    size_list_converted = []
    for size in args.size_list.split(','):      
        try:
            size = baseSize_convert(size)
        except:
            sys.exit('One or more basesizes in -sl is not given as a number '
                     '(K, M and G suffixes excluded)')
        size_list_converted.append(size)        
    size_list_converted.sort()
    input_list = args.size_list.split(",")


# Convert input genome size to int
if args.genome_size_string:
    try:
        genome_size = baseSize_convert(args.genome_size_string)
    except:
        sys.exit("Genomesize is not given as a number (K, M and G suffixes excluded)")

    
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

# Check for placement of timestamp in read headers, due to potential difference in header versions.
# Albacore and new Guppy versions has the time at index 4 (0-based indexing)
# Older Guppy versions had the time at index 5.
with gzip.open(args.input_filename, 'r') as infile:           
    for line in infile:
        if line.startswith(b'@') and b'start_time' in line: 
            headerSplit = line.split(b" ")
            break
    for i, element in enumerate(headerSplit):
        if b'start_time' in element:
            time_index = i
        
            
# Go through each read, save the timestamp, number of bases and the total read including header
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
        # Save line for output later
        reads[-1][2] += line
 
# Sort reads based on timestamp        
reads.sort()

###########################################################################
# Figure out the subsets of reads to include in each output file
###########################################################################  
reads_output = list()   # Fill with the number of reads required to reach each value in the input list
counter = 0 

# Coverage or Size
if mode == "coverage" or mode == "size":
    # Create list of required number of bases for each element in input list
    if mode == "coverage":
        required_bp = [i * genome_size for i in args.coverage_list]
    elif mode == "size":
        required_bp = size_list_converted
        
    total_bp = 0
    # Go through each read     
    for i in range(len(reads)): 
        read_bp = reads[i][1]
        total_bp += read_bp
    
        # When reads containing enough bases to reach the requirement is found,
        # then continue to next element in the list.
        if total_bp > required_bp[counter]:
            reads_output.append([input_list[counter],i])
            counter += 1
        
            # Reads to complete all elements in the input list are found, no need to continue
            if total_bp > required_bp[-1]:
                break 

# Time
if mode == "time":
    start_time = reads[0][0]
    # Go through each read     
    for i in range(len(reads)): 
        read_time = reads[i][0]
        total_time = read_time-start_time   # Time since start
        # Convert to minutes
        if total_time.days == 0:
            total_time_minutes = total_time.seconds/60
        else:
            total_time_minutes = (total_time.days*24*60)+(total_time.seconds)/60
    
        # When time since start reaches a time in the time interval list, continue to next element in the list
        if total_time_minutes > input_list[counter]:
            reads_output.append([input_list[counter],i-1])
            counter += 1
        
            # Reads for the given time intervals are found, no need to continue
            if total_time_minutes > input_list[-1]:
                break 
  
###########################################################################
# Create output files
###########################################################################

# Create list of output filenames
isolate_filename = os.path.basename(args.input_filename)
isolate = '.'.join(isolate_filename.split('.')[0:-2]) # To not include .fastq.gz
output_filenames = list()
affix_list = list() # Append 0's in front of values to make all values equal in length
affix_len = len(str(reads_output[-1][0]))
for i in range(len(reads_output)):
    value = reads_output[i][0]
    value_len = len(str(value))
    affix = (affix_len-value_len)*"0"
    output_filenames.append('{}{}.{}_{}{}.fastq.gz'.
                                    format(args.output_folder,isolate,mode,affix,value))
    affix_list.append(affix)
          
# Write to first file
with gzip.open(output_filenames[0], 'w') as outfile:
    first_read_index = 0
    last_read_index = reads_output[0][1]
    for i in range(first_read_index,last_read_index+1):
        outfile.write(reads[i][2])
        
# Copy last file and append to it            
for i in range(1,len(reads_output)):
    shutil.copyfile(output_filenames[i-1], output_filenames[i])
    first_read_index = last_read_index+1
    last_read_index = reads_output[i][1]
    with gzip.open(output_filenames[i], 'a') as outfile:
        for j in range(first_read_index,last_read_index+1):
            outfile.write(reads[j][2])
          
# Write summary information file
outfilename = ('{}{}.{}_stats.txt'.format(args.output_folder, isolate, mode))
with open(outfilename, 'w') as outfile:
    
    # Output general information
    outfile.write('Input file: ' + args.input_filename)
    first_read_time = reads[0][0]
    last_read_time = reads[-1][0]
    total_time = (datetime.datetime.min + (last_read_time-first_read_time)).time()
    outfile.write('\nTotal sequencing time: {}'.format(total_time))
    outfile.write('\nTotal bases in file: {}'.format(total_bases))
    
    if mode == "coverage":
        outfile.write('\nGenome size: {}'.format(args.genome_size_string))
        outfile.write('\nTotal coverage: {}'.format(total_bases/genome_size))
        outfile.write('\nTime to reach coverage (Coverage: Hours:Minutes:Seconds):\n')
        
    elif mode == "time":
        outfile.write('\nCoverage after sequencing time (Minutes: Coverage):\n')
        
    elif mode == "size":
        outfile.write('\nTime to reach bases (Bases: Hours:Minutes:Seconds):\n')
    
    # Output time for to reach coverages or bases, depending on mode
    if mode == "coverage" or mode == "size":
        for i in range(len(reads_output)):
            read_index = reads_output[i][1] # Index of last included read
            last_read_time = reads[read_index][0]   # Timestamp for last included read     
            time_since_start = (datetime.datetime.min + (last_read_time-first_read_time)).time()
            outfile.write('{}{}: {} \n'.format(affix_list[i],input_list[i],time_since_start))
    
    # Output number of bases for each given time
    if mode == "time":
        basecount = 0    
        for i in range(len(reads_output)): 
            # Count and output bases for each time in the time list
            if i > 0:
                first_read_index = reads_output[i-1][1]+1
            else:
                first_read_index = 0
            last_read_index = reads_output[i][1]
            # Calculate bases in reads added between previous time and now
            bases = sum([reads[x][1] for x in range(first_read_index,last_read_index+1)])
            basecount += bases
            basecount_string = str(basecount/1000000)+"Mbp"
            outfile.write('{}{}: {} \n'.format(affix_list[i],input_list[i],basecount_string))
