# NanoporeRealTimeSimulator
Python script for simulating real-time data acquisition of Nanopore sequencing reads. The script sorts the reads in the input Nanopore fastq.gz file based on their timestamp, and creates new files based on data amount or sequencing time. These new files can then be used to evaluate how the amount of sequencing data influences the results of various analyses, similar to how the results of a real-time analysis would change as more data is generated.

## Requirements

- Input file: Nanopore reads in fastq.gz format

## Usage
Data amount (Number of bases)
```
python3 NanoporeRealTimeSimulator.py [-i <Input file>] [-o <Output folder>]  [-sl <Size list>]
```
Data amount (Average genome coverage, based on genome size)
```
python3 NanoporeRealTimeSimulator.py [-i <Input file>] [-o <Output folder>] [-gs <Genome size>] [-cl <Coverage list>]
```
Sequencing time (Minutes)
```
python3 NanoporeRealTimeSimulator.py [-i <Input file>] [-o <Output folder>] [-tl <Time list>]
```

## Input arguments

Required arguments
```
-i <Path to input file. Nanopore reads in .fastq.gz format>
-o <Path to output folder>
```
Mutually exclusive arguments (Either -sl, -cl or -tl must be included).
```
-sl <Size list. Used to decide which reads to include in the smaller Nanopore fastq files, based on the data amount. Data amount is given as the amount of bases. Values can be given as the full number (5000000), or as Kilo- (5000K), Mega- (5M) or Gigabases (0.005G). Values need to be comma-delimetered, eg. 1M,2M,3M,4M,5M>
-cl <Coverage list. Used to decide which reads to include in the smaller Nanopore fastq files, based on the data amount. Data amount is given as the average genome coverage (-gs), so a coverage of one corresponds to the amount of basepairs equal to the genome size given. Values need to be comma-delimetered, eg. 1,2,3,4,5>
-tl <Time list. Used to decide which reads to include in the smaller Nanopore fastq files, based on the sequencing time. Values need to be comma-delimetered, eg. 10,20,30,40,50,60. These values correspond to the amount of Minutes since the first read is sequenced.>
```
Required arguments with -cl
```
-gs <Genome size. As an example, an approximate Genome size for Gram-negative bacteria is 5.000.000 bp. Genome size can be given as the full number (5000000), or as Kilo- (5000K), Mega- (5M) or Gigabases (0.005G).>
```
## Output
Data amount (Number of bases) - Example: -sl 1M,2M,3M,4M,5M
- size_1M.fastq.gz file containing the first 1 Mbp of sequencing reads.
- size_2M.fastq.gz file containing the first 2 Mbp of sequencing reads.
- size_3M.fastq.gz file containing the first 3 Mbp of sequencing reads.
- size_4M.fastq.gz file containing the first 4 Mbp of sequencing reads.
- size_5M.fastq.gz file containing the first 5 Mbp of sequencing reads.
- size_stats.txt file containing general information as well as the time required to reach the given number of bases.
- 
Data amount (Average genome coverage) - Example: -gs 5M -cl 1,2,3,4,5
- coverage1.fastq.gz file containing the first 5 Mbp of sequencing reads.
- coverage2.fastq.gz file containing the first 10 Mbp of sequencing reads.
- coverage3.fastq.gz file containing the first 15 Mbp of sequencing reads.
- coverage4.fastq.gz file containing the first 20 Mbp of sequencing reads.
- coverage5.fastq.gz file containing the first 25 Mbp of sequencing reads.
- coverage_stats.txt file containing general information as well as the time required to reach the given average genome coverages.

Sequencing time (Minutes) - Example: -tl 10,20,30,40,50
- time10.fastq.gz file containing the reads sequenced in the first 10 minutes after the first read is sequenced.
- time20.fastq.gz file containing the reads sequenced in the first 20 minutes after the first read is sequenced.
- time30.fastq.gz file containing the reads sequenced in the first 30 minutes after the first read is sequenced.
- time40.fastq.gz file containing the reads sequenced in the first 40 minutes after the first read is sequenced.
- time50.fastq.gz file containing the reads sequenced in the first 50 minutes after the first read is sequenced.
- time_stats.txt file containing general information as well as the number of bases sequenced at the given sequencing times.
