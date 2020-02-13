# NanoporeRealTimeSimulator
Python script for simulating real-time data acquisition of Nanopore sequencing reads. The script sorts the reads in the input Nanopore fastq.gz file based on their timestamp, and creates new files based on data amount or sequencing time. These new files can then be used to evaluate how the amount of sequencing data influences the results of various analyses, similar to how the results of a real-time analysis would change as more data is generated.

## Requirements

- Input file: Nanopore reads in fastq.gz format

## Usage
Data amount (Average genome coverage, based on genome size)
```
python3 NanoporeRealTimeSimulator.py [-i <Input file>] [-o <Output folder>] [-gs <Genome size>] [-cl <Coverage list>]
```
Sequencing time
```
python3 NanoporeRealTimeSimulator.py [-i <Input file>] [-o <Output folder>] [-gs <Genome size>] [-tl <Time list>]
```

## Input arguments
Required arguments
```
-i <Path to input file. Nanopore reads in .fastq.gz format>
-o <Path to output folder>
-gs <Genome size. As an example, an approximate Genome size for Gram-negative bacteria is 5.000.000 bp. Genome size can be given as the full number (5000000), or as Giga- (0.005G), Mega- (5M) or Kilobases (5000K).>
```
Optional arguments (Either -cl or -tl must be included).
```
-cl <Coverage list. Used to decide which reads to include in the smaller Nanopore fastq files, based on the data amount. Data amount is given as the average genome coverage, so a coverage of one corresponds to the amount of basepairs equal to the genome size given. Values should be integers ordered from lowest to highest and comma-delimetered, eg. 1,2,3,4,5>
-tl <Time list. Used to decide which reads to include in the smaller Nanopore fastq files, based on the sequencing time. Values should be integers ordered from lowest to highest and comma-delimetered, eg. 10,20,30,40,50,60. These values correspond to the amount of Minutes since the first read is sequenced.>
-alb <Script is written with Guppy basecalling in mind. If the input reads are basecalled using Albacore, use '-alb yes'.>
```
