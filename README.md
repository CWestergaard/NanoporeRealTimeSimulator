# NanoporeRealTimeSimulator
Python script for simulating real-time data acquisition of Nanopore sequencing reads. The script sorts the reads in the input Nanopore fastq.gz file based on their timestamp, and creates new files based on data amount or sequencing time. These new files can then be used to evaluate how the amount of sequencing data influences the results of various analyses, similar to how the results of a real-time analysis would change as more data is generated.

## Requirements

- Input file: Nanopore reads in fastq.gz format

## Usage
Data amount (Average genome coverage, based on genome size)
```
python3 NanoporeRealTimeSimulator.py [-i <Input file>] [-o <Output folder>] [-gs <Genome size>] [-cl <Coverage list>]
```
Sequencing time (Minutes)
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

## Output

Data amount (Average genome coverage) - Example: -gs 5M -cl 1,2,3,4,5
- Coverage1.fastq.gz file containing the first 5 Mbp of sequencing reads.
- Coverage2.fastq.gz file containing the first 10 Mbp of sequencing reads.
- Coverage3.fastq.gz file containing the first 15 Mbp of sequencing reads.
- Coverage4.fastq.gz file containing the first 20 Mbp of sequencing reads.
- Coverage5.fastq.gz file containing the first 25 Mbp of sequencing reads.
- Cov_stats.txt file containing general information as well as the time required to reach the given average genome coverages.

Sequencing time (Minutes) - Example: -gs 5M -tl 10,20,30,40,50
- Time10.fastq.gz file containing the reads sequenced in the first 10 minutes after the first read is sequenced.
- Time20.fastq.gz file containing the reads sequenced in the first 20 minutes after the first read is sequenced.
- Time30.fastq.gz file containing the reads sequenced in the first 30 minutes after the first read is sequenced.
- Time40.fastq.gz file containing the reads sequenced in the first 40 minutes after the first read is sequenced.
- Time50.fastq.gz file containing the reads sequenced in the first 50 minutes after the first read is sequenced.
- Time_stats.txt file containing general information as well as the average genome coverage reached at the given sequencing times.
