# NanoporeRealTimeSimulator
Python script for simulating real-time data acquisition of Nanopore sequencing reads. The script sorts the reads in the input Nanopore fastq.gz file based on their timestamp, and creates new files based on data amount or sequencing time. These new files can then be used to evaluate how the amount of sequencing data influences the results of various analyses, similar to how a real-time analysis would change as more data is generated.
