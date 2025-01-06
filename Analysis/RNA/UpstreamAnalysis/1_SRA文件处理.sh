#!/bin/bash

parallel-fastq-dump --sra-id SRR26521328.lite.1 --split-files --gzip --threads 20
parallel-fastq-dump --sra-id SRR26521329.lite.1 --split-files --gzip --threads 20
parallel-fastq-dump --sra-id SRR26521330.lite.1 --split-files --gzip --threads 20
parallel-fastq-dump --sra-id SRR26521331.lite.1 --split-files --gzip --threads 20
parallel-fastq-dump --sra-id SRR26521332.lite.1 --split-files --gzip --threads 20

