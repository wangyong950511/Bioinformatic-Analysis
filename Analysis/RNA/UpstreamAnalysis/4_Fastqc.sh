#!/bin/bash

# FastQC
mkdir fastqcoutcom
ls *.fq.gz | parallel fastqc -o ./fastqcoutcom {}
