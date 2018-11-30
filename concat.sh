#!/bin/bash

for i in {118282..118289}; do cat `find /home/jingquma/Run_2509/ -name "$i*\_R1_001.fastq.gz"|sort` > /home/jingquma/Run_2509/JM/$i\_R1.fastq.gz; done

for i in {118282..118289}; do cat `find /home/jingquma/Run_2509/ -name "$i*\_R2_001.fastq.gz"|sort` > /home/jingquma/Run_2509/JM/$i\_R2.fastq.gz; done
