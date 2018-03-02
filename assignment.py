#! /groups/quantitativegenomics/home/qguser/bin/anaconda3/bin/python

import sys,gzip,os
from collections import defaultdict

fastq_dir = sys.argv[1]
ssPath = sys.argv[2]
platePath = sys.argv[3]
org = sys.argv[4]
output_dir = sys.argv[5]

scriptDir = os.path.abspath('/'.join(sys.argv[0].split('/')[:-1]))
star = scriptDir + '/dep/STAR/source/STAR'
samtools = scriptDir + '/dep/samtools-1.5/samtools'
trimmomatic = scriptDir + '/dep/Trimmomatic-0.36/trimmomatic-0.36.jar'
starcode = scriptDir + '/dep/starcode/starcode'
pigz = scriptDir + '/dep/pigz/pigz'
output_dir = os.path.abspath(output_dir)

stardb = scriptDir + '/db/' + org
gtfPath = scriptDir + '/db/' + org + '.gtf'
qPath = scriptDir + '/Q.index'

qFile = open(qPath)
ssData = open(ssPath).read().strip()
plateFile = open(platePath)
qIndex = {}
for line in qFile:
	sampBC, qid, wid = line.strip().split('\t')
	qIndex[(qid,wid)] = sampBC

libBC = {}

libs = [x.split(',') for x in ssData.split('[Data]')[-1].strip().split('\n')\
if x.split(',')[0].strip() != ''][1:]

pidIndex = {}
for lib in libs:
	libName,_,_,_,_,bc,_,_ = lib
	pid = libName.split('_')[1]
	pidIndex[pid] = libName
	libBC[libName] = bc

bcIndex = {}
wellIndex = defaultdict(set)
sampIndex = {}
for line in plateFile:
	wid, pid, qid, name = line.strip().split('\t')
	name = '_'.join(name.split())

	if pid in pidIndex:
		lib = pidIndex[pid]
		lib_bc = libBC[lib]
		well_bc = qIndex[(qid,wid)]

		bcIndex[lib_bc] = (pid,qid)
		wellIndex[well_bc].add(lib_bc)
		sampIndex[(lib_bc,well_bc)] = (name,wid)

qFile.close()
plateFile.close()

libs = set(['_'.join(x.split('_')[:-1]) for x in os.listdir(fastq_dir) if x.find('fastq.gz') != -1])
libs = list(libs)

indexFile = open(output_dir + '/index.cluster')
indexCluster = {}
for line in indexFile:
	rep, count, members = line.strip().split()
	members = members.split(',')
	for member in members:
		indexCluster[member] = rep
indexFile.close()

bcFile = open(output_dir + '/bc.cluster')
bcCluster = {}
for line in bcFile:
	rep, count, members = line.strip().split()
	members = members.split(',')
	for member in members:
		bcCluster[member] = rep
bcFile.close()

logFile = open(output_dir + '/reads.unassigned.log','w')
statFile = open(output_dir + '/reads.assignment.stat','w')

stats = defaultdict(lambda : defaultdict(int))
for lib in libs:
	bcPath = os.path.abspath(sys.argv[1] + '/' + lib + '_1.fastq.gz')
	seqPath = os.path.abspath(sys.argv[1] + '/' + lib + '_2.fastq.gz')

	bcFile = gzip.open(bcPath,'r')
	seqFile = gzip.open(seqPath,'r')

	for header in bcFile:
		header = header.decode('utf8').strip().split()
		seq = next(bcFile).decode('utf8').strip()
		next(bcFile)
		next(bcFile)

		next(seqFile)
		read_seq = next(seqFile).decode('utf8').strip()
		next(seqFile)
		read_qual = next(seqFile).decode('utf8').strip()

		rid = header[0]
		iindex = header[1].split(':')[3]
		filtered = header[1].split(':')[1]
		bc = seq[1:9]
		umi = seq[9:19]
		polyt = seq[19:]

		rep_iindex = 'null'
		if iindex in indexCluster:
			rep_iindex = indexCluster[iindex]

		rep_lib = 'null'
		rep_q = 'null'
		if rep_iindex in bcIndex:
			rep_lib,re_q = bcIndex[rep_iindex]

		rep_bc = 'null'
		if bc in bcCluster:
			rep_bc = bcCluster[bc]

		contamination = 0
		lib_exist = 0
		samp_exist = 0

		if rep_lib != 'null':
			lib_exist =  1

		rep_samp = 'null'
		rep_wid = 'null'
		if rep_bc in wellIndex:
			samp_exist = 1
			if not rep_iindex in wellIndex[rep_bc]:
				contamination = 1
			else:
				rep_samp, rep_wid = sampIndex[(rep_iindex,rep_bc)]

		header = '\t'.join([rid, filtered, \
		str(lib_exist), iindex, rep_iindex, rep_lib, \
		str(samp_exist), bc, rep_bc, rep_samp, \
		str(contamination), umi, polyt])

		if lib_exist == 1 and samp_exist == 1 and contamination == 0:
			print(':'.join([rid, rep_lib, rep_wid, rep_samp, umi, polyt])\
			 + '\n' + read_seq + '\n+\n' + read_qual)

			stats[rep_lib][rep_wid + ':' + rep_samp] += 1
		else:
			logFile.write('\t'.join([rid, filtered, \
		str(lib_exist), iindex, rep_iindex, rep_lib, \
		str(samp_exist), bc, rep_bc, rep_wid, rep_samp, \
		str(contamination), umi, polyt]) + '\n')

			if lib_exist == 1:
				if samp_exist == 1:
					stats[rep_lib]['contamination'] += 1
				else:
					stats[rep_lib]['unknown'] += 1
			else:
				if samp_exist == 1:
					stats['unknown']['contamination'] += 1
				else:
					stats['unknown']['unknown'] += 1

logFile.close()
for lib_stat, lib_data in stats.items():
	for lib_bc, count in lib_data.items():
		statFile.write('\t'.join([lib_stat,lib_bc,str(count)]) + '\n')
statFile.close()


