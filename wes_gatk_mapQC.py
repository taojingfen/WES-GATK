#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   wes_gatk_mapQC.py
@Time    :   2020/10/27 10:27:52
@Author  :   tao_jingfen
@Contact :   taojingfen@foxmail.com
@Desc    :   Do the QC of sequencing data
'''

import os
import sys
import argparse
import pysam
import time
import load_config
import cmd_logger
import multiprocessing as mlt

def run_mismatchrate(finalbam):
	read_NOmismatch = 0
	f = pysam.AlignmentFile(finalbam, 'rb')
	for line in f:
		if not line.is_unmapped:
			if line.get_tag('NM') > 0:
				pass
			else:
				read_NOmismatch += 1
	return read_NOmismatch

def run_samstats(finalbam, outfolder, sample, logger):
	cmd = ("samtools stats {finalbam} > {outfolder}/{sample}_samtools_stats"
	      ).format(finalbam=finalbam, outfolder=outfolder, sample=sample)
	tag = cmd_logger.run_cmd(cmd, logger)
	if not tag:
		logger.error("Failed samtools stats")
		sys.exit(1)
	else:
		logger.info("Finished samtools stats")

def run_basedepth(finalbam, region, genomefile, outfolder, sample, logger):
	cmd = ("bedtools coverage -a {region} -b {finalbam} -g {genomefile} -sorted -hist > "
	      "{outfolder}/{sample}_base_depth.txt").format(finalbam=finalbam, 
		  region=region, genomefile=genomefile,outfolder=outfolder, sample=sample)
	tag = cmd_logger.run_cmd(cmd, logger)
	if not tag:
		logger.error("Failed bedtools coverage")
		sys.exit(1)
	else:
		logger.info("Finished bedtools coverage")

def run_regionreads(finalbam, region, outfile, logger):
	cmd = ("samtools view --threads 5 -F 260 -L {region} -c {finalbam} > {outfile}"
	      ).format(finalbam=finalbam, region=region, outfile=outfile)
	tag = cmd_logger.run_cmd(cmd, logger)
	if not tag:
		logger.error("Failed samtools view")
		sys.exit(1)
	else:
		logger.info("Finished samtools view")

def info_samstats(infile):
	fin = open(infile)
	reads=0
	mapped=0
	proper_mapped=0
	error_rate=0
	reads_len=0
	insert_size=0
	for i in fin:
		if i.startswith('SN'):
			i = i.strip().split('\t')
			if i[1] == 'raw total sequences:':
				reads = int(i[2])
			elif i[1] == 'reads mapped:':
				mapped = int(i[2])
			elif i[1] == 'reads properly paired:':
				proper_mapped = int(i[2])
			elif i[1] == 'error rate:':
				error_rate = float(i[2])
			elif i[1] == 'average length:':
				reads_len = int(i[2])
			elif i[1] == 'insert size average:':
				insert_size = float(i[2])
	return reads,mapped,proper_mapped,error_rate,reads_len,insert_size

def summary_depth(infile):
	fin = open(infile)
	cover_list = []
	for i in range(21):
		cover_list.append(0)
	total_Base = 0
	map_Base = 0
	for line in fin:
		if line.startswith("all"):
			cols = line.split()
			map_Base += int(cols[1])*int(cols[2])
			if int(cols[1])<=20:
				cover_list[int(cols[1])]=float(cols[-1])
			else:
				cover_list.append(float(cols[-1]))
			total_Base = int(cols[-2])
	cover1 = sum(cover_list[1:])
	cover4 = sum(cover_list[4:])
	cover10 = sum(cover_list[10:])
	cover20 = sum(cover_list[20:])
	return [cover1,cover4,cover10,cover20],map_Base,total_Base

def run_qc(finalbam: str, markdupmetric: str, sample: str,
           outfolder: str, config: str, logger: str):
	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
	para_database = load_config.get_config(config, 'database')
	targetbed = para_database['targetbed']
	genomefile = para_database['genomefile']

	try:
		pool = mlt.Pool(6)
		x = pool.apply_async(run_mismatchrate, args=(finalbam,))
		read_NOmismatch = x.get()
		pool.apply_async(run_samstats, args=(finalbam, outfolder, sample, logger))
		time.sleep(1)
		pool.apply_async(run_basedepth,args=(finalbam, targetbed, genomefile, outfolder, sample, logger))
		time.sleep(1)
		outfile = '%s/%s_targetReadsCover' % (outfolder,sample)
		pool.apply_async(run_regionreads,args=(finalbam, targetbed, outfile, logger))
		time.sleep(1)
		pool.close()
		pool.join()
		
		# generate QC table
		infile = '%s/%s_samtools_stats' % (outfolder, sample)
		reads,mapped,proper_mapped,error_rate,reads_len,insert_size = info_samstats(infile)	
		infile = markdupmetric
		duplicate_rate = float(open(infile).readlines()[7].strip().split('\t')[-2])
		fout = open('%s/%s_mapQC.txt'%(outfolder,sample),'w')
		print('QC Statistics Sample\t%s' % (sample), file=fout)
		print('Paired-end read length\t%d*2' % (reads_len), file=fout)
		print('Total effective data yield(Gb)\t%.2f' % (reads/1000000000.0*reads_len), file=fout)
		print('Total reads number(M)\t%.2f' % (reads/1000000.0), file=fout)
		print('Reads mapping rate\t%.2f%%' % (mapped*1.0/reads*100), file=fout)
		print('Properly paired mapping reads rate\t%.2f%%' % (proper_mapped*1.0/reads*100), file=fout)
		print('Average insert size\t%.2f' % (insert_size), file=fout)
		print('No-mismatch mapping reads rate:\t%.2f%%' % (read_NOmismatch*1.0/reads*100), file=fout)
		print('Mismatch alignment bases rate:\t%.2f%%' % (error_rate*100), file=fout)
		infile='%s/%s_targetReadsCover' % (outfolder,sample)
		reads0 = int(open(infile).readline().strip())
		print('Capture efficiency rate on target regions\t%.2f%%' % (reads0*100.0/reads), file=fout)
		infile = '%s/%s_base_depth.txt' % (outfolder,sample)
		covers, map_Base, total_Base = summary_depth(infile)
		print('PCR duplication rate\t%.2f%%' % (duplicate_rate*100), file=fout)
		print('Mean coverage sequencing depth on official target\t%.2f'%(map_Base*1.0/total_Base), file=fout)
		print('Fraction of official target covered\t%.2f%%' % (covers[0]*100), file=fout)
		print('Fraction of official target covered with at least 4X\t%.2f%%' % (covers[1]*100), file=fout)
		print('Fraction of official target covered with at least 10X\t%.2f%%' % (covers[2]*100), file=fout)
		print('Fraction of official target covered with at least 20X\t%.2f%%' % (covers[3]*100), file=fout)
		fout.close()
		os.mknod("%s/%s_mapQC.success" % (outfolder, sample))
	except:
		os.mknod("%s/%s_mapQC.fail" % (outfolder, sample))

def main():
	parser=argparse.ArgumentParser()
	config_default_file=os.path.join(os.path.dirname(os.path.realpath(__file__)),"wes_config.ini")
	parser.add_argument('--finalbam',help='The final bam file')
	parser.add_argument('--markdupmetric',help='The marked duplicate metric')
	parser.add_argument('--sample',help='The sample name')
	parser.add_argument('--output',help='The result folder path')
	parser.add_argument('--config',default=config_default_file,help='The configure file. Default: %s'%(config_default_file))
	args=parser.parse_args()

	if (args.finalbam is None  or args.markdupmetric is None
	    or args.sample is None or args.output is None):
		parser.print_usage()
		sys.exit(1)
	args.output = args.output.rstrip('/')
	logger = cmd_logger.create_logger(args.sample, "%s/%s_mapQC.log" % (args.output,args.sample))
	run_qc(args.finalbam, args.markdupmetric, args.sample,
	       args.output, args.config, logger)

if __name__=="__main__":
	main()
