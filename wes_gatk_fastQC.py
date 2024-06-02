#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@File    :   wes_gatk_fastQC.py
@Time    :   2020/10/15 16:51:10
@Author  :   tao_jingfen
@Contact :   taojingfen@foxmail.com
@Desc    :   Run fastqc
'''

import os
import sys
import re
import glob
import argparse
import cmd_logger
import load_config
from multiprocessing.dummy import Pool as ThreadPool
from collections import Counter, OrderedDict

def run_fastqc(fastq, outfolder, logger): 
    fastqc_cmd = ("fastqc {fastq} -o {outfolder}"
                  ).format(fastq=fastq, outfolder=outfolder)
    tag = cmd_logger.run_cmd(fastqc_cmd, logger)

def run_trim(reads1, reads2, outfolder, prefix, truseq, logger):
    trim_cmd = ("trimmomatic PE "
                "-threads 20 "
                "-phred33 "
                "{reads1} "
                "{reads2} "
                "{outfolder}/{prefix}.R1.clean.fastq.gz "
                "{outfolder}/{prefix}.R1.unpair.fastq.gz " 
                "{outfolder}/{prefix}.R2.clean.fastq.gz "
                "{outfolder}/{prefix}.R2.unpair.fastq.gz "
                "ILLUMINACLIP:{truseq}:2:30:10:8:true "
                "LEADING:20 "
                "TRAILING:20 "
                "SLIDINGWINDOW:4:10 "
                "MINLEN:35;"
                ).format(reads1=reads1, reads2=reads2,
                outfolder=outfolder, prefix=prefix, truseq=truseq)
    tag = cmd_logger.run_cmd(trim_cmd, logger)

def parseFastqcHtml(htmlfile):
    with open(htmlfile) as f:
        s = ''.join(f.readlines())
        m1 = re.search(r'Total Sequences</td><td>(\d+)</td>', s)
        readCnt = int(m1.group(1))
        m2 = re.search(r'Sequence length</td><td>(?:\d+-)?(\d+)</td>', s)
        readLen = int(m2.group(1))
    return readCnt, readCnt * readLen

def fqStat(fastq, phred=33):
    Q20 = chr(20 + phred)
    Q30 = chr(30 + phred)
    with os.popen('zcat %s' % fastq) as fqFile:
        lineCnt = baseCnt = 0
        gCnt = cCnt = 0
        q20Cnt = q30Cnt = 0
        for line in fqFile:
            lineCnt += 1
            if lineCnt % 4 == 0:
                for qual in line.strip():
                    if qual >= Q30:
                        q30Cnt += 1
                        q20Cnt += 1
                    elif qual >= Q20:
                        q20Cnt += 1
            elif lineCnt % 2 == 0:
                seq = line.strip()
                baseCnt += len(seq)
                cnt = Counter(seq)
                gCnt += cnt['G']
                cCnt += cnt['C']
    readCnt = lineCnt / 4
    q20ratio = 1.0 * q20Cnt / baseCnt
    q30ratio = 1.0 * q30Cnt / baseCnt
    gcContent = 1.0 * (gCnt + cCnt) / baseCnt
    gcSkew = 1.0 * (gCnt - cCnt) / (gCnt + cCnt)
    return readCnt, baseCnt, q20ratio, q30ratio, gcContent, gcSkew

def cmbStat(outfolder, prefix, read='R1'):
    rawReadCnt, rawBaseCnt = parseFastqcHtml('%s/%s_%s_001_fastqc.html' % (outfolder, prefix, read))
    cleanReadCnt, cleanBaseCnt, q20ratio, q30ratio, gcContent, gcSkew = fqStat(
        '%s/%s.%s.clean.fastq.gz ' % (outfolder, prefix, read))
    return OrderedDict([
        ('sampleName',   "%s_%s" % (prefix, read)),
        ('rawReads',     rawReadCnt),
        ('rawBases',     rawBaseCnt),
        ('cleanReads',   cleanReadCnt),
        ('cleanBases',   cleanBaseCnt),
        ('cleanReadsRatio(%)',  round(100.0 * cleanReadCnt / rawReadCnt, 2)),
        ('cleanBasesRatio(%)',  round(100.0 * cleanBaseCnt / rawBaseCnt, 2)),
        ('Q20',          round(q20ratio, 4)),
        ('Q30',          round(q30ratio, 4)),
        ('GC-Content',   round(gcContent, 4)),
        ('GC-Skew',      round(gcSkew, 4)),
    ])

def qc_report(fqfolder: str, sample: str, lanes: str,
               outfolder: str, config: str, logger: str):
    para_database = load_config.get_config(config,'database')
    truseq = para_database['adapter']
    pool = ThreadPool(20)
    pool.starmap(run_fastqc, [(fq, outfolder, logger) for fq in 
       ['%s/%s_%s_%s_001.fastq.gz' % (fqfolder, sample, lane, R) 
       for lane in lanes for R in ['R1','R2']]])
    pool.close()
    pool.join()
    pool = ThreadPool(20)
    pool.starmap(run_trim,
      [('%s/%s_%s_R1_001.fastq.gz' % (fqfolder, sample, lane), 
        '%s/%s_%s_R2_001.fastq.gz' % (fqfolder, sample, lane),
        outfolder, 
        '%s_%s' % (sample, lane), 
        truseq, logger
        ) for lane in lanes])
    pool.close()
    pool.join()
    os.system("rm %s/*.unpair.fastq.gz" % outfolder)
    start = True
    with open('%s/%s_fastQC.txt' % (outfolder, sample), 'w') as out:
            for lane in lanes:
                for read in ['R1','R2']:
                    h = cmbStat(outfolder,'%s_%s' % (sample,lane),read)
                    if start:
                        out.write('\t'.join(h.keys()) + '\n')
                        start = False
                    out.write('\t'.join(map(str, h.values())) + '\n')
    os.mknod("%s/%s_fastQC.success" % (outfolder, sample))

def main():
    parser = argparse.ArgumentParser()
    config_default_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),"wes_config.ini")
    parser.add_argument('--fqfolder', type=str, help=("The fastq file folder," 
      "e.g. /data_obs/NGS_Data/Wet_Lab/201024_A01200_0048_BHTMKJDMXX/Project_UNDEFINED/tZTJCYCBhb-NP22FM0398/FASTQ"))
    parser.add_argument('--sample', type=str, help=("The prefix sample name,"
      "e.g. tZTJCYCBhb-NP22FM0398"))
    parser.add_argument('--output', type=str, help =("The result folder path, "
      "e.g. /home/taojingfen/projects/20200927_WES/test/tZTJCYCBhb-NP22FM0398"))
    parser.add_argument('--config',default = config_default_file,
        help = 'The configure file. Default: %s' % (config_default_file))
    args = parser.parse_args()

    if (args.fqfolder is None or args.sample is None 
        or args.output is None):
        parser.print_usage()
        sys.exit(1)
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    args.fqfolder = args.fqfolder.rstrip('/')
    args.output = args.output.rstrip('/')
    logger = cmd_logger.create_logger(args.sample, "%s/%s_fastQC.log" % (
             args.output, args.sample))
    lanes = [re.split(args.sample+'_(.+)_R1_001.fastq.gz',os.path.basename(i))[1] 
             for i in glob.glob(args.fqfolder + '/*_R1_001.fastq.gz')]
    qc_report(args.fqfolder, args.sample, lanes,
               args.output, args.config, logger)

if __name__ == '__main__':
    main()
