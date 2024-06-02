#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   wes_gatk_mapping.py
@Time    :   2020/10/13 16:00:08
@Author  :   tao_jingfen
@Contact :   taojingfen@foxmail.com
@Desc    :   do the bwa mapping,
             map the fastq into reference
'''

import os
import sys
import argparse
import re
import glob
import load_config
import cmd_logger
from multiprocessing.dummy import Pool as ThreadPool


def bwa_mapping_sort(reads1: str, reads2: str, sample: str, lane: str,
                    outfolder: str, config: str, logger: str):
    """BWA mapping and alignment sorting, 
       do the <bwa mem> and <gatk SortSam>

    - Args：
        - reads1: the first fastq file, 
          e.g. /home/taojingfen/projects/20200927_WES/demo/input/sample1_1.fq
        - reads2: the second fastq file,
          e.g. /home/taojingfen/projects/20200927_WES/demo/input/sample1_2.fq
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, all the file will be placed in,
          e.g. ./result
        - config: the config file with ref and software information
          e.g. ./wes_config.ini
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - out_bam: the sort bamfile, e.g. {sample}.sort.bam
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    refdict = para_database['refdict']
    para_parameter = load_config.get_config(config,'parameter')
    cpu = para_parameter['cpu']
    cmd_dict = {}
    for key in ['sample','ref','cpu','reads1','reads2', 'outfolder', 'lane', 'refdict']:
        cmd_dict[key] = eval(key)
    mapping_success_file = "{outfolder}/{sample}_{lane}_mapping.success".format(**cmd_dict)
    if not os.path.exists(mapping_success_file):
        map_cmd = ("bwa mem -M -R \"@RG\\tID:{sample}\\tSM:{sample}\\tLB:{lane}\\tPL:ILLUMINA\" "
                "-t {cpu} {ref} {reads1} {reads2} | "
                "samtools view -Sb -o {outfolder}/{sample}_{lane}.raw.bam -"
                ).format(**cmd_dict)
        tag = cmd_logger.run_cmd(map_cmd,logger)
        if not tag:
            logger.error("Failed bwa mem")
            os.mknod("{outfolder}/{sample}_{lane}_mapping.fail".format(**cmd_dict))
            return tag
        else:
            logger.info("Finished bwa mem")
        sort_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                "SortSam -SO coordinate --VALIDATION_STRINGENCY SILENT "
                "-I {outfolder}/{sample}_{lane}.raw.bam "
                "-O {outfolder}/{sample}_{lane}.sort.bam "
                ).format(**cmd_dict)
        tag = cmd_logger.run_cmd(sort_cmd, logger)
        if not tag:
            logger.error("Failed gatk SortSam")
            os.mknod("{outfolder}/{sample}_{lane}_mapping.fail".format(**cmd_dict))
            return tag
        else:
            logger.info("Finished gatk SortSam")
        reorder_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                    "ReorderSam --VALIDATION_STRINGENCY SILENT "
                    "-I {outfolder}/{sample}_{lane}.sort.bam "
                    "-SD {refdict} "
                    "-O {outfolder}/{sample}_{lane}.reorder.bam"
                    ).format(**cmd_dict)
        tag = cmd_logger.run_cmd(reorder_cmd, logger)
        if not tag:
            logger.error("Failed gatk ReorderSam")
            os.mknod("{outfolder}/{sample}_{lane}_mapping.fail".format(**cmd_dict))
        else:
            logger.info("Finished gatk ReorderSam")
            os.system(("mv {outfolder}/{sample}_{lane}.reorder.bam "
                    "{outfolder}/{sample}_{lane}.sort.bam").format(**cmd_dict))
            os.system("rm {outfolder}/{sample}_{lane}.raw.bam".format(**cmd_dict))
            os.mknod("{outfolder}/{sample}_{lane}_mapping.success".format(**cmd_dict))
    else:
        tag = True
    return tag

def merge_bam(bamfolder: str, sample: str, lanes: str, 
              outfolder: str, config: str, logger: str):
    input = ' '.join(map(lambda x: '-I %s/%s_%s.sort.bam' % (bamfolder,sample,x),lanes))
    cmd_dict = {}
    for key in ['sample','input', 'outfolder']:
        cmd_dict[key] = eval(key)
    merge_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                 "MergeSamFiles "
                 "{input} "
                 "-O {outfolder}/{sample}.sort.merge.bam "
                ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(merge_cmd, logger)
    if not tag:
        logger.error("Failed gatk MergeSamFiles")
        os.mknod("{outfolder}/{sample}_merge.fail".format(**cmd_dict))
    else:
        logger.info("Finished gatk MergeSamFiles")
        os.mknod("{outfolder}/{sample}_merge.success".format(**cmd_dict))
    return tag

def map_merge(fqfolder: str, sample: str,
              outfolder: str, config: str, logger: str):
    lanes = [re.split(sample+'_(.+).R1.clean.fastq.gz',os.path.basename(i))[1]
             for i in glob.glob(fqfolder + '/*.R1.clean.fastq.gz')]
    pool = ThreadPool(20)
    result = pool.starmap(bwa_mapping_sort,
      [('%s/%s_%s.R1.clean.fastq.gz' % (fqfolder, sample, lane),
        '%s/%s_%s.R2.clean.fastq.gz' % (fqfolder, sample, lane),
        sample, 
        lane,
        outfolder,
        config,
        logger) for lane in lanes])
    pool.close()
    pool.join()
    if all(result):
        tag = merge_bam(outfolder, sample, lanes, outfolder, config, logger)
    else:
        tag = False
    out_bam = "%s/%s.sort.merge.bam" % (outfolder, sample)
    return tag, out_bam

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

    if args.fqfolder is None or args.sample is None or args.output is None:
        parser.print_usage()
        sys.exit(1)
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    args.fqfolder = args.fqfolder.rstrip('/')
    args.output = args.output.rstrip('/')
    logger = cmd_logger.create_logger(args.sample, "%s/%s_mapping.log" % (
             args.output, args.sample))
    tag, out_bam = map_merge(args.fqfolder, args.sample,
                   args.output, args.config, logger)
    if not tag:
        print("Failed bwa map & sort & merge")

if __name__=="__main__":
    main()
