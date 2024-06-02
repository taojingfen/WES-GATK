#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   wes_gatk_BQSR.py
@Time    :   2020/10/15 16:50:10
@Author  :   tao_jingfen
@Contact :   zb-taojingfen@kingmed.com.cn
@Desc    :   Recalibration Base Quality Score
'''

import os
import sys
import argparse
import load_config
import cmd_logger

def bam_bqsr(inbam: str, sample: str, outfolder: str,
             config: str, logger: str):
    """run BQSR, Recalibrate base quality scores

    - Args: 
        - inbam: the bamfile that mark duplicate,
          e.g. {sample}.sort.merge.markdup.bam
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - config: the config file with ref and software information
          e.g. ./wes_config.ini
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - out_bam: the bamfile after BQSR, 
          e.g. {sample}.final.bam
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    dbsnp = para_database['dbsnp']
    mills_indel = para_database['mills_indel']
    db1000g_phase1_snp = para_database['db1000g_phase1_snp']
    cmd_dict = {}
    for key in ['inbam', 'sample', 'ref','dbsnp', 'mills_indel',
                'db1000g_phase1_snp', 'outfolder']:
        cmd_dict[key] = eval(key)
    baserecal_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
               "BaseRecalibrator "
               "-I {inbam} "
               "-R {ref} "
               "-O {outfolder}/{sample}_BQSR.table "
               "--known-sites {dbsnp} "
               "--known-sites {mills_indel} "
               "--known-sites {db1000g_phase1_snp} "
               ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(baserecal_cmd,logger)
    if not tag:
        logger.error("Failed gatk BaseRecalibrator ")
        os.mknod("{outfolder}/{sample}_BQSR.fail".format(**cmd_dict))
    else:
        logger.info("Finished gatk BaseRecalibrator ")
    bqsr_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                "ApplyBQSR "
                "-I {inbam} "
                "-R {ref} "
                "--bqsr-recal-file {outfolder}/{sample}_BQSR.table "
                "-O {outfolder}/{sample}.final.bam").format(**cmd_dict)
    tag = cmd_logger.run_cmd(bqsr_cmd,logger)
    if not tag:
        logger.error("Failed gatk ApplyBQSR ")
        os.mknod("{outfolder}/{sample}_BQSR.fail".format(**cmd_dict))
    else:
        logger.info("Finished gatk ApplyBQSR ")
        os.mknod("{outfolder}/{sample}_BQSR.success".format(**cmd_dict)) 
    out_bam = "{outfolder}/{sample}.final.bam".format(**cmd_dict)
    return tag, out_bam

def main():
    parser = argparse.ArgumentParser()
    config_default_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),"wes_config.ini")
    parser.add_argument('--inbam',help = 'The markduplicate bam file')
    parser.add_argument('--sample',help = 'The sample name')
    parser.add_argument('--output',help ='The result folder path')
    parser.add_argument('--config',default = config_default_file,
                        help = 'The configure file. Default: %s' % (config_default_file))
    args = parser.parse_args()

    if args.inbam is None is None or args.sample is None or args.output is None:
        parser.print_usage()
        sys.exit(1)
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    logger = cmd_logger.create_logger(args.sample, "%s/%s_BQSR.log" % (
                    args.output, args.sample))
    tag,out_bam = bam_bqsr(args.inbam, args.sample,args.output.rstrip('/'),
                           args.config,logger)

if __name__=="__main__":
    main()
