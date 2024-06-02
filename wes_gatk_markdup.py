#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   wes_gatk_markdup_QC.py
@Time    :   2020/10/15 14:42:01
@Author  :   tao_jingfen
@Contact :   taojingfen@foxmail.com
@Desc    :   Mark the duplicate in the bamfile
'''

import os
import sys
import argparse
import cmd_logger

def bam_markdup(inbam: str, sample: str, outfolder: str, logger: str):
    """mark the duplicate in the bamfile

    - Args:
        - inbam: the merge bamfile, e.g. {sample}.sort.merge.bam
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - out_bam: the bamfile after mark duplicated, 
          e.g. {sample}.sort.merge.markdup.bam
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    cmd_dict = {}
    for key in ['inbam','sample','outfolder']:
        cmd_dict[key] = eval(key)
    markdup_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                 "MarkDuplicates "
                 "--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 " 
                 "-I {inbam} "
                 "-O {outfolder}/{sample}.sort.merge.markdup.bam "
                 "-M {outfolder}/{sample}_markdup_metrics.txt "
                 ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(markdup_cmd,logger)
    if not tag:
        logger.error("Failed gatk MarkDuplicates")
        os.mknod("{outfolder}/{sample}_markdup.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished gatk MarkDuplicates")
        os.mknod("{outfolder}/{sample}_markdup.success".format(**cmd_dict))  
    out_bam = "{outfolder}/{sample}.sort.merge.markdup.bam".format(**cmd_dict)
    out_metrics = "{outfolder}/{sample}_markdup_metrics.txt".format(**cmd_dict)
    return tag, out_bam, out_metrics

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inbam',help = 'The merge bam file')
    parser.add_argument('--sample',help = 'The sample name')
    parser.add_argument('--output',help ='The result folder path')
    args = parser.parse_args()

    if args.inbam is None is None or args.sample is None or args.output is None:
        parser.print_usage()
        sys.exit(1)
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    logger = cmd_logger.create_logger(args.sample, "%s/%s_markdup.log" % (
                    args.output, args.sample))
    tag,out_bam,out_metrics = bam_markdup(args.inbam, args.sample,
							args.output.rstrip('/'),logger)

if __name__=="__main__":
    main()
