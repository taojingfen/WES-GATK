#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   wes_gatk_hardFilter.py
@Time    :   2020/10/20 13:31:04
@Author  :   tao_jingfen
@Contact :   zb-taojingfen@kingmed.com.cn
@Desc    :   Filter Variants by Variant (Quality Score) Recalibration 
'''

import os
import sys
import argparse
import load_config
import cmd_logger

def snp_recal(invcf: str, sample: str, outfolder: str,
              config: str, logger: str):
    """run SNP hardFilter

    - Args: 
        - invcf: the vcf file after variant calling,
          e.g. {sample}.raw.vcf
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - config: the config file with ref and software information
          e.g. ./wes_config.ini
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - outvcf: the snp vcf file after filtering, 
          e.g. {sample}.snp.filter.vcf
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    cmd_dict = {}
    for key in ['invcf', 'ref', 'sample','outfolder']:
        cmd_dict[key] = eval(key)
    snp_select_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                      "SelectVariants "
                      "-V {invcf} "
                      "-R {ref} "
                      "--select-type-to-include SNP "
                      "-O {outfolder}/{sample}.snp.vcf"
                      ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(snp_select_cmd, logger)
    if not tag:
        logger.error("Failed SNP gatk SelectVariants")
        os.mknod("{outfolder}/{sample}_SNP_hardFilter.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished SNP gatk SelectVariants")
    snp_filter_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                      "VariantFiltration "
                      "-V {outfolder}/{sample}.snp.vcf "
                      "-R {ref} "
                      "-O {outfolder}/{sample}.snp.filter.vcf "
                      "-filter \"DP<=4 || QD<2.0 || MQ<40.0 || FS>60.0 || "
                      "HaplotypeScore>60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" " 
                      "--filter-name \"StandardFilterSNV\"").format(**cmd_dict)
    tag = cmd_logger.run_cmd(snp_filter_cmd, logger)
    if not tag:
        logger.error("Failed SNP gatk VariantFiltration")
        os.mknod("{outfolder}/{sample}_SNP_hardFilter.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished SNP gatk VariantFiltration")
        os.mknod("{outfolder}/{sample}_SNP_hardFilter.success".format(**cmd_dict))
    outvcf = "{outfolder}/{sample}.snp.filter.vcf".format(**cmd_dict)
    return tag, outvcf

def indel_recal(invcf: str, sample: str, outfolder: str,
                config: str, logger: str):
    """run Indel hardFilter

    - Args: 
        - invcf: the vcf file after variant calling,
          e.g. {sample}.raw.vcf
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - config: the config file with ref and software information
          e.g. ./wes_config.ini
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - outvcf: the Indel vcf file after filtering, 
          e.g. {sample}.indel.filter.vcf
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    cmd_dict = {}
    for key in ['invcf', 'ref', 'sample','outfolder']:
        cmd_dict[key] = eval(key)
    indel_select_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                    "SelectVariants "
                    "-V {invcf} "
                    "-R {ref} "
                    "--select-type-to-include INDEL "
                    "-O {outfolder}/{sample}.indel.vcf"
                    ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(indel_select_cmd, logger)
    if not tag:
        logger.error("Failed INDEL gatk SelectVariants")
        os.mknod("{outfolder}/{sample}_INDEL_hardFilter.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished INDEL gatk SelectVariants")
    indel_filter_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                      "VariantFiltration "
                      "-V {outfolder}/{sample}.indel.vcf "
                      "-R {ref} "
                      "-O {outfolder}/{sample}.indel.filter.vcf "
                      "-filter \"DP<=4 || QD<2.0 || FS>200.0 || ReadPosRankSum < -20.0\" "
                      "--filter-name \"StandardFilterINDEL\"").format(**cmd_dict)
    tag = cmd_logger.run_cmd(indel_filter_cmd, logger)
    if not tag:
        logger.error("Failed INDEL gatk VariantFiltration")
        os.mknod("{outfolder}/{sample}_INDEL_hardFilter.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished INDEL gatk VariantFiltration")
        os.mknod("{outfolder}/{sample}_INDEL_hardFilter.success".format(**cmd_dict))
    outvcf = "{outfolder}/{sample}.indel.filter.vcf".format(**cmd_dict)
    return tag, outvcf

def merge_vcf(invcf: str, snpvcf: str, indelvcf: str, sample: str,
              outfolder: str, config: str, logger: str):
    """merge the vcffile

    - Args: 
        - invcf: the vcf file after variant calling,
          e.g. {sample}.raw.vcf
        - snpvcf: the snp vcf file after hardFilter
          e.g. {sample}.snp.filter.vcf
        - indelvcf: the indel vcf file after hardFilter
          e.g. {sample}.indel.filter.vcf
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - outvcf: the Indel vcf file after filtering, 
          e.g. {sample}.final.vcf
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    cmd_dict = {}
    for key in ['invcf', 'ref', 'outfolder', 
                'sample', 'snpvcf', 'indelvcf']:
        cmd_dict[key] = eval(key)
    mix_select_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                      "SelectVariants "
                      "-V {invcf} "
                      "-R {ref} "
                      "--select-type-to-include MIXED "
                      "-O {outfolder}/{sample}.mix.vcf"
                     ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(mix_select_cmd, logger)
    if not tag:
        logger.error("Failed mix gatk SelectVariants")
        os.mknod("{outfolder}/{sample}_MergeVcfs.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished mix gatk SelectVariants")
    merge_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                 "MergeVcfs "
                 "-I {outfolder}/{sample}.mix.vcf "
                 "-I {snpvcf} "
                 "-I {indelvcf} "
                 "-O {outfolder}/{sample}.final.vcf "
                ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(merge_cmd, logger)
    if not tag:
        logger.error("Failed gatk MergeVcfs")
        os.mknod("{outfolder}/{sample}_MergeVcfs.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished gatk MergeVcfs")
        os.mknod("{outfolder}/{sample}_MergeVcfs.success".format(**cmd_dict))  
    outvcf = "{outfolder}/{sample}.final.vcf".format(**cmd_dict)
    return tag, outvcf

def main():
    parser = argparse.ArgumentParser()
    config_default_file = os.path.join(os.path.dirname(
                          os.path.realpath(__file__)),"wes_config.ini")
    parser.add_argument('--invcf',help = 'The vcf file after variant calling')
    parser.add_argument('--sample',help = 'The sample name')
    parser.add_argument('--output',help = 'The result folder path')
    parser.add_argument('--config',default = config_default_file,
          help = 'The configure file. Default: %s' % (config_default_file))
    args = parser.parse_args()

    if args.invcf is None is None or args.sample is None or args.output is None:
        parser.print_usage()
        sys.exit(1)
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    logger = cmd_logger.create_logger(args.sample, "%s/%s_hardFilter.log" % (
                    args.output, args.sample))
    tag, snpvcf = snp_recal(args.invcf, args.sample, args.output.rstrip('/'),
                          args.config,logger)
    tag, indelvcf = indel_recal(args.invcf, args.sample, args.output.rstrip('/'),
                          args.config,logger)
    tag, outvcf = merge_vcf(args.invcf, snpvcf, indelvcf, args.sample, 
                  args.output.rstrip('/'), args.config, logger)

if __name__=="__main__":
    main()

