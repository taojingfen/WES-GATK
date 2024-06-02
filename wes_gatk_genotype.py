#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   wes_gatk_VQSR.py
@Time    :   2020/10/20 13:31:04
@Author  :   tao_jingfen
@Contact :   taojingfen@foxmail.com
@Desc    :   Filter Variants by Variant (Quality Score) Recalibration 
'''

import os
import sys
import argparse
import load_config
import cmd_logger

def snp_recal(invcf: str, sample: str, outfolder: str,
              config: str, logger: str):
    """run SNP VQSR

    - Args: 
        - invcf: the vcf file after variant calling,
          e.g. {sample}.vcf.gz
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - config: the config file with ref and software information
          e.g. ./wes_config.ini
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - outvcf: the snp vcf file after filtering, 
          e.g. {sample}.VQSR.snp.vcf
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    hapmap = para_database['hapmap']
    omni = para_database['omni']
    db1000g_phase1_snp = para_database['db1000g_phase1_snp']
    dbsnp = para_database['dbsnp']
    cmd_dict = {}
    for key in ['invcf', 'ref', 'hapmap', 'omni', 'dbsnp', 
                'db1000g_phase1_snp', 'sample','outfolder']:
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
        os.mknod("{outfolder}/{sample}_VQSR_SNP.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished SNP gatk SelectVariants")
    snp_recal_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                    "VariantRecalibrator "
                    "-V {outfolder}/{sample}.snp.vcf "
                    "-R {ref} "
                    "-mode SNP "
                    "--resource:hapmap,known=false,training=true,"
                    "truth=true,prior=15.0 {hapmap} "
                    "--resource:omni,known=false,training=true,"
                    "truth=false,prior=12.0 {omni} "
                    "--resource:1000G,known=false,training=true,"
                    "truth=false,prior=10.0 {db1000g_phase1_snp} "
                    "--resource:dbsnp,known=true,training=false,"
                    "truth=false,prior=2.0 {dbsnp} "
                    "-an QD -an MQ -an MQRankSum -an ReadPosRankSum "
                    "-an FS -an SOR --max-gaussians 4 "
                    "--tranches-file {outfolder}/{sample}.VQSR.snp.tranches "
                    "-O {outfolder}/{sample}.VQSR.snp.recal"
                    ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(snp_recal_cmd, logger)
    if not tag:
        logger.error("Failed SNP gatk VariantRecalibrator")
        os.mknod("{outfolder}/{sample}_VQSR_SNP.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished SNP gatk VariantRecalibrator")
    snp_apply_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                    "ApplyVQSR "
                    "-V {outfolder}/{sample}.snp.vcf "
                    "-R {ref} "
                    "-mode SNP "
                    "--recal-file {outfolder}/{sample}.VQSR.snp.recal "
                    "-ts-filter-level 99.0 "
                    "--tranches-file {outfolder}/{sample}.VQSR.snp.tranches "                   
                    "-O {outfolder}/{sample}.VQSR.snp.vcf"
                    ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(snp_apply_cmd, logger)
    if not tag:
        logger.error("Failed SNP gatk ApplyVQSR")
        os.mknod("{outfolder}/{sample}_VQSR_SNP.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished SNP gatk ApplyVQSR")
        os.mknod("{outfolder}/{sample}_VQSR_SNP.success".format(**cmd_dict))  
    outvcf = "{outfolder}/{sample}.VQSR.snp.vcf".format(**cmd_dict)
    return tag, outvcf

def indel_recal(invcf: str, sample: str, outfolder: str,
                config: str, logger: str):
    """run Indel VQSR

    - Args: 
        - invcf: the vcf file after variant calling,
          e.g. {sample}.vcf.gz
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - config: the config file with ref and software information
          e.g. ./wes_config.ini
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - outvcf: the Indel vcf file after filtering, 
          e.g. {sample}.VQSR.indel.vcf
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    dbsnp = para_database['dbsnp']
    mills_indel = para_database['mills_indel']
    cmd_dict = {}
    for key in ['invcf', 'ref', 'dbsnp', 
                'mills_indel', 'sample','outfolder']:
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
        os.mknod("{outfolder}/{sample}_VQSR_INDEL.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished INDEL gatk SelectVariants")
    indel_recal_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                    "VariantRecalibrator "
                    "-V {outfolder}/{sample}.indel.vcf "
                    "-R {ref} "
                    "-mode INDEL "
                    "--resource:mills,known=false,training=true,"
                    "truth=true,prior=12.0 {mills_indel} "
                    "--resource:dbsnp,known=true,training=false,"
                    "truth=false,prior=2.0 {dbsnp} "
                    "-an QD -an DP -an MQRankSum -an ReadPosRankSum "
                    "-an FS -an SOR -an InbreedingCoeff --max-gaussians 4 "
                    "--tranches-file {outfolder}/{sample}.VQSR.indel.tranches "
                    "-O {outfolder}/{sample}.VQSR.indel.recal"
                    ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(indel_recal_cmd, logger)
    if not tag:
        logger.error("Failed INDEL gatk VariantRecalibrator")
        os.mknod("{outfolder}/{sample}_VQSR_INDEL.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished INDEL gatk VariantRecalibrator")
    indel_apply_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                    "ApplyVQSR "
                    "-V {outfolder}/{sample}.indel.vcf "
                    "-R {ref} "
                    "-mode INDEL "
                    "--recal-file {outfolder}/{sample}.VQSR.indel.recal "
                    "-ts-filter-level 99.0 "
                    "--tranches-file {outfolder}/{sample}.VQSR.inde.tranches "                   
                    "-O {outfolder}/{sample}.VQSR.indel.vcf"
                    ).format(**cmd_dict)
    tag = cmd_logger.run_cmd(indel_apply_cmd, logger)
    if not tag:
        logger.error("Failed INDEL gatk ApplyVQSR")
        os.mknod("{outfolder}/{sample}_VQSR_INDEL.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished INDEL gatk ApplyVQSR")
        os.mknod("{outfolder}/{sample}_VQSR_INDEL.success".format(**cmd_dict))  
    outvcf = "{outfolder}/{sample}.VQSR.indel.vcf".format(**cmd_dict)
    return tag, outvcf

def merge_vcf(invcf: str, snpvcf: str, indelvcf: str, sample: str,
              outfolder: str, config: str, logger: str):
    """merge the vcffile

    - Args: 
        - invcf: the vcf file after variant calling,
          e.g. {sample}.vcf.gz
        - snpvcf: the snp vcf file after VQSR
          e.g. {sample}.VQSR.snp.vcf
        - indelvcf: the indel vcf file after VQSR
          e.g. {sample}.VQSR.indel.vcf
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - outvcf: the final vcf file after VQSR, 
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
    logger = cmd_logger.create_logger(args.sample, "%s/%s_VQSR.log" % (
                    args.output, args.sample))
    tag, snpvcf = snp_recal(args.invcf, args.sample, args.output.rstrip('/'),
                          args.config,logger)
    tag, indelvcf = indel_recal(args.invcf, args.sample, args.output.rstrip('/'),
                          args.config,logger)
    tag, outvcf = merge_vcf(args.invcf, snpvcf, indelvcf, args.sample, 
                  args.output.rstrip('/'), args.config, logger)

if __name__=="__main__":
    main()

