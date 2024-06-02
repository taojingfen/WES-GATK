#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   wes_gatk_hcaller.py
@Time    :   2020/10/16 15:33:30
@Author  :   tao_jingfen
@Contact :   zb-taojingfen@kingmed.com.cn
@Desc    :   Do the variant calling use HaplotypeCaller and GenotypeGVCFs
'''

import os
import sys
import argparse
import load_config
import cmd_logger

def variant_haplotypecaller(inbam: str, sample: str, outfolder: str,
                 config: str, logger: str):
    """run HaplotypeCaller

    - Args: 
        - inbam: the bamfile after BQSR,
          e.g. {sample}.final.bam
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - config: the config file with ref and software information
          e.g. ./wes_config.ini
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - outgvcf: the gvcffile after calling, 
          e.g. {sample}.g.vcf.gz
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    bedfile = para_database['bedfile']
    cmd_dict = {}
    for key in ['inbam', 'ref', 'bedfile', 
                'sample','outfolder']:
        cmd_dict[key] = eval(key)
    hcaller_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                 "HaplotypeCaller " 
                 "-I {inbam} "
                 "-R {ref} "
                 "-L {bedfile} "
                 "-ERC GVCF "
                 "-O {outfolder}/{sample}.g.vcf.gz").format(**cmd_dict)
    tag = cmd_logger.run_cmd(hcaller_cmd, logger)
    if not tag:
        logger.error("Failed gatk HaplotypeCaller")
        os.mknod("{outfolder}/{sample}_HaplotypeCaller.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished gatk HaplotypeCaller")
        os.mknod("{outfolder}/{sample}_HaplotypeCaller.success".format(**cmd_dict))  
    outgvcf = "{outfolder}/{sample}.g.vcf.gz".format(**cmd_dict)
    return tag, outgvcf

def variant_genotype(gvcf: str, sample: str, outfolder: str,
                 config: str, logger: str):
    """run GenotypeGVCFs

    - Args: 
        - gvcf: the gvcffile that remove duplicate,
          e.g. {sample}.g.vcf.gz
        - sample: the sample name, e.g. NA12878
        - outfolder: the output path, e.g. ./result
        - config: the config file with ref and software information
          e.g. ./wes_config.ini
        - logger: instance of logging
    - Returns:
        - tag: True success else False
        - outvcf: the vcffile after genotyping, 
          e.g. {sample}.raw.vcf
    """
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    para_database = load_config.get_config(config,'database')
    ref = para_database['ref']
    bedfile = para_database['bedfile']
    dbsnp = para_database['dbsnp']
    cmd_dict = {}
    for key in ['gvcf', 'ref', 'bedfile', 'dbsnp', 
                'sample','outfolder']:
        cmd_dict[key] = eval(key)
    genotype_cmd = ("gatk --java-options \"-Xmx20g -Djava.io.tmpdir={outfolder}\" "
                   "GenotypeGVCFs "
                   "-V {gvcf} "
                   "-R {ref} "
                   "--dbsnp {dbsnp} "
                   "-O {outfolder}/{sample}.raw.vcf").format(**cmd_dict)
    tag = cmd_logger.run_cmd(genotype_cmd, logger)
    if not tag:
        logger.error("Failed gatk GenotypeGVCFs")
        os.mknod("{outfolder}/{sample}_GenotypeGVCF.fail".format(**cmd_dict))
        sys.exit(1)
    else:
        logger.info("Finished gatk GenotypeGVCFs")
        os.mknod("{outfolder}/{sample}_GenotypeGVCF.success".format(**cmd_dict))  
    outvcf = "{outfolder}/{sample}.raw.vcf".format(**cmd_dict)
    return tag, outvcf

def main():
    parser = argparse.ArgumentParser()
    config_default_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),"wes_config.ini")
    parser.add_argument('--inbam',help = 'The BQSR bam file')
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
    logger = cmd_logger.create_logger(args.sample, "%s/%s_hcaller.log" % (
                    args.output, args.sample))
    tag, outgvcf = variant_haplotypecaller(args.inbam, args.sample,
                   args.output.rstrip('/'), args.config, logger)
    tag, outvcf = variant_genotype(outgvcf, args.sample, 
                  args.output.rstrip('/'), args.config, logger)

if __name__=="__main__":
    main()
