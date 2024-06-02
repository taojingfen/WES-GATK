#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   wes_gatk.py
@Time    :   2020/10/23 11:27:52
@Author  :   tao_jingfen
@Contact :   zb-taojingfen@kingmed.com.cn
@Desc    :   Do the WES SNP and INDEl analysis
             Based on GATK best practice
'''

import os
import sys
import argparse
import cmd_logger

import wes_gatk_mapping
import wes_gatk_markdup
import wes_gatk_BQSR
import wes_gatk_hcaller
import wes_gatk_hardFilter

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
	# 1. reads quality control
	fastQC_status_file = "%s/%s_fastQC.success" % (args.output, args.sample)
	qcscript = os.path.join(os.path.dirname(os.path.realpath(__file__)),"wes_gatk_fastQC.py")
	if not os.path.exists(fastQC_status_file):
		if os.system(("python %s "
		"--fqfolder %s --sample %s --output %s --config %s") % (
		qcscript, args.fqfolder, args.sample, args.output, args.config)):
			print("Fail fastQC")
			sys.exit(1)
	# 2. reads mapping and bam merge
	mapping_status_file = "%s/%s_merge.success" % (args.output, args.sample)
	if not os.path.exists(mapping_status_file):
		logger = cmd_logger.create_logger(args.sample+'_mapping', 
			 "%s/%s_mapping.log" % (args.output, args.sample))
		tag, sortmergebam = wes_gatk_mapping.map_merge(args.output, args.sample,
		                    args.output, args.config, logger)
		if not tag:
			print("Failed bwa map & sort & merge")
			sys.exit(1)
	else:
		sortmergebam = "%s/%s.sort.merge.bam" % (args.output, args.sample)
	# 3. bam mark duplicate
	markdup_status_file = "%s/%s_markdup.success" % (args.output, args.sample)
	if not os.path.exists(markdup_status_file):
		logger = cmd_logger.create_logger(args.sample+'_markdup', 
				"%s/%s_markdup.log" % (args.output, args.sample))
		tag,markdupbam,markdupmetric = wes_gatk_markdup.bam_markdup(sortmergebam, 
		            args.sample, args.output, logger)
		if not tag:
			print("Failed bam mark duplicate")
			sys.exit(1)
	else:
		markdupbam = "%s/%s.sort.merge.markdup.bam" % (args.output, args.sample)
		markdupmetric = "%s/%s_markdup_metrics.txt" % (args.output, args.sample)
	# 4. BQSR
	bqsr_status_file = "%s/%s_BQSR.success" % (args.output, args.sample)
	if not os.path.exists(bqsr_status_file):
		logger = cmd_logger.create_logger(args.sample+'_BQSR',
				 "%s/%s_BQSR.log" % (args.output, args.sample)) 		
		tag,finalbam = wes_gatk_BQSR.bam_bqsr(markdupbam, args.sample,
		              args.output, args.config, logger)
		if not tag:
			print("Failed BQSR")
			sys.exit(1)
	else:
		finalbam = "%s/%s.final.bam" % (args.output, args.sample)
	# 5. mapping quality
	mapQC_status_file = "%s/%s_mapQC.success" % (args.output, args.sample)
	qcscript = os.path.join(os.path.dirname(os.path.realpath(__file__)),"wes_gatk_mapQC.py")
	if not os.path.exists(mapQC_status_file):
		if os.system(("python %s "
		"--finalbam %s --markdupmetric %s --sample %s --output %s --config %s") % (
		qcscript, finalbam, markdupmetric, args.sample, args.output, args.config)):
			print("Fail mapQC")
			sys.exit(1)
	# 6. variant calling
	hcall_status_file = "%s/%s_HaplotypeCaller.success" % (args.output, args.sample)
	logger = cmd_logger.create_logger(args.sample+'_hcaller',
			 "%s/%s_hcaller.log" % (args.output, args.sample))
	if not os.path.exists(hcall_status_file):
		tag, outgvcf = wes_gatk_hcaller.variant_haplotypecaller(finalbam,
		               args.sample, args.output, args.config, logger)
		if not tag:
			print("Failed HaplotypeCaller")
			sys.exit(1)
	else:
		outgvcf = "%s/%s.g.vcf.gz" % (args.output, args.sample)	
	gtype_status_file = "%s/%s_GenotypeGVCF.success" % (args.output, args.sample)
	if not os.path.exists(gtype_status_file):
		tag, outvcf = wes_gatk_hcaller.variant_genotype(outgvcf, args.sample,
		              args.output, args.config, logger)
		if not tag:
			print("Failed GenotypeGVCFs")
			sys.exit(1)
	else:
		outvcf = "%s/%s.raw.vcf" % (args.output, args.sample)
	# 7. hardFilter
	snpvqsr_status_file = "%s/%s_SNP_hardFilter.success" % (
						  args.output, args.sample)
	logger = cmd_logger.create_logger(args.sample+'_hardFilter',
			 "%s/%s_hardFilter.log" % (args.output, args.sample))
	if not os.path.exists(snpvqsr_status_file):
		tag, snpvcf = wes_gatk_hardFilter.snp_recal(outvcf, args.sample,
		              args.output, args.config, logger)
		if not tag:
			print("Failed SNP hardFilter")
			sys.exit(1)
	else:
		snpvcf = "%s/%s.snp.filter.vcf" % (args.output, args.sample)
	indelvqsr_status_file = "%s/%s_INDEL_hardFilter.success" % (
						  args.output, args.sample)
	if not os.path.exists(indelvqsr_status_file):
		tag, indelvcf = wes_gatk_hardFilter.indel_recal(outvcf, args.sample,
		              args.output, args.config, logger)
		if not tag:
			print("Failed INDEL hardFilter")
			sys.exit(1)
	else:
		indelvcf = "%s/%s.indel.filter.vcf" % (args.output, args.sample)
	merge_status_file = "%s/%s_MergeVcfs.success" % (
						  args.output, args.sample)
	if not os.path.exists(merge_status_file):
		tag, finalvcf = wes_gatk_hardFilter.merge_vcf(outvcf, snpvcf, indelvcf, 
	                args.sample, args.output, args.config, logger)
		if tag:
			status_file = "%s/%s_wes_gatk.success" % (args.output,args.sample)
			if not os.path.exists(status_file):
				os.mknod('%s/%s_wes_gatk.success' % (args.output,args.sample))

if __name__=="__main__":
	main()
