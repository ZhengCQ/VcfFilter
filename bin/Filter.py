#!/usr/bin/env python
#coding=utf-8"

import sys
import argparse
import Read as Read
import time

ARGS = argparse.ArgumentParser(description="VCF过滤")
ARGS.add_argument(
	'-v', '--vcf', dest='vcf', required=True, help='vcf')
ARGS.add_argument(
	'-w', '--work_dir', dest='work_dir', default='.', help='工作目录，默认：当前目录。')
ARGS.add_argument(
	'-o', '--output', dest='output_pre', default='recode_filter', help='输入文件名前缀, 默认: recode_filter')
ARGS.add_argument(
	'-groups',  dest='group_name', 
	help='Groups的名称，多个样本以:分隔，如group1:group2')
ARGS.add_argument(
	'-samples',  dest='samples_name', 
	help='具体样本名称，多个样本以逗号分隔，如sample1,sample2.多个组别用:分隔,\
	如group1_sample1,group1_sample2:group2_sample1,group2_sample2')

class HandleGroup(object):
	"""
	@return grp_dict: {'A': {'index': [0, 3, 4, 6], 'name': ['AB00001802', 'GCT_4libs', 'GCT_AB0000353', 'GCT_AB0001103']}, 
	'B': {'index': [1, 5, 6, 7], 'name': ['AB0003001M1', 'GCT_AB0001004', 'GCT_AB0001103', 'GCT_AB0002974']}}
	"""
	def __init__(self, invcf, grp_dict, grp_name="", samples_lst=[]):
		self.invcf = invcf
		self.all_samples_lst = Read.Readvcf(self.invcf).samples.split('|')
		self.grp_dict = grp_dict
		if grp_name:
			self.grp_name = grp_name
		else:
			self.grp_name = 'All'
		if len(samples_lst) > 1:
			self.samples_lst = samples_lst
		else:
			self.samples_lst = self.all_samples_lst
		self.handle_grp()

	def handle_grp(self):
	#将存在于vcf中的样本名称和idx存于字典
		for idx, val in enumerate(self.all_samples_lst):
			if val in self.samples_lst:
				self.grp_dict.setdefault(self.grp_name, {}).setdefault('name', []).append(val)
				self.grp_dict.setdefault(self.grp_name, {}).setdefault('index', []).append(idx)

		#检查样本在vcf是否存在
		for each in self.samples_lst:
			try:
				if each not in self.grp_dict[self.grp_name]['name']:
					print ('Warnings: %s not in vcf'%(each))
			except:
				print ('Warnings: %s haven\'t sample in vcf'%(self.grp_name))

class GenotypeCount(object):
	def __init__(self, gtinfo, grp_dict, grp_name=""):
		self.gtinfo = gtinfo
		self.grp_dict = grp_dict
		self.grp_name_lst = []
		if grp_name:
			if grp_name in grp_dict:
				self.grp_name_lst.append(grp_name)
			else:
				print ("Error {} not in {}".format(grp_name, grp_dict))
		else:
			self.grp_name_lst = grp_dict.keys()
		self.gt_mut = {}
		self.main()

	def add_num(self, grp_name, sample_name, mut_type):
		#初始化每个样本及每个突变基因型的值，
		for i in ['miss', 'ref_hom', 'mut_hom', 'mut_mut', 'mut_het']:
			self.gt_mut.setdefault(sample_name, {}).setdefault(i, 0)
			try:
				self.gt_mut[grp_name][i]
			except:
				self.gt_mut.setdefault(grp_name, {}).setdefault(i, 0)
		#指定类型的值增加
		self.gt_mut[grp_name][mut_type] += 1
		self.gt_mut[sample_name][mut_type] += 1

	def group_count(self, grp_name):
		grp_gt_list = [self.gtinfo.split('|')[i] for i in self.grp_dict[grp_name]['index']]

		for idx, val in enumerate(grp_gt_list):
			sample_name = self.grp_dict[grp_name]['name'][idx]
			gt = val.split('/')
			if gt[0] == '.' and gt[1] == '.':
				self.add_num(grp_name, sample_name, 'miss') #./.				
			elif gt[0] == '0' and gt[1] == '0':#字符为0,非数字
				self.add_num(grp_name, sample_name, 'ref_hom') #0/0
			elif gt[0] != '0' and gt[1] != '0':
				if gt[0] == gt[1]:
					self.add_num(grp_name, sample_name, 'mut_hom') #1/1
				else:
					add_num(grp_name, sample_name, 'mut_mut') #类似突变为1/2
			elif gt[0] != '0' or gt[1] != '0':
				self.add_num(grp_name, sample_name, 'mut_het') #0/1，1/0, 2/0等
	
	def main(self):
		for grp_name in self.grp_name_lst:
			self.group_count(grp_name)


class Filter(object):
	def __init__(self, eachinfo, grp_dict, lowdepth = 4, highdepth=1000):
		self.__each = eachinfo
		self.__grp_dict = grp_dict
		self.__lowdpth = int(lowdepth)
		self.__highdepth = int(highdepth)
		self.__infos = eachinfo.infos
		self.line = ''
		self.dp_num = 0
		self.homogeneous_num = 0
		self.run()

	def depth(self, gtdetail):
		"""
		Depth少于或者大于阈值，该基因型将变为缺失
		"""
		mark_num = 0
		for idx, val in enumerate(self.__each.dp.split('|')):
			if int(val) <= self.__lowdpth or int(val) >= self.__highdepth:
				tmp = gtdetail[idx].split(':')[1:]
				tmp.insert(0, './.')
				gtdetail[idx] = ':'.join(tmp)
				mark_num = mark_num + 1
		return gtdetail, mark_num

	def homo_mark(self, group_gt, group_name):
		"""
		该group中的基因型一致，且不为ref，添加标记
		"""
		filter_mark = ''
		if len(set(group_gt)) == 1 and len(set(['0|0','0/0']) & set(group_gt)) == 0:
			filter_mark = 'GROUP {} homogeneous mutation'.format(group_name)
		return filter_mark

	def homogeneous(self, gtdetail):
		"""
		添加每个group的标记
		"""
		newdetail = []
		mark_num = 0
		for group_name in self.__grp_dict:
			group_gt = []
			for i in self.__grp_dict[group_name]['index']:
				newdetail.append(gtdetail[i])
				group_gt.append(gtdetail[i].split(':')[0])
			filter_mark = self.homo_mark(group_gt, group_name)
			if filter_mark:
				self.__infos[6] = self.__infos[6].replace('PASS', '')
				self.__infos[6] = self.__infos[6] + filter_mark + ';'
				mark_num = mark_num + 1
		return newdetail, mark_num

	def run(self):
		#低/高深度标记
		newgtdetail, mark_num = self.depth(self.__infos[9:])
		self.dp_num = mark_num
		self.__infos[9:] = newgtdetail
        #相同基因型标记
		newgtdetail, mark_num = self.homogeneous(self.__infos[9:])
		self.homogeneous_num = mark_num
		self.__infos[9:] = newgtdetail
		self.line = '\t'.join(self.__infos)


def main():
	args = ARGS.parse_args()
	invcf = args.vcf
	outvcf = invcf.replace('.vcf', '').replace('.vcf.gz', '') + args.output_pre + '.vcf'
	outfile = open(outvcf, 'w')

	print """Input vcf file: {}
Output vcf file: {}/{}""".format(invcf, args.work_dir, outvcf)

	if not invcf:
		print('Use --help for command line help')
		return
	try:
		os.makedirs(args.work_dir)
	except:
		pass
		#print ('%s exists' %(args.work_dir)) 
    
    #处理gourp
	grp_dict = {}
	if args.group_name and args.samples_name:	
		groups = args.group_name
		samples = args.samples_name
		for idx, val in enumerate(groups.split(':')):
			HandleGroup(invcf, grp_dict, val, samples.split(':')[idx].split(','))
	else:
		"""
		{'All': {'index': [0, 1, 2, 3], 'name': ['AB00001802', 'GCT_4libs', 'GCT_AB0000353', 'GCT_AB0001103']},
		"""
		HandleGroup(invcf, grp_dict)

	##header
	outfile.write("{}".format(Read.Readvcf(invcf).header))

	#samples行
	samples = Read.Readvcf(invcf).samples.split('|')
	newsamples = []
	for group_name in grp_dict:
		for i in grp_dict[group_name]['index']:
			newsamples.append(samples[i])
	outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(newsamples)))

    #具体数值
	vcfinfo = Read.Readvcf(invcf).extract
	total_dp_num = 0
	total_homogeneous_num = 0
	for each in vcfinfo:
		neweach = Filter(each, grp_dict)
		total_dp_num = total_dp_num + neweach.dp_num
		total_homogeneous_num = total_homogeneous_num + neweach.homogeneous_num
		outfile.write("{}\n".format(neweach.line))
	print "低/高深度标记次数: {}".format(total_dp_num)
	print "Homogeneous标记次数: {}".format(total_homogeneous_num)

if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print "时间总计{:0.2f}S".format(end - start)
	
