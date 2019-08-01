#!/usr/bin/env python
#coding=utf-8"
# Author: ChenQing Zheng
# Date: 2019.07.30

import sys
import os
import argparse
import time
import multiprocessing as mp

import Read as Read

ARGS = argparse.ArgumentParser(description="VCF过滤")
ARGS.add_argument(
	'-v', '--vcf', dest='vcf', required=True, help='vcf')
ARGS.add_argument(
	'-w', '--work_dir', dest='work_dir', default='.', help='工作目录，默认：当前目录。')
ARGS.add_argument(
	'-o', '--output', dest='output_pre', default='_recode_filter', help='输入文件名前缀, 默认: recode_filter')
ARGS.add_argument(
	'-groups',  dest='group_name', 
	help='Groups的名称，多个样本以:分隔，如group1:group2')
ARGS.add_argument(
	'-samples',  dest='samples_name', 
	help='具体样本名称，多个样本以逗号分隔，如sample1,sample2.多个组别用:分隔,\
	如group1_sample1,group1_sample2:group2_sample1,group2_sample2')
ARGS.add_argument(
	'-dp',  dest='depth_cut', default='4,1000',
	help='>=最低深度,<=最高深度,如:4,1000')
#ARGS.add_argument(
#    '--n_core', dest='n_core', type=int, default=1,
#    help='多进程数目, 默认为1')

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
			if int(val) < self.__lowdpth or int(val) > self.__highdepth:
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
		if len(set(group_gt)) == 1 and len(set(['0|0', '0/0', './.', '.|.']) & set(group_gt)) == 0:
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


class Run(object):
	def __init__(self, args):
		self.__invcf = args.vcf
		self.__outvcf = os.path.basename(self.__invcf).replace('.vcf', '').replace('.gz', '') + args.output_pre + '.vcf'
		self.__group_name = args.group_name
		self.__samples_name = args.samples_name
		self.__depth = args.depth_cut
		self.__total_dp_num = 0
		self.__total_homogeneous_num = 0
		self.__grp_dict = {}
		self.main()

	def run_filter(self, each, grp_dict, low_dp, high_dp):
		neweach = Filter(each, grp_dict, low_dp, high_dp)
		self.__total_dp_num = self.__total_dp_num + neweach.dp_num
		self.__total_homogeneous_num  = self.__total_homogeneous_num  + neweach.homogeneous_num
		self.wirte_file(neweach)
		#return neweach

	def wirte_file(self, each):
		with open(self.__outvcf, 'a+') as fi:
			fi.write('%s\n'%(each.line))

	def get_group(self):
		if self.__group_name and self.__samples_name:
			for idx, val in enumerate(self.__group_name.split(':')):
				HandleGroup(self.__invcf, self.__grp_dict, val, self.__samples_name.split(':')[idx].split(','))
		else:
			"""
			{'All': {'index': [0, 1, 2, 3], 'name': ['AB00001802', 'GCT_4libs', 'GCT_AB0000353', 'GCT_AB0001103']},
			"""
			HandleGroup(self.__invcf, self.__grp_dict)

	def write_head(self):
		with open(self.__outvcf, 'w') as fi:
			##header
			fi.write("{}".format(Read.Readvcf(self.__invcf).header))
			#samples行
			samples = Read.Readvcf(self.__invcf).samples.split('|')
			newsamples = []
			for group_name in self.__grp_dict:
				for i in self.__grp_dict[group_name]['index']:
					newsamples.append(samples[i])
			fi.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(newsamples)))

	def main(self):
		if not self.__invcf:
			print('Use --help for command line help')
			return
		try:
			os.makedirs(args.work_dir)
		except:
			pass
			#print ('%s exists' %(args.work_dir))

		self.get_group()
		self.write_head()
		low_dp, high_dp = self.__depth.split(',')
		vcfinfo = Read.Readvcf(self.__invcf).extract
		#pool = mp.Pool(1) #启动多线程池
		for each in vcfinfo:
			self.run_filter(each, self.__grp_dict, low_dp, high_dp)
			#pool.apply_async(self.run_filter, args=(each, self.__grp_dict, low_dp, high_dp))#函数写入到多线程池

		"""
		print('Waiting for all subprocesses done...')
		pool.close()
		pool.join()
		print('All subprocesses done.')
		pool.terminate()
		"""

		print ("低/高深度标记次数: {}".format(self.__total_dp_num))
		print ("Homogeneous标记次数: {}".format(self.__total_homogeneous_num))

if __name__ == '__main__':
	start = time.time()
	args = ARGS.parse_args()
	Run(args)
	end = time.time()
	print "时间总计{:0.2f}S".format(end - start)	
