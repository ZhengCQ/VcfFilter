# VcfFilter
#### 过滤低深度或高深度，标记群体与参考序列间的homogeneous
#### 适用
经GATK合并过滤后的vcf或者snpEff注释后的vcf

#### 示范 
```
python ../bin/Filter.py -v combine_cohort.filtered_pass.snpIndel.vcf.gz -groups A:B -samples R1095636,R17172210:R17172212,R180739,R180741 -dp 4,1000
```
1. 针对每一组group进行homogenous标记,要求这一组的基因型一致并且为突变基因型。
2. 深度不满足要求的，会被标记为'./.'miss基因型。
3. 标记后的vcf根据实际需求过滤。或直接取全部PASS
```
egrep '#|PASS' out.vcf >filter.vcf
```