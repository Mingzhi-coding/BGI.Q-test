简介：编程小白首次接触生信，花了两周时间算是刚刚入门python，经老师，师兄的指点才完成。在网上搜寻资料过程中发现有关资料较少，特此分享出来，欢迎大家批评指点。

题目要求：
拆分带有barcode标记的DNA序列（fasta格式），去除拆分后每个样本的冗余序列，并绘制去冗余序列后的GC content分布图。

具体要求：
1、拆分样本。已有fasta序列文件reads_of_insert.fa，其含有两个样本的序列。barcode（16个base）信息位于每条序列的开头或结尾。要求根据barcode将reads_of_insert.fasta中的两个样本序列拆分出来，2.存放在相应的文件中，结果文件命名sample1.fa、sample2.fa。
注：最终目标是得到去除barcode的两个与原来fasta文件类型一样的fasta文件。（若首尾同时存在barcode则同时去除）

2、去冗余序列。将第1步得到的两个结果文件中的冗余序列（长度、碱基顺序完全一致）去除（即冗余序列只保留一条），生成去冗余后的文件，ID行加入冗余数目信息，并将序列按照冗余数目倒序排列，结果文件命名sample1.uniq.fa、sample2.uniq.fa。
注：去冗余只是单个文件内的去冗余，即分别对第一步得到的两个fasta文件去冗余。还要求出重复次数并加在fasta第一行ID后。最后排序，导出到两个fasta文件中。

3、 GC content。 绘制第2步两个结果文件的GC content频数分布图（按100bp、500bp、整条序列三个bin），结果文件命名sample1. $bin.gc_content.png、sample2.$bin. gc_content.png。
横坐标：GC含量    纵坐标：GC含量频数
先滑动窗口获得GC含量获得列表，再对得到的列表进行count（估计此处应使用counter，先导入包），统计列表中GC含量出现的频数。

补充：
(1)   reads_of_insert.fa 路径：/zfsqd1/ST_OCEAN/USRS/shichch/training/Q1/

(2)   barcode 信息
sample1 GTACACGCTGTGACTA
sample2 TCTATGTCTCAGTAGT

（3） 待拆分序列的barcode，既可能位于待拆分序列开头、也可能位于待拆分序列结尾；即可能是提供的barcode序列，也可能是提供的barcode的反向互补序列。

（4）拆分容错：在barcode16个碱基的比较中，至多允许一个碱基错配。

1， 集群操作：熟悉Linux命令，程序必须全程在集群上操作，需要现场运行脚本跑出结果

2， 程序规范：

a)   程序不能报错

b)   程序需要有使用说明

c)   输入文件使用变量或者传参，不能写死

d)   输出结果使用变量或者传参输出到使用者规定的文件中，不要输出到屏幕或者写死

e)   请不要擅自改变固定的输出文件的格式（如：不建议将输出fasta/vcf/gff文件等命名为.txt）


其他细节问题可由技术导师们指点补充~
