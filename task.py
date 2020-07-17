#-*- coding: UTF-8 -*-
import numpy as py
import matplotlib
import matplotlib.pyplot as plt
import sys
import argparse as arg

# 1.1 读取fasta文件(将fasta文件的每一条序列以字典的形式储存起来，字典的键为fasta序列的名称，值为其碱基序列)
def read_fa(path):
	with open(path , "r") as f:
		seqs = {}
		for line in f:
			if line.startswith(">"):#判断是否以>开头
				name = line.strip("\n")#是，去掉换行符
				seqs[name] = ""  # 定义一个字符串，以字符串形式存放元素
			else:
				seqs[name]+=line.replace("\n","")#去掉行尾的换行符(\n)等后按行存入seqs
	f.close()
	return seqs


# 1.2 fuzzy matching/模糊匹配
barcode1="GTACACGCTGTGACTA"
rev_barcode1="TAGTCACAGCGTGTAC"
barcode2="TCTATGTCTCAGTAGT"
rev_barcode2="ACTACTGAGACATAGA"

sam1 = open(r"C:\Users\hp\Desktop\test-master\sample1.fa", "w", encoding="utf-8")
sam2 = open(r"C:\Users\hp\Desktop\test-master\sample2.fa", "w", encoding="utf-8")
sam3 = open(r"C:\Users\hp\Desktop\test-master\sample3.fa", "w", encoding="utf-8")

#自定义mismatch函数进行评分
def mismatch(x,y):
    xscore = 0
    for i in range(len(x)):
        if x[i] == y[i]:
            xscore += 1
        else:
            xscore += 0
    return(xscore)

#根据得分对序列进行划分并输出到相应fasta文件中
a = read_fa(r"C:\Users\hp\Desktop\test-master\reads_of_insert.fa")
for key,val in a.items():
	if mismatch(val[:16],barcode1) >= 15 or mismatch(val[:16],rev_barcode1) >= 15:
		sam1.write(key+"\n"+val[16:]+"\n")
	elif mismatch(val[:16],barcode2) >= 15 or mismatch(val[:16],rev_barcode2) >= 15:
		sam2.write(key+"\n"+val[16:]+"\n")
	elif mismatch(val[-16:],barcode1) >= 15 or mismatch(val[-16:],rev_barcode1) >= 15:
		sam1.write(key+"\n"+val[:-16]+"\n")
	elif mismatch(val[-16:],barcode2) >= 15 or mismatch(val[-16:],rev_barcode2) >= 15:
		sam2.write(key+"\n"+val[:-16]+"\n")
	else:
		sam3.write(key+"\n"+val+"\n")

sam1.close()
sam2.close()
sam3.close()#不关闭上一次的文件始终打开，

b=read_fa(r"C:\Users\hp\Desktop\test-master\sample1.fa")
c=read_fa(r"C:\Users\hp\Desktop\test-master\sample2.fa")

#按一定序列长度（70bp）分割序列，按行排列后输出新的fasta文件
def length_limit(path,dict,length):
	newsam1 = open(path,"w")
	for key,val in dict.items():
		print(key,file=newsam1)
		while len(val)>length:
			print(val[0:length],file=newsam1)
			val=val[length:len(val)]
		print(val,file=newsam1)
	newsam1.close()
length_limit(r"C:\Users\hp\Desktop\test-master\newsample1.fa",b,70)
length_limit(r"C:\Users\hp\Desktop\test-master\newsample2.fa",c,70)



#2.去冗余序列并统计
# 2.1 先对字典进行两次键值反转去除重复序列
def remove_dup(dict,path):
	dict1 = {dict[key]:key for key in dict}
	dict2 = {dict1[key]:key for key in dict1}
	a1 = []
	for i in dict2.keys():
		a1.append(i)	# 将新字典b1的key存入列表a1

# 2.2 对原字典value进行冗余数的统计，并通过zip成元组列表去冗余
	b1 = []
	c1 = []
	for i in dict.values():
		b1.append(i)
	for i in b1:
		c1.append(b1.count(i) - 1)
	d1 = list(map(str, c1))  # 将整型转化为字符串形式
	_merge1 = zip(d1, b1)  # 生成一个打包成元组的列表
	# 去除重复元素
	_unique = []
	for (j, k) in _merge1:
		if (j, k) not in _unique:
			_unique.append((j, k))

# 2.3 将去重后的序列ID与冗余数，序列相联系，并根据冗余数进行排序
	_merge2 = zip(a1,_unique)
	h1 = sorted(_merge2, key=lambda i: (int(i[1][0]), i[0]), reverse=True)
# 2.4 导出文件
	sam_uniq = open(path, "w")
	for i, tup in h1:
		sam_uniq.write(i + "|" + 'duplicate number is' + tup[0] + "\n" + tup[1] + "\n")

remove_dup(b,r"C:\Users\hp\Desktop\test-master\sample1.uniq.fa")
remove_dup(c,r"C:\Users\hp\Desktop\test-master\sample2.uniq.fa")

d=read_fa(r"C:\Users\hp\Desktop\test-master\sample1.uniq.fa")
e=read_fa(r"C:\Users\hp\Desktop\test-master\sample2.uniq.fa")

length_limit(r"C:\Users\hp\Desktop\test-master\newsample1.uniq.fa",d,70)
length_limit(r"C:\Users\hp\Desktop\test-master\newsample1.uniq.fa",e,70)



#3.绘制GC含量频数分布图
# 3.1 计算GC含量
#将fa文件中序列合并，并以一定窗口切割，最终得到切割后的序列列表
def cut1(dict,window):
    l=[]
    for i in dict.values():
        l.append(i)
    L="".join(l)
    return[L[i:i+window]for i in range(0,len(L),window)]
L1=cut1(d,100)
L2=cut1(d,500)
L4=cut1(e,100)
L5=cut1(e,500)

def cut2(dict):
	l=[]
	for i in dict.values():
		l.append(i)
	return l
L3=cut2(d)
L6=cut2(e)

#对列表中每段窗口求GC含量，并以列表形式输出
def GC_content(list):
	L2=[]
	for i in list:
		GC_count=i.count("C")+i.count("G")
		GC_content = round(GC_count / len(i), 2) * 1000 / 10
		L2.append(GC_content)
	return sorted(L2)

#绘图
def plt_hist(list,path):
	plt.hist(list,bins=40,histtype="bar")#bins默认为10
	plt.xlabel("GC content (%)")
	plt.ylabel("GC Frequency")
	plt.title("sample.$bin.GC_content distribution")
	plt.savefig(path)#保存图片
	return plt.show()

plt_hist(GC_content(L1),"C:/Users/hp/Desktop/test-master/sample1.1$bin.gc_content.png")
plt_hist(GC_content(L2),"C:/Users/hp/Desktop/test-master/sample1.2$bin.gc_content.png")
plt_hist(GC_content(L3),"C:/Users/hp/Desktop/test-master/sample1.3$bin.gc_content.png")
plt_hist(GC_content(L4),"C:/Users/hp/Desktop/test-master/sample2.1$bin.gc_content.png")
plt_hist(GC_content(L5),"C:/Users/hp/Desktop/test-master/sample2.2$bin.gc_content.png")
plt_hist(GC_content(L6),"C:/Users/hp/Desktop/test-master/sample2.3$bin.gc_content.png")