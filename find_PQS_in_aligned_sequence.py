# -*- coding: utf-8 -*-
# """
# Created on 2022.12.12
#
# @author: Max
# """
#find PQS in aligned sequence
import os
import numpy as np
genome_dict={};AJ1102={}
basecount=0;gapcount=0
os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/2022.12/第一周/analysis-PQS-in-PEDV/")
with open("sequence.txt","w")as f:
    data = np.loadtxt('MK584552.1.csv', delimiter=",",usecols=[1,2,],skiprows=1)
    with open("PEDV_rawdata.aln.fasta","r")as F:
        for line in F:
            if line.startswith(">"):
                key = line.rstrip()
                genome_dict[key] = ""
            else:
                genome_dict[key] += line.rstrip()
    AJ1102={key: value for key, value in genome_dict.items() if key==">MK584552.1"}
    for key,value in AJ1102.items():
        for i in range(len(data)):
            #print(i)#0,1,..,14
            raw_start,raw_end=int(data[i][0]),int(data[i][1])
            #print(int(data[i][0]),int(data[i][1]))
            for base in value:
                basecount+=1
                if base=="-":
                    gapcount+=1
                if basecount-gapcount==raw_start:
                    real_start=basecount
                if basecount-gapcount==raw_end:
                    real_end=basecount
                    print(real_start, real_end)
                    f.write("PQS" + str(i+1)+"\n")
                    for key,value in genome_dict.items():
                        f.write(key+"\n")
                        f.write(value[real_start:real_end]+"\n")
                        print(value[real_start:real_end])
            basecount=0;gapcount=0