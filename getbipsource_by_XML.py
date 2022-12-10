# -*- coding: utf-8 -*-
"""
Created on 2022.12.9

@author: Max
"""
#firststep-get_sequence_id
import os
genome_dict = {}
os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/2022.12/第一周/analysis-PQS-in-PEDV/raw_data/")
with open("PEDV_rawdata_id.txt","w")as f:
    with open("PEDV_rawdata.fasta","r")as F:
        for line in F:
            if line.startswith(">"):
                key=line[1:].rstrip()
                genome_dict[key] = ""
            else:
                genome_dict[key] += line
        for key, value in genome_dict.items():
            f.write(key+"\n")

#secondstep-get_biosource_by_XML
import os
import  xml.dom.minidom
dict_country={};dict_time={}
os.chdir("D:/360MoveData/Users/dgwei/Desktop/")
dom=xml.dom.minidom.parse("PEDV.gbc.xml")
root=dom.documentElement
seq=root.getElementsByTagName("INSDSeq")
for l in seq:
    id=l.getElementsByTagName("INSDSeq_accession-version")
    seqid=id[0].firstChild.data
    info=l.getElementsByTagName("INSDQualifier")
    for line in info:
        name=line.getElementsByTagName("INSDQualifier_name")
        list=name[0].firstChild.data
        if list=="country":
            key1=seqid
            dict_country[key1]=line.getElementsByTagName("INSDQualifier_value")[0].firstChild.data
        if list == "collection_date":
            key2=seqid
            dict_time[key2] = line.getElementsByTagName("INSDQualifier_value")[0].firstChild.data
print(dict_country)
print(dict_time)
with open("PEDV_info.txt","w")as F:
    country1=dict_country.keys()&dict_time.keys()
    for i in country1:
        print(i+"\t"+dict_country[i]+"\t"+dict_time[i])
        F.write(i+"\t"+dict_country[i]+"\t"+dict_time[i]+"\n")
    country2=dict_country.keys()-country1
    for ii in country2:
        print(ii+"\t"+dict_country[ii]+"\t"+"unkown")
        F.write(ii+"\t"+dict_country[ii]+"\t"+"unkown"+"\n")