#########################################################
#2022-12-04
#first write by Max
#test virus conservation analysis pipeline 
#########################################################
import os,argparse
import numpy as np

def Conservation(genomefile):
    genome_dict = {}
    with open("PEDV_rawdata.fasta","w")as f:
        with open(genomefile,"r")as F:
            for line in F:
                #print(line)
                if line.startswith(">"):
                    key = line.split()[0].rstrip()
                    genome_dict[key] = ""
                else:
                    genome_dict[key] += line.rstrip()
        for key, value in genome_dict.items():
            if len(value)>=20000:
                f.write(key+"\n")
                f.write(value+"\n")
    count=0; genome={};seqdif={}
    os.system("mafft --auto --thread 20 PEDV_rawdata.fasta  > PEDV_rawdata.aln.fasta")
    os.system("muscle -maketree -in PEDV_rawdata.aln.fasta -out PEDV_rawdata.aln.nwk -cluster neighborjoining")
    os.system("phyloFit PEDV_rawdata.aln.fasta --tree PEDV_rawdata.aln.nwk -out PEDV_rawdata.aln.mod")
    os.system("cat PEDV_rawdata.aln.nwk | tr -d '\n' > test.nwk")
    os.system("mkdir gerpdir")
    os.system("mkdir wigdir")
    with open("PEDV_rawdata.fasta","r")as F:
        for line in F:
            if line.startswith(">"):
                key = line.rstrip()
                genome[key] = ""
            else:
                genome[key] += line.rstrip()
    for key,value in genome.items():
        count+=1
        os.system("phastCons --target-coverage 0.9 --expected-length 12 --refidx %s PEDV_rawdata.aln.fasta ut.mod > %s.wig" % (count, key[1:]))
        os.system("mv %s.wig wigdir/" % (key[1:]))
        os.system("gerpcol -a -f  PEDV_rawdata.aln.fasta -t test.nwk  -e %s -x %s.rates"%(key[1:],key[1:]))
        os.system("mv *rates %s.rates"%(key[1:]))
        os.system("mv %s.rates gerpdir"%(key[1:]))
    os.system("python detection_G4.py -f PEDV_rawdata.fasta -g 2 -o PEDV_rawdata.g24.txt")
    os.system("python G4phastcons.py -wd wigdir/ -g PEDV_rawdata.g24.txt -o PEDV_rawdata.phastcons.txt")
    os.system("python conservation.py -AS PEDV_rawdata.aln.fasta -e 5 -grun 2 -o PEDV_rawdata.aln.patcon5.txt")
    os.system("Rscript g4hunterscore_real.r PEDV_rawdata.g24.txt PEDV_g4hunter.txt")
    os.system("python G4GerpScore.py -f PEDV_rawdata.fasta -gerp gerpdir -g PEDV_rawdata.g24.txt -o PEDV_rawdata.gerp.txt > gerp++_output.txt")

    with open("gerp++_output.txt","r")as F:
        for line in F:
            l=line.split(",")
            key=l[0].split()[0]
            seqdif[key]=int(l[0].split()[-1])-int(l[1].split()[1])
    for key,value in seqdif.items():
        file=os.path.join("gerpdir",key+".rates")
        df=np.loadtxt(file);df1=df[:value]
        df2=np.vstack((df,df1))
        np.savetxt(file,df2,fmt="%.3f",delimiter="\t")
    os.system("python G4GerpScore.py -f PEDV_rawdata.fasta -gerp gerpdir -g PEDV_rawdata.g24.txt -o PEDV_rawdata.gerp.txt > gerp++_output1.txt ")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To calculate the conservation. ")
    parser.add_argument("-f", "--genome", type=str,required=True,help="Please type in the genome file path.")
    Args = parser.parse_args()
    try:
        Conservation(os.path.abspath(Args.genome))
    except:
        pass




