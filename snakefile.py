#########################################################
#2022-12-04
#first write by Max
#test virus conservation analysis pipeline with snakemake
#########################################################

rep_index={"PEDV"}


rule all:
    input:
        ""

rule clean_data:
    input:
        fa="rawdata/{sample}.fasta"
    output:
        fa="rawdata/{sample}.rawdata.fasta"
    threads: 20
    run:
        genome_dict = {}
        with open({output},w)as f:
            with open({input},r)as F:
                for line in F:
                    if line.startswith(">"):
                        key = line.split()[0].rstrip()
                        genome_dict[key] = ""
                    else:
                        genome_dict[key] += line.rstrip()
            for key, value in genome_dict.items():
                if len(value)>=20000:
                    f.write(key+"\n")
                    f.write(value+"\n")

rule sequence_conservation:
    input:
        fa="rawdata/{sample}.rawdata.fasta"
    output:
        fa="rawdata/{sample}.rawdata.aln.fasta"
        nwk="rawdata/{sample}.rawdata.aln.nwk"
        mod="seqcons/{sample}.rawdata.aln.mod"
        txt="rawdata/{sample}.rawdata.g24.txt"
        txt="seqcons/{sample}.rawdata.phastcons.txt"
    threads: 20
    run:
        import os
        genome_dict = {};count=0
        os.system("mafft --auto {input}  > {output[0]}")
        os.system("muscle -maketree -in {output[0]} -out {output[1]} -cluster neighborjoining")
        os.system("phyloFit {output[0]} --tree output[1] -out output[2]")
        os.system("mkdir wigdir")
        with open({input},"r")as F:
            for line in F:
                if line.startswith(">"):
                    key = line.rstrip()
                    genome_dict[key] = ""
                else:
                    genome_dict[key] += line.rstrip()
        for key,value in genome_dict.items():
            count+=1
            os.system("phastCons --target-coverage 0.9 --expected-length 12 --refidx %s {output[0]} ut.mod > %s.wig" % (count, key[1:]))
            os.system("mv %s.wig wigdir/" % (key[1:]))
        os.system("python detection_G4.py -f {input} -g 2 -o {output[3]}")
        os.system("python G4phastcons.py -wd wigdir/ -g {output[3]} -o {output[4]}")

rule pattern_conservation:
    input:
        fa="rawdata/{sample}.rawdata.aln.fasta"
    output:
        txt="patcons/{sample}.rawdata.aln.patcons.txt"
    threads: 20
    run:
        import os
        os.system("python conservation.py -AS {input} -e 5 -grun 2 -o {output}")

rule g4hunter:
    input:
        txt="rawdata/{sample}.rawdata.g24.txt"
    output:
        txt="g4hunter/{sample}.rawdata.aln.g4hunter.txt"
    threads: 20
    run:
        Rscript g4hunterscore_real.r {input} {output}

rule gerp:
    input:
        nwk="rawdata/{sample}.rawdata.aln.nwk"
        fa="rawdata/{sample}.rawdata.fasta"
        fa="rawdata/{sample}.rawdata.aln.fasta"
        txt="rawdata/{sample}.rawdata.g24.txt"
    output:
        nwk="rawdata/{sample}.rawdata.aln.test.nwk"
        txt="gerp/{sample}.rawdata.aln.gerp.txt"
    threads: 20
    log:
        log="log/{sample}.rawdata.aln.first.log"
        log="log/{sample}.rawdata.aln.second.log"
    run:
        import os
        import numpy as np
        os.system("cat {input[0]} | tr -d '\n' > {output[0]}")
        os.system("mkdir gerpdir")
        genome_dict = {}
        with open('{input[1]}', 'r') as F:
            for line in F:
                if line.startswith(">"):
                    key = line.rstrip()
                    genome_dict[key] = ""
                else:
                    genome_dict[key] += line.rstrip()
        a=0
        for key,value in genome_dict.items():
            a+=1
            os.system("gerpcol -a -f  {input[2]} -t test.nwk  -e %s -x %s.rates"%(key[1:],key[1:]))
            os.system("mv *rates %s.rates"%(key[1:]))
            os.system("mv %s.rates gerpdir"%(key[1:]))
        os.system("python G4GerpScore.py -f {input[1]} -gerp gerpdir -g {input[3]} -o {output[1]} > {log[0]} 2>&1")
        seqdif={}
        with open({log[0]},"r")as F:
            for line in F:
                l=line.split(",")
                key=l[0].split()[0]
                seqdif[key]=int(l[0].split()[-1])-int(l[1].split()[1])
        for key,value in seqdif.items():
            file=os.path.join("gerpdir",key+".rates")
            df=np.loadtxt(file);df1=df[:value]
            df2=np.vstack((df,df1))
            np.savetxt(file,df2,fmt="%.3f",delimiter="\t")
        os.system("python G4GerpScore.py -f {input[1]} -gerp gerpdir -g {input[3]} -o {output[1]} >  {log[1]} 2>&1 ")




