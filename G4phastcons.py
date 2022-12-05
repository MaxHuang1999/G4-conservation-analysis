import re,argparse,os


def SeqConservation(wigdir,seqlocation_file,opfile):
    def returnloc(start,end, seq, nrun,gtype):
        gruns = re.finditer('%s{%s}.' % (gtype,nrun), seq)
        gruns_location = [(start + g.start(), start + g.end() - 2) for g in gruns]
        gruns_location[-1]=(end-nrun+1,end)
        loops_location = [(gruns_location[i][1]+1, gruns_location[i + 1][0]-1) for i in range(len(gruns_location)-1)]
        return (gruns_location, loops_location)
    wigdict={}
    for file in os.listdir(wigdir):
        key=file.split('.')[0]
        with open("/".join([wigdir,file]),'r') as F:
            wigscore=[float(i.strip()) for i in F if re.match('\d',i)]
            wigdict[key]=wigscore
    oplines=[]
    with open(seqlocation_file,'r') as F:
        for line in F:
            splitlines=line.rstrip().split(',')
            shortkey=splitlines[0].split('.')[0]
            if shortkey not in wigdict:
                continue
            start,end=int(splitlines[2]),int(splitlines[3])
            meanscore=round(sum(wigdict[shortkey][start-1:end])/(end-start+1),3)
            if splitlines[1] not in ['GGG','CCC','GG','CC']:
                splitlines.extend([str(meanscore),'None','None'])
                oplines.append(','.join(splitlines)+'\n')
                continue
            seq=splitlines[4]
            gtype=seq[0]
            grun=len(splitlines[1])
            gruns_loc,loops_loc=returnloc(start,end,seq,grun,gtype)
            gruns_score=sum([sum(wigdict[shortkey][loc[0]-1:loc[1]]) for loc in gruns_loc])
            gruns_length=sum([loc[1]-loc[0]+1 for loc in gruns_loc])
            mean_gruns=round(gruns_score/gruns_length,3)
            loops_score=sum([sum(wigdict[shortkey][loc[0]-1:loc[1]]) for loc in loops_loc])
            loops_length=sum([loc[1]-loc[0]+1 for loc in loops_loc])
            mean_loops=round(loops_score/loops_length,3)
            splitlines.extend([str(meanscore),str(mean_gruns),str(mean_loops)])
            oplines.append(','.join(splitlines)+'\n')
    with open(opfile,'w') as F:
        F.writelines(oplines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To calculate the phastcons score. First-writtern@20190731 by ZhenLi.")
    parser.add_argument("-wd", "--wigdir", type=str,required=True,help="Please type in the directory of your wigfiles.")
    parser.add_argument("-g","--g4file",type=str,required=True,help="Please type in your g4file.")
    parser.add_argument("-o", "--opfile", required=True, type=str, help="The output file name.")
    Args = parser.parse_args()
    try:
        SeqConservation(os.path.abspath(Args.wigdir),os.path.abspath(Args.g4file),os.path.abspath(Args.opfile))
    except:
        pass
