import numpy as np
def get_resscore(filename,window):
    p={}
    with open(filename) as result:
        for l in result:
            if(l.startswith('ATOM')):
                resn=l[22:26].replace(' ','')
                #print('|'+l[30:38]+'|'+l[38:46]+'|'+l[46:54]+'|',resn)
                res_name = l[17:20]
                x=float(l[30:38])
                y=float(l[38:46])
                z=float(l[46:54])
                cd=([x,y,z])
                sco=float(l[61:67])
                p[resn]=[cd,sco,res_name]
                #print('Res',resn,p[resn])
        #Window Score
        for resn in p:
            cnt=1
            sco=p[resn][1]
            c=str(int(resn)-window)
            if c in p:
                sco=sco+p[c][1]
                cnt=cnt+1
            c=str(int(resn)+window)
            if c in p:
                sco=sco+p[c][1]
                cnt=cnt+1
            p[resn].append(sco/float(cnt))
    return p
def get_pdb(filename):
    p={}
    with open(filename) as result:
        for l in result:
            if(l.startswith('ATOM') and l[13:16] == 'CA '):
                resn=l[22:26].replace(' ','')
                #print('|'+l[30:38]+'|'+l[38:46]+'|'+l[46:54]+'|',resn)
                x=float(l[30:38])
                y=float(l[38:46])
                z=float(l[46:55])
                cd=[x,y,z]
                sco=float(l[61:67])
                AA = l[17:20]
                p[resn]=[cd,AA]
    return p

def save_pdb_with_score(p,ch,filename):
    
    # p: pdb checking dict [coord,amino acid]
    # d: score matching dict
    output = open(filename, 'w')
    Natm=1
    for resn in p:
        #sco = d[resn][window+1]
        sco = -p[resn][3] #Opposit!! Lower is better for pymol
        line='ATOM{:7d}  CA  {:3} {:1}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}\n'.format(Natm,p[resn][2] ,ch,int(resn),p[resn][0][0],p[resn][0][1],p[resn][0][2],sco)
        Natm = Natm+1
        output.write(line)
    output.close()