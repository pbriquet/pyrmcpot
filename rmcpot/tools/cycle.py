import os
import subprocess as sp
import shlex
from ButterflyChaser import ButterflyChaser

hops=200
hopstart=1
bestparamfile='bestparamset.dat'
cmd='rmcpot &'

bc=ButterflyChaser()
bc.rundirs=['run1','run2','run3','run4','run5','run6']
bc.paramfiles=[bestparamfile for i in range(len(bc.rundirs))]
bc.changeration=0.0001

for i in range(hopstart,hops):
    print("Step %d running" % i)

    bc.readparams()
    bc.updateparams()
    bc.writeparams('paramset.dat.0')
    
    j=0
    
    for rundir in bc.rundirs:      
        os.chdir(rundir)
        
        with open(bestparamfile,'r') as f:
           line=f.readline().split()
           val=float(line[0])

        with open('errorhistory.dat','a') as f:
            f.write("%d %f \n" % (i,val))

        os.chdir('..')

        j+=1

    j=0
    args=shlex.split(cmd)

    for rundir in bc.rundirs:
        os.chdir(rundir)

        if j!=bc.__indexmin__:
            os.system('mv '+bestparamfile+' '+bestparamfile+'.bkp')
            
            process=sp.Popen(args,stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE)
            run=process.communicate()

            if not os.path.exists(bestparamfile):
                os.system('mv '+bestparamfile+'.bkp '+bestparamfile)    

        if os.path.exists('paramset.dat.0'):
            os.system('rm paramset.dat.0')
            
        os.chdir('..')

        j+=1

print("Done!")
