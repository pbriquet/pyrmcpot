import subprocess as sp
import shlex
import io
import math as m

class rmcpot_prob:
    def fitness(self,x):
        try:
            sp.Popen(["rm bestparamset.dat"]).wait()
        except:
            pass
        
        with open("paramset.dat.0","w") as f:
            f.write("1000000000000.0 %d\n" % len(x))
            
            for i in range(len(x)):
                f.write("%f\n" % x[i])
        
        finished=False
        cmd='rmcpot &'
        args=shlex.split(cmd)
        run=sp.Popen(args,stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE).communicate()
        out=io.StringIO(run[0].decode())
        lines=out.readlines()
        
        for line in lines:
            l=line.rstrip().split()
            
            if l[0]==">>>>>>>":
                finished=True
                val=float(l[1])
                
                break
        
        if finished:
            return [val]
        else:
            return [m.inf]                

    def get_bounds(self):
        func=[[],[]]
        indlow=1
        indhigh=2
        
        with open("paramlist2",'r') as f:
            lines=f.readlines()
            
        for line in lines:            
            l=line.split()
            lbound=float(l[indlow])
            hbound=float(l[indhigh])
            
            func[0].append(lbound)
            func[1].append(hbound)
            
        return tuple(func)
            