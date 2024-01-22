import pygmo as pg
import numpy as np

import sys, os 
_libraries = ["\\problem"]
for _l in _libraries:
    s = "e:\PBC\GitHub\pyrmcpot\\rmcpot\\tools\\rmcpot_pygmo\\" + _l
    sys.path.append(os.path.abspath(s))


from rmcpot_prob import rmcpot_prob

def addprevchampion(population,index=None):
    try:
        g=open("besttillnow","r")
        lines=g.readlines()
    
        g.close()
    except:
        return population
    
    x=[]
    
    for line in lines:
        l=line.split()
        
        if len(l)==2:
            f=float(l[0])
        else:
            x.append(float(l[0]))
            
    f=np.array(f)
    x=np.array(x)
    
    if index is None:
        population.push_back(x,f)
    else:
        population.set_xf(index,x,f)
    
    return population

def savechampion(population):
    g=open("besttillnow","w")
    
    g.write("%f %d\n" % (population.champion_f[0],len(population.champion_x)))
    
    for i in range(len(population.champion_x)):
        g.write("%f\n" % population.champion_x[i])
        
    g.close()

def run(mintype="global",minalgo="pso",ngen=1,nind=1):
    tmp = rmcpot_prob()
    prob=pg.problem(tmp)
    
    if mintype=="global":
        print("Performing global minimization:\n-------------------------------\n")
        
        pop=pg.population(prob,nind)
        pop=addprevchampion(pop)
        
        if minalgo=="pso":
            print("\nUsing Particle Swarm Optimization algorithm...\n")
            
            algo=pg.algorithm(pg.pso(gen=ngen))
        else:
            print("Not an accepted global minimization method!")
    elif mintype=="local":
        print("Performing local minimization:\n------------------------------\n")
        
        pop=pg.population(prob,1)
        pop=addprevchampion(pop,0)
        
        if minalgo=="powell":
            print("\nUsing Powell algorithm...\n")
            
            algo=pg.algorithm(pg.scipy_optimize(method="Powell"))
        else:
            print("Not an accepted local minimization method!")
    else:
        print("Minimization type must be either 'global' or 'local'!")
        
        return
            
    pop=algo.evolve(pop)
    
    print("\nBest fitness: %f\n" % pop.champion_f[0])
    print("\nParameters:\n")
    
    for i in range(len(pop.champion_x)):
        print("   %f\n" % pop.champion_x[i])
        
    savechampion(pop)
    
    print("\nDone!")