import sys, os 
_libraries = ["\driver"]
for _l in _libraries:
    s = "e:\PBC\GitHub\pyrmcpot\\rmcpot\\tools\\rmcpot_pygmo\\" + _l
    sys.path.append(os.path.abspath(s))

import rmcpot_driver as driver

nind=9
ngen=500
mintype="global"
minalgo="pso"

driver.run(mintype,minalgo,ngen,nind) 