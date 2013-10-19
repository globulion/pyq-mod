#!/usr/bin/env python
from pylab import *
from PyQuante import SCF,Molecule


def energy(R=1.217,method='UHF'):
    O2 = Molecule('O2',atomlist=[('O',(0,0,0)),('O',(R,0,0))],\
                  units='Angstrom',multiplicity=3)
    O2abinitio = SCF(O2,method=method,basis='6-31G**')
    O2abinitio.iterate()
    return O2abinitio.energy

def Escan(method='UHF'):
    data = []
    for R in [0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2]:
        E = energy(R,method)
        print R,E
        data.append((R,E))
    data = array(data)
    plot(data[:,0],data[:,1],'bo-')
    show()
    savefig('o2-rohf.png')
    return

def Edot(method='ROHF'):
    O = Molecule('O',atomlist=[('O',(0,0,0))],multiplicity=3)
    job = SCF(O,method=method,basis='6-31G**')
    job.iterate()
    return job.energy

if __name__ == '__main__':
    method='UHF'
    EO = Edot(method)
    print "O energy: ",EO,2*EO
    Escan(method)
