"""\
 Ints.py Basic routines for integrals in the PyQuante framework

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from CGBF import CGBF,coulomb
#from contracted_gto import coulomb
from NumWrap import zeros,dot,reshape,array
from PyQuante.cints import ijkl2intindex as intindex
from PyQuante.cints import ijkl2intindexBTF as intindexBTF  # ADDED FOR COULOMB.py!
from PyQuante.Basis.Tools import get_basis_data
import settings
import logging

logger = logging.getLogger("pyquante")

sym2powerlist = {
    'S' : [(0,0,0)],
    'P' : [(1,0,0),(0,1,0),(0,0,1)],
    'D' : [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],  #GLOBULION MODIFICATION ! ---> to set this list same as in g09!!!
    'F' : [(3,0,0),(0,3,0),(0,0,3),(1,2,0),(2,1,0),(2,0,1),
           (1,0,2),(0,1,2),(0,2,1),(1,1,1)]
   #'F' : [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),(1,2,0),
   #       (0,2,1),(1,0,2),(0,1,2),(1,1,1)]
   #'D' : [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)], ---> old D-func sym2powerlist
   #'F' : [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),
   #       (0,3,0),(0,2,1),(0,1,2),(0,0,3)]
    }

sorted = True
jints = {}
kints = {}
jintsBTF_12 = {}  # GLOBULION ADD
kintsBTF_12 = {}  # GLOBULION ADD
jintsBTF_21 = {}  # GLOBULION ADD
kintsBTF_21 = {}  # GLOBULION ADD

def getbasis(atoms,basis_data=None,**opts):
    """\
    bfs = getbasis(atoms,basis_data=None)
    
    Given a Molecule object and a basis library, form a basis set
    constructed as a list of CGBF basis functions objects.
    """
    from PyQuante.Basis.basis import BasisSet
    return BasisSet(atoms, basis_data, **opts)
    # Option to omit f basis functions from imported basis sets
    omit_f = opts.get('omit_f',False)
    if not basis_data:
        from PyQuante.Basis.p631ss import basis_data
    elif type(basis_data) == type(''):
        # Assume this is a name of a basis set, e.g. '6-31g**'
        #  and import dynamically
        basis_data = get_basis_data(basis_data)
    bfs = []
    for atom in atoms:
        bs = basis_data[atom.atno]
        for sym,prims in bs:
            if omit_f and sym == "F": continue
            for power in sym2powerlist[sym]:
                bf = CGBF(atom.pos(),power,atom.atid)
                for expnt,coef in prims:
                    bf.add_primitive(expnt,coef)
                bf.normalize()
                bfs.append(bf)
    return bfs


def getints(bfs,atoms):
    logger.info("Calculating Integrals...")
    S,h = get1ints(bfs,atoms)
    Ints = get2ints(bfs)
    logger.info("Integrals Calculated.")
    return S,h,Ints

def get1ints(bfs,atoms):
    "Form the overlap S and h=t+vN one-electron Hamiltonian matrices"
    nbf = len(bfs)
    S = zeros((nbf,nbf),'d')
    h = zeros((nbf,nbf),'d')

    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(nbf):
            bfj = bfs[j]
            S[i,j] = bfi.overlap(bfj)
            h[i,j] = bfi.kinetic(bfj)
            for atom in atoms:
                h[i,j] = h[i,j] + atom.atno*bfi.nuclear(bfj,atom.pos())
    return S,h

def getT(bfs):
    "Form the kinetic energy matrix"
    nbf = len(bfs)
    T = zeros((nbf,nbf),'d')

    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(nbf):
            bfj = bfs[j]
            T[i,j] = bfi.kinetic(bfj)
    return T

def getS(bfs):
    "Form the overlap matrix"
    nbf = len(bfs)
    S = zeros((nbf,nbf),'d')

    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(nbf):
            bfj = bfs[j]
            S[i,j] = bfi.overlap(bfj)
    return S

def getSAB(bfsA,bfsB):
    """get the S_ij between two molecules A and B"""
    nbfA = len(bfsA)
    nbfB = len(bfsB)
    S = zeros((nbfA,nbfB),'d')
    
    for i in xrange(nbfA):
        bfi = bfsA[i]
        for j in xrange(nbfB):
            bfj = bfsB[j]
            S[i,j] = bfi.overlap(bfj)
    return S

def getTAB(bfsA,bfsB):
    """get the T_ij between two molecules A and B"""
    nbfA = len(bfsA)
    nbfB = len(bfsB)
    T = zeros((nbfA,nbfB),'d')
    
    for i in xrange(nbfA):
        bfi = bfsA[i]
        for j in xrange(nbfB):
            bfj = bfsB[j]
            T[i,j] = bfi.kinetic(bfj)
    return T


def getSA1B(bfsA,bfsB):
    """get the derivatives of S_ij between two molecules A and B.
The derivative refers to molecule A basis functions"""
    nbfA = len(bfsA)
    nbfB = len(bfsB)
    S = zeros((nbfA,nbfB,3),'d')
    
    for i in xrange(nbfA):
        bfi = bfsA[i]
        for j in xrange(nbfB):
            bfj = bfsB[j]
            S[i,j,:] =  bfi.overlap_1(bfj)
    return S

def getTA1B(bfsA,bfsB):
    """get the derivatives of T_ij between two molecules A and B.
The derivative refers to molecule A basis functions"""
    nbfA = len(bfsA)
    nbfB = len(bfsB)
    T = zeros((nbfA,nbfB,3),'d')
    
    for i in xrange(nbfA):
        bfi = bfsA[i]
        for j in xrange(nbfB):
            bfj = bfsB[j]
            T[i,j,:] = array( bfi.kinetic_1(bfj) , 'd' )
    return T

#  GLOBULION ADD
def getM(bfs, hexadecapoles=False):
    "Form the multipole intragral supermatrixes"
    nbf = len(bfs)
    D = zeros((3,      nbf,nbf),'d')
    Q = zeros((3,3,    nbf,nbf),'d')
    O = zeros((3,3,3,  nbf,nbf),'d')
    H = zeros((3,3,3,3,nbf,nbf),'d')

    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(i+1):
            bfj = bfs[j]
            # dipole integral matrix
            D[0,i,j] = bfi.multipole(bfj,1,0,0)
            D[1,i,j] = bfi.multipole(bfj,0,1,0)
            D[2,i,j] = bfi.multipole(bfj,0,0,1)
            D[:,j,i] = D[:,i,j]
            # quadrupole integral matrix
            Q[0,0,i,j] = bfi.multipole(bfj,2,0,0)
            Q[1,1,i,j] = bfi.multipole(bfj,0,2,0)
            Q[2,2,i,j] = bfi.multipole(bfj,0,0,2)
                 
            a =  bfi.multipole(bfj,1,1,0)           
            Q[0,1,i,j] = a
            Q[1,0,i,j] = a
            a = bfi.multipole(bfj,1,0,1)
            Q[0,2,i,j] = a
            Q[2,0,i,j] = a
            a = bfi.multipole(bfj,0,1,1)
            Q[1,2,i,j] = a
            Q[2,1,i,j] = a
            Q[:,:,j,i] = Q[:,:,i,j]
            # octupole integral matrix
            O[0,0,0,i,j] = bfi.multipole(bfj,3,0,0)
            O[1,1,1,i,j] = bfi.multipole(bfj,0,3,0)
            O[2,2,2,i,j] = bfi.multipole(bfj,0,0,3)

            a = bfi.multipole(bfj,2,1,0)
            O[0,0,1,i,j] = a
            O[0,1,0,i,j] = a
            O[1,0,0,i,j] = a
 
            a = bfi.multipole(bfj,2,0,1)
            O[0,0,2,i,j] = a
            O[0,2,0,i,j] = a
            O[2,0,0,i,j] = a

            a = bfi.multipole(bfj,1,2,0)
            O[1,1,0,i,j] = a
            O[1,0,1,i,j] = a
            O[0,1,1,i,j] = a

            a = bfi.multipole(bfj,0,2,1)
            O[1,1,2,i,j] = a
            O[1,2,1,i,j] = a
            O[2,1,1,i,j] = a

            a = bfi.multipole(bfj,0,1,2)
            O[1,2,2,i,j] = a
            O[2,1,2,i,j] = a
            O[2,2,1,i,j] = a

            a = bfi.multipole(bfj,1,0,2)
            O[0,2,2,i,j] = a
            O[2,0,2,i,j] = a
            O[2,2,0,i,j] = a

            a = bfi.multipole(bfj,1,1,1)
            O[0,1,2,i,j] = a
            O[0,2,1,i,j] = a
            O[1,2,0,i,j] = a
            O[1,0,2,i,j] = a
            O[2,0,1,i,j] = a
            O[2,1,0,i,j] = a
            O[:,:,:,j,i] = O[:,:,:,i,j]

            if hexadecapoles:
               a = bfi.multipole(bfj,1,0,3)
               # XZZZ
               H[0,2,2,2,i,j] = a
               # ZXZZ
               H[2,0,2,2,i,j] = a
               # ZZXZ
               H[2,2,0,2,i,j] = a
               # ZZZX
               H[2,2,2,0,i,j] = a
               
               a = bfi.multipole(bfj,2,0,2)
               # XXZZ
               H[0,0,2,2,i,j] = a
               # XZXZ
               H[0,2,0,2,i,j] = a
               # XZZX
               H[0,2,2,0,i,j] = a
               # ZXXZ
               H[2,0,0,2,i,j] = a
               # ZXZX
               H[2,0,2,0,i,j] = a
               # ZZXX
               H[2,2,0,0,i,j] = a
               
               a = bfi.multipole(bfj,0,4,0)
               # YYYY
               H[1,1,1,1,i,j] = a
               
               a = bfi.multipole(bfj,1,2,1)
               # XYYZ
               H[0,1,1,2,i,j] = a
               # XYZY
               H[0,1,2,1,i,j] = a
               # XZYY
               H[0,2,1,1,i,j] = a
               # YXYZ
               H[1,0,1,2,i,j] = a
               # YXZY
               H[1,0,2,1,i,j] = a
               # YYXZ
               H[1,1,0,2,i,j] = a
               # YYZX
               H[1,1,2,0,i,j] = a
               # YZXY
               H[1,2,0,1,i,j] = a
               # YZYX
               H[1,2,1,0,i,j] = a
               # ZXYY
               H[2,0,1,1,i,j] = a
               # ZYXY
               H[2,1,0,1,i,j] = a
               # ZYYX
               H[2,1,1,0,i,j] = a
               
               a = bfi.multipole(bfj,3,1,0)
               # XXXY
               H[0,0,0,1,i,j] = a
               # XXYX
               H[0,0,1,0,i,j] = a
               # XYXX
               H[0,1,0,0,i,j] = a
               # YXXX
               H[1,0,0,0,i,j] = a
               
               a = bfi.multipole(bfj,0,1,3)
               # YZZZ
               H[1,2,2,2,i,j] = a
               # ZYZZ
               H[2,1,2,2,i,j] = a
               # ZZYZ
               H[2,2,1,2,i,j] = a
               # ZZZY
               H[2,2,2,1,i,j] = a
               
               a = bfi.multipole(bfj,1,3,0)
               # XYYY
               H[0,1,1,1,i,j] = a
               # YXYY
               H[1,0,1,1,i,j] = a
               # YYXY
               H[1,1,0,1,i,j] = a
               # YYYX
               H[1,1,1,0,i,j] = a
               
               a = bfi.multipole(bfj,0,3,1)
               # YYYZ
               H[1,1,1,2,i,j] = a
               # YYZY
               H[1,1,2,1,i,j] = a
               # YZYY
               H[1,2,1,1,i,j] = a
               # ZYYY
               H[2,1,1,1,i,j] = a
               
               a = bfi.multipole(bfj,1,1,2)
               # XYZZ
               H[0,1,2,2,i,j] = a
               # XZYZ
               H[0,2,1,2,i,j] = a
               # XZZY
               H[0,2,2,1,i,j] = a
               # YXZZ
               H[1,0,2,2,i,j] = a
               # YZXZ
               H[1,2,0,2,i,j] = a
               # YZZX
               H[1,2,2,0,i,j] = a
               # ZXYZ
               H[2,0,1,2,i,j] = a
               # ZXZY
               H[2,0,2,1,i,j] = a
               # ZYXZ
               H[2,1,0,2,i,j] = a
               # ZYZX
               H[2,1,2,0,i,j] = a
               # ZZXY
               H[2,2,0,1,i,j] = a
               # ZZYX
               H[2,2,1,0,i,j] = a
               
               a = bfi.multipole(bfj,0,2,2)
               # YYZZ
               H[1,1,2,2,i,j] = a
               # YZYZ
               H[1,2,1,2,i,j] = a
               # YZZY
               H[1,2,2,1,i,j] = a
               # ZYYZ
               H[2,1,1,2,i,j] = a
               # ZYZY
               H[2,1,2,1,i,j] = a
               # ZZYY
               H[2,2,1,1,i,j] = a
               
               a = bfi.multipole(bfj,4,0,0)
               # XXXX
               H[0,0,0,0,i,j] = a
               
               a = bfi.multipole(bfj,2,1,1)
               # XXYZ
               H[0,0,1,2,i,j] = a
               # XXZY
               H[0,0,2,1,i,j] = a
               # XYXZ
               H[0,1,0,2,i,j] = a
               # XYZX
               H[0,1,2,0,i,j] = a
               # XZXY
               H[0,2,0,1,i,j] = a
               # XZYX
               H[0,2,1,0,i,j] = a
               # YXXZ
               H[1,0,0,2,i,j] = a
               # YXZX
               H[1,0,2,0,i,j] = a
               # YZXX
               H[1,2,0,0,i,j] = a
               # ZXXY
               H[2,0,0,1,i,j] = a
               # ZXYX
               H[2,0,1,0,i,j] = a
               # ZYXX
               H[2,1,0,0,i,j] = a
               
               a = bfi.multipole(bfj,2,2,0)
               # XXYY
               H[0,0,1,1,i,j] = a
               # XYXY
               H[0,1,0,1,i,j] = a
               # XYYX
               H[0,1,1,0,i,j] = a
               # YXXY
               H[1,0,0,1,i,j] = a
               # YXYX
               H[1,0,1,0,i,j] = a
               # YYXX
               H[1,1,0,0,i,j] = a
               
               a = bfi.multipole(bfj,0,0,4)
               # ZZZZ
               H[2,2,2,2,i,j] = a
               
               a = bfi.multipole(bfj,3,0,1)
               # XXXZ
               H[0,0,0,2,i,j] = a
               # XXZX
               H[0,0,2,0,i,j] = a
               # XZXX
               H[0,2,0,0,i,j] = a
               # ZXXX
               H[2,0,0,0,i,j] = a

               H[:,:,:,:,j,i] = H[:,:,:,:,i,j]


    return D,Q,O,H

def getV(bfs,atoms):
    "Form the nuclear attraction matrix V"
    nbf = len(bfs)
    V = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(nbf):
            bfj = bfs[j]
            for atom in atoms:
                V[i,j] = V[i,j] + atom.atno*bfi.nuclear(bfj,atom.pos())
    return V

def getVEFP(bfs1,bfs2,q,r):
    """potential matrix elements. q - charges, r - their positions"""
    nbf1 = len(bfs1)
    nbf2 = len(bfs2)
    V = zeros((nbf1,nbf2),'d')
    for i in xrange(nbf1):
        bfi = bfs1[i]
        for j in xrange(nbf2):
            bfj = bfs2[j]
            v = 0
            for a in xrange(len(q)):
                v += q[a] * bfi.nuclear(bfj,tuple(r[a]))
            V[i,j] = v
    return V

if settings.libint_enabled == True:
    # Libint Integrals
    import numpy as np
    import clibint
    
    def get2ints(basis):
        lenbasis = len(basis.bfs)
        
        Ints = np.zeros((lenbasis**4),dtype=np.float64)

        for i,a in enumerate(basis.shells):
            for j,b in enumerate(basis.shells[:i+1]):
                for k,c in enumerate(basis.shells):
                    for l,d in enumerate(basis.shells[:k+1]):
                        if (i+j)>=(k+l):
                            clibint.shell_compute_eri(a,b,c,d,Ints)
        if sorted:
            sortints(lenbasis,Ints)
        return Ints

    def get2intsBTF(basis1,basis2):
        """this modified version of get2ints computes 2-el ints 
           only for interacting basis functions"""
        lenbasis1 = len(basis1.bfs)
        lenbasis2 = len(basis2.bfs)
        
        #Ints = np.zeros((lenbasis1**2*lenbasis2**2),dtype=np.float64)
        Ints = np.zeros((100000000),dtype=np.float64)

        # compute only (ab|cd) integrals where a,b belong to basis1 
        #                                  and c,d belong to basis2
        for i,a in enumerate(basis1.shells):
            for j,b in enumerate(basis1.shells[:i+1]):
                for k,c in enumerate(basis2.shells):
                    for l,d in enumerate(basis2.shells[:k+1]):
                        if (i>=j) and (k>=l):
                            clibint.shell_compute_eri(a,b,c,d,Ints)
        if sorted:
           sortintsBTF_12(lenbasis1,lenbasis2,Ints)
           sortintsBTF_21(lenbasis1,lenbasis2,Ints)
        return Ints

else:
    # PyQuante Integrals
    def get2ints(bfs):
        """Store integrals in a long array in the form (ij|kl) (chemists
        notation. We only need i>=j, k>=l, and ij <= kl"""

        from array import array
        nbf = len(bfs)
        totlen = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
        Ints = array('d',[0]*totlen)

        for i in xrange(nbf):
            for j in xrange(i+1):
                ij = i*(i+1)/2+j
                for k in xrange(nbf):
                    for l in xrange(k+1):
                        kl = k*(k+1)/2+l
                        if ij <= kl:
                            Ints[intindex(i,j,k,l)] = coulomb(bfs[i],bfs[j],
                                                              bfs[k],bfs[l])

        if sorted:
            sortints(nbf,Ints)
        return Ints

def sortints(nbf,Ints):
    for i in xrange(nbf):
        for j in xrange(i+1):
            jints[i,j] = fetch_jints(Ints,i,j,nbf)
            kints[i,j] = fetch_kints(Ints,i,j,nbf)
    return

# THIS FUNCTION IS ADDED FOR COULOMB.py PURPOSES
def sortintsBTF_12(nbf1,nbf2,Ints):
    for i in xrange(nbf1):
        for j in xrange(i+1):
            jintsBTF_12[i,j] = fetch_jintsBTF_12(Ints,i,j,nbf1,nbf2)
            kintsBTF_12[i,j] = fetch_kintsBTF_12(Ints,i,j,nbf1,nbf2)
    return

# THIS FUNCTION IS ADDED FOR COULOMB.py PURPOSES
def sortintsBTF_21(nbf1,nbf2,Ints):
    for i in xrange(nbf2):
        for j in xrange(i+1):
            jintsBTF_21[i,j] = fetch_jintsBTF_21(Ints,i,j,nbf1,nbf2)
            kintsBTF_21[i,j] = fetch_kintsBTF_21(Ints,i,j,nbf1,nbf2)
    return

# THIS FUNCTION IS ADDED FOR COULOMB.py PURPOSES
def fetch_jintsBTF_12(Ints,i,j,nbf1,nbf2):
    temp = zeros(nbf2*nbf2,'d')
    kl = 0
    for k in xrange(nbf2):
        for l in xrange(nbf2):
            index = intindexBTF(i,j,k,l)
            temp[kl] = Ints[index]
            kl+=1
    return temp

# THIS FUNCTION IS ADDED FOR COULOMB.py PURPOSES
def fetch_kintsBTF_12(Ints,i,j,nbf1,nbf2):
    temp = zeros(nbf2*nbf2,'d')
    kl = 0
    for k in xrange(nbf2):
        for l in xrange(nbf2):
            index = intindexBTF(i,k,j,l)
            temp[kl] = Ints[index]
            kl+=1
    return temp

# THIS FUNCTION IS ADDED FOR COULOMB.py PURPOSES
def fetch_jintsBTF_21(Ints,i,j,nbf1,nbf2):
    temp = zeros(nbf1*nbf1,'d')
    kl = 0
    for k in xrange(nbf1):
        for l in xrange(nbf1):
            index = intindexBTF(i,j,k,l)
            temp[kl] = Ints[index]
            kl+=1
    return temp

# THIS FUNCTION IS ADDED FOR COULOMB.py PURPOSES
def fetch_kintsBTF_21(Ints,i,j,nbf1,nbf2):
    temp = zeros(nbf1*nbf1,'d')
    kl = 0
    for k in xrange(nbf1):
        for l in xrange(nbf1):
            index = intindexBTF(i,k,j,l)
            temp[kl] = Ints[index]
            kl+=1
    return temp

def fetch_jints(Ints,i,j,nbf):
    temp = zeros(nbf*nbf,'d')
    kl = 0
    for k in xrange(nbf):
        for l in xrange(nbf):
            index = intindex(i,j,k,l)
            temp[kl] = Ints[index]
            kl += 1
    return temp

def fetch_kints(Ints,i,j,nbf):
    temp = zeros(nbf*nbf,'d')
    kl = 0
    for k in xrange(nbf):
        for l in xrange(nbf):
            temp[kl] = Ints[intindex(i,k,j,l)]
            kl += 1
    return temp

# THIS FUNCTION IS ADDED FOR COULOMB.py PURPOSES
def fetch_jintsBTF(Ints,i,j,nbf1,nbf2):
    temp = zeros(nbf1*nbf1,'d')
    kl = 0
    for k in xrange(nbf2):
        for l in xrange(nbf2):
            index = intindexBTF(i,j,k,l)
            temp[kl] = Ints[index]
            kl+=1
    return temp

# THIS FUNCTION IS ADDED FOR COULOMB.py PURPOSES
def fetch_kintsBTF(Ints,i,j,nbf1,nbf2):
    temp = zeros(nbf1*nbf1,'d')
    kl = 0
    for k in xrange(nbf1):
        for l in xrange(nbf1):
            index = intindexBTF(i,k,j,l)
            temp[kl] = Ints[index]
            kl+=1
    return temp

def getJ(Ints,D):
    "Form the Coulomb operator corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    J = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        for j in xrange(i+1):
            if sorted:
                temp = jints[i,j]
            else:
                temp = fetch_jints(Ints,i,j,nbf)
            J[i,j] = dot(temp,D1d)
            J[j,i] = J[i,j]
    return J

# THIS FUNCTION HAS BEEN ADDED FOR COULOMB.py PURPOSES
def getJ_BTF(Ints,D,nbf1,nbf2,frag_number):
    if frag_number == 1:
       D1d = reshape(D,(nbf2*nbf2,)) #1D version of Dens
       J = zeros((nbf1,nbf1),'d')
       for i in xrange(nbf1):
           for j in xrange(i+1):
               if sorted:
                  temp = jintsBTF_12[i,j]
               else:
                  temp = fetch_jintsBTF_12(Ints,i,j,nbf1,nbf2)
               J[i,j] = dot(temp,D1d)
               J[j,i] = J[i,j]
    elif frag_number == 2:
       D1d = reshape(D,(nbf1*nbf1,)) #1D version of Dens
       J = zeros((nbf2,nbf2),'d')
       for i in xrange(nbf2):
           for j in xrange(i+1):
               if sorted:
                  temp = jintsBTF_21[i,j]
               else:
                  temp = fetch_jintsBTF_21(Ints,i,j,nbf1,nbf2)
               J[i,j] = dot(temp,D1d)
               J[j,i] = J[i,j]
    return J

# THIS FUNCTION HAS BEEN ADDED FOR COULOMB.py PURPOSES
def getK_BTF(Ints,D,nbf1,nbf2,frag_number):
    if frag_number == 1:
       D1d = reshape(D,(nbf2*nbf2,)) #1D version of Dens
       K = zeros((nbf1,nbf1),'d')
       for i in xrange(nbf1):
           for j in xrange(i+1):
               if sorted:
                  temp = kintsBTF_12[i,j]
               else:
                  temp = fetch_kintsBTF_12(Ints,i,j,nbf1,nbf2)
               K[i,j] = dot(temp,D1d)
               K[j,i] = K[i,j]
    elif frag_number == 2:
       D1d = reshape(D,(nbf1*nbf1,)) #1D version of Dens
       K = zeros((nbf2,nbf2),'d')
       for i in xrange(nbf2):
           for j in xrange(i+1):
               if sorted:
                  temp = kintsBTF_21[i,j]
               else:
                  temp = fetch_kintsBTF_21(Ints,i,j,nbf1,nbf2)
               K[i,j] = dot(temp,D1d)
               K[j,i] = K[i,j]
    return K

def getK(Ints,D):
    "Form the exchange operator corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    K = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        for j in xrange(i+1):
            if sorted:
                temp = kints[i,j]
            else:
                temp = fetch_kints(Ints,i,j,nbf)
            K[i,j] = dot(temp,D1d)
            K[j,i] = K[i,j]
    return K

def get2JmK(Ints,D):
    "Form the 2J-K integrals corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    G = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        for j in xrange(i+1):
            if sorted:
                temp = 2*jints[i,j]-kints[i,j]
            else:
                temp = 2*fetch_jints(Ints,i,j,nbf)-fetch_kints(Ints,i,j,nbf)
            G[i,j] = dot(temp,D1d)
            G[j,i] = G[i,j]
    return G

