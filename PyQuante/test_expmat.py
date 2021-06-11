
from rohf import expmat
from scipy.linalg.matfuncs import expm
from numpy import zeros,identity,array

n = 2
I = identity(n,'d')
#A = zeros((n,n),'d')
#A = I/3.
#A = array([[-49,24],[-64,31]],'d') # Example from Moler/Van Loan, doesn't work
#A = array([[0,0.5],[0.5,0]],'d')
A = array([[1,0.5],[0.5,1]],'d')
print(A)

E = expmat(A,nmax=20)
print(E)
print(expm(A))
