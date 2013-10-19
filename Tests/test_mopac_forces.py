from time import time
from PyQuante.Molecule import h2o,h2,oh
from PyQuante.MINDO3 import scf,get_energy_forces,numeric_forces


h2o_energy = scf(h2o,verbose=True)
oh_energy = scf(oh,verbose=True)

print numeric_forces(h2o)
print numeric_forces(oh)
