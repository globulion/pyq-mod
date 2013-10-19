import unittest

from PyQuante.TestMolecules import h,li,oh
from PyQuante import SCF

class ODFTTests(unittest.TestCase):
    def testH(self):
        solv = SCF(h,method='DFT')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-7.431364,4)

    def testLi(self):
        solv = SCF(li,method='DFT')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-7.431364,4)

    def testOH(self):
        solv = SCF(oh,method='DFT')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-7.431364,4)

def runsuite(verbose=True):
    #import psyco; psyco.full() # Uncomment to use psyco
    if verbose: verbosity=2
    else: verbosity=1
    # If you want more output, uncomment this line:
    #logging.basicConfig(format="%(message)s",level=logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(ODFTTests)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)
    # Running without verbosity is equivalent to replacing the above
    # two lines with the following:
    #unittest.main()
    return

def debugsuite():
    import cProfile,pstats
    cProfile.run('runsuite()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__':
    import sys
    if "-d" in sys.argv:
        debugsuite()
    else:
        runsuite()
    
    
