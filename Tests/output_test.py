from PyQuante.hartree_fock import rhf
from PyQuante import SCF,Molecule,configure_output
import unittest

class Stream(object):
    def __init__(self):
        self.written = False
    def write(self,text):
        self.written = True

class OutputTest(unittest.TestCase):
    def setUp(self):
        molecule = Molecule("h2",[(1,(0,0,0)),
                                  (1,(1,0,0))]
                            )
        self.m = molecule
    def runapp(self):
        molecule = self.m
        res = rhf(molecule)
        SCF(molecule).iterate()
    def out(self):
        print
        print "="*40
        print "Must produce output"
        print "="*40
        print
    def no_out(self):
        print
        print "="*40
        print "No output"
        print "="*40
        print
    def testVoid(self):
        configure_output()
        self.runapp()
        self.out()
    def testStreaFalse(self):
        configure_output(stream=None)
        self.runapp()
        self.no_out()
    def testStream(self):
        stream = Stream()
        configure_output(stream=stream)
        self.runapp()
        self.assert_(stream.written)
        self.no_out()
    def testFile(self):
        configure_output(filename="/tmp/example.log")
        self.runapp()
        self.out()
    def testFileStream(self):
        configure_output(filename="/tmp/example1.log",stream=Stream())
        self.runapp()
        self.no_out()

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(OutputTest)
if __name__ == '__main__':
    unittest.main()
