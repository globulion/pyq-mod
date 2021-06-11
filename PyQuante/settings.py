'''
This module contains the settings related to PyQuante, useful to
selecting things as backends. In this way the settings can also be
modifed at runtime.
'''
import sys

try:
    import openbabel
    openbabel_enabled = True
except ImportError:
    openbabel_enabled = False
    print("openbabel not found in path, switching to PyQuante backend", file=sys.stderr)

try:
    import clibint
    libint_enabled = True
except ImportError:
    libint_enabled = False
    print("libint extension not found, switching to normal ERI computation", file=sys.stderr)
