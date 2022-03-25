#!/usr/bin/env python
"""symcone.py

Class definition of SymCone
This is the module of sdpap.

September 2010, Kenta KATO
"""

__all__ = ['SymCone']

class SymCone(object):
    """Description of a Cartesian product of symmetric cones.

    Attributes:
      f: A column vector space to denote free variables or equality constraints
      l: An LP cone
      q: An SOCP cone or a product of SOCP cones
      s: A column vector space corresponding to
           an SDP cone or a product of SDP cones
    """
    def __init__(self, f=0, l=0, q=(), s=()):
        """Constructor

        Args:
          f: Size of free variables or equality constraints
          l: Size LP cone
          q: A tuple of size of SOCP cone
          s: A tuple of size of SDP cone
        """
        if isinstance(f, int):
            self.f = f
        else:
            print("SymCone(): SymCone.f must be integer.\n"
                  "SymCone.f = 0 is set.")
            self.f = 0

        if isinstance(l, int):
            self.l = l
        else:
            print("SymCone(): SymCone.l must be integer.\n"
                  "SymCone.l = 0 is set.")
            self.l = 0

        if isinstance(q, tuple):
            self.q = q
        else:
            print("SymCone(): SymCone.q must be tuple of integers.\n"
                  "SymCone.q = () is set.")
            self.q = ()

        if isinstance(s, tuple):
            self.s = s
        else:
            print("SymCone(): SymCone.s must be tuple of integers\n"
                  "SymCone.s = () is set.")
            self.s = ()

    def __str__(self):
        """Convert into strings for print()."""
        retStr = ""
        retStr += "       f: " + str(self.f) + "\n"
        retStr += "       l: " + str(self.l) + "\n"
        retStr += "  len(q): " + str(len(self.q)) + "\n"
        retStr += " size(q): " + str(sum(self.q)) + "\n"
        if len(self.q) > 0:
            retStr += "range(q): " + \
                      str(min(self.q)) + "-" + str(max(self.q)) + "\n"
        else:
            retStr += "range(q): 0\n"

        retStr += "  len(s): " + str(len(self.s)) + "\n"
        retStr += " size(s): " + str(sum(k ** 2 for k in self.s)) + "\n"
        if len(self.s) > 0:
            retStr += "range(s): " + \
                      str(min(self.s)) + "-" + str(max(self.s)) + "\n"
        else:
            retStr += "range(s): 0\n"
        return retStr

    def debugprint(self):
        """Print SymCone for debug"""
        retStr = str(self)
        retStr += "       q: " + str(self.q) + "\n"
        retStr += "       s: " + str(self.s) + "\n"
        print(retStr)

    def check_validity(self):
        """Check validity for SDPA input

        Returns:
          True: Valid
          False: Invalid
        """
        if not isinstance(self.f, int):
            print("SymCone.check_validity(): SymCone.f must be integer.")
            return False
        if not isinstance(self.l, int):
            print("SymCone.check_validity(): SymCone.l must be integer.")
            return False
        if not isinstance(self.q, tuple):
            print("SymCone.check_validity(): "
                  "SymCone.q must be tuple of integers.")
            return False
        if not isinstance(self.s, tuple):
            print("SymCone.check_validity(): "
                  "SymCone.s must be tuple of integers")
            return False

        return True

    def fromdict(self, otherDict):
        """Convert from canonical dictionary
        See also todict().

        Args:
          otherDict: A dictionary which keys are ('f', 'l', 'q', 's')
        """
        self.f = int(otherDict["f"])
        self.l = int(otherDict["l"])
        self.q = tuple(otherDict["q"])
        self.s = tuple(otherDict["s"])
        return

    def todict(self):
        """Convert to canonical dictionary
        See also fromdict().

        Returns:
          otherDict: A dictionary which keys are ('f', 'l', 'q', 's')
        """
        otherDict = {}
        otherDict["f"] = int(self.f)
        otherDict["l"] = int(self.l)
        otherDict["q"] = tuple(self.q)
        otherDict["s"] = tuple(self.s)
        return otherDict

    def copy(self):
        """Deepcopy of SymCone for SymCone()."""
        import copy
        return copy.deepcopy(self)

