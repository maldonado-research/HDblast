"""Change only the momentum integration domain in the sibling native-sector implementation."""
from pathlib import Path
import importlib.util
import math
import os

import numpy as np

base=Path(__file__).resolve().parents[1]/'quantum_shell.py'
spec=importlib.util.spec_from_file_location('base_quantum_shell',base)
native=importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)


class QuantumShell(native.QuantumShell):
    def __init__(self,rb,loga_ref,nodes,feedback,decay):
        super().__init__(rb,loga_ref,nodes,feedback,decay)
        self.edge=float(os.environ['HDBLAST_COMOVING_BAND_EDGE'])
        if self.edge not in (12.,14.):
            raise ValueError('Unknown wider-band experiment')
        x,w=np.polynomial.legendre.leggauss(nodes)
        self.k=.5*self.edge*(1+x)
        self.weight=.5*self.edge*w*self.k**2/(2*np.pi**2)

    def extras(self,m,Ab,phib,AT,phiT,eB,T=None):
        values=super().extras(m,Ab,phib,AT,phiT,eB,T)
        if values:
            a=math.exp(Ab-self.ref)
            values['max_physical_frequency']=math.hypot(self.edge/a,native.G*(phib-native.PHI_STAR))
        return values
