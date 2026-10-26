#!/usr/bin/env python3
"""Encode the selected-screen H'' formula supplied by the theory reviewers.

This is a numeric transcription, not an independent derivation. The parent
task supplied the displayed formula after two independent geometric reviews.
Geometric validity is established by those separate derivations and controls.
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def term(coefficient, **powers):
    return {'coefficient':str(coefficient),'powers':powers}


def build():
    # R=rho*c has been expanded so both branches share only rho,c variables.
    q = [
      term(1,D=1,Mt=1,rho=-2,c=-1),
      term(4,M=1,Mt=1,rho=-2,c=-1),
      term(-1,L=1,Mt=1,rho=-2,c=-1),
      term(-2,M=1,Lt=1,rho=-2,c=-1),
      term(-2,E=1,D=1,rho=-2,c=-2),
      term(-2,E=1,L=1,rho=-2,c=-2),
      term(2,C=1,Ct=1,rho=-2,c=-1),
      term(-2,theta=1,C=2,rho=-2,c=-1),
      term(-2,C=1,Ey=1,rho=-1,c=-1),
      term(2,B=1,C=1,rho=-2,c=-2),
      term(-1,My=1,Ct=1,rho=-1),
      term(1,theta=1,C=1,My=1,rho=-1),
      term(1,My=1,Ey=1),
      term(-1,B=1,My=1,rho=-1,c=-1),
    ]
    linear = [term('3/2',u_x=1,ell=-1,rho=-2,c=-1)]
    total = linear+q
    negate = lambda terms:[{'coefficient':str(-__import__('fractions').Fraction(t['coefficient'])),
                            'powers':t['powers']} for t in terms]
    return {'schema':'HDBLAST_CONSENSUS_SELECTED_SCREEN_MAGNETIC_MAP_V1',
      'status':'NUMERIC_TRANSCRIPTION_OF_SUPPLIED_THEORY_FORMULA',
      'formula_provenance':'Parent supplied the formula after Astra and an independent auditor agreed; geometric proof artifacts are separate.',
      'conventions':{'R':'rho*c','theta':'tanh(t)','parameter_order':'Second parameter derivative with g(eps)=g0+eps*h+eps^2*k/2.',
                     'selected_screen':'Conditional chosen slicing normal U and spatial y normal; no observer-invariance claim.',
                     'full_selected_coefficient':'h=3*u_x/(2*ell*rho^2*c)+4*F/(rho*c); coefficients of Tplus,Tminus are h,-h.'},
      'supplied_grouped_formula': 'F=((D+4*M-L)*Mt-2*M*Lt)/(4*rho)-E*(D+L)/(2*R)+(C/R-My/2)*(R*(Ct-theta*C)/(2*rho)-R*Ey/2+B/2)',
      'components':{'linear':linear,'quadratic':q,'Tplus_total':total,'Tminus_total':negate(total)}}


def main():
    path = HERE/'CONSENSUS_SELECTED_MAGNETIC_MAP.json'
    if path.exists() or path.is_symlink():
        raise ValueError('FRESH_OUTPUT_REQUIRED')
    with path.open('x') as handle:
        json.dump(build(),handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'output':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
