"""Exact moment regrouping of independent Pair inventory; no physical evaluations."""
import numpy as np
LD=np.longdouble

def weighted_subtractions(L,moments,jet):
    h0,h1,h2,h3,h4,h5=jet
    m1,m3,m5,m7,m9,m11,m13,m15=moments
    l2=L*L
    l3=l2*L
    l4=l3*L
    l5=l4*L
    l6=l5*L
    l7=l6*L
    l8=l7*L
    l9=l8*L
    l10=l9*L
    l11=l10*L
    l12=l11*L
    l13=l12*L
    l14=l13*L
    l15=l14*L
    l16=l15*L
    b0=((LD(1)/LD(8)))*h3
    b1=l2
    b2=b1*h0*m1
    b3=b1*h2*m3
    b4=l3
    b5=b4*h1*m3
    b6=((LD(1)/LD(4)))*b4*h3*m5
    b7=l4
    b8=((LD(1)/LD(2)))*b7
    b9=h2*m5
    b10=l5
    b11=b10*h1*m5
    b12=b10*m7
    b13=l6
    b14=b13*h2*m7
    b15=l7
    b16=b15*h1*m7
    b17=l8
    b18=b17*h0*m7
    b19=l9
    b20=l10
    b21=((LD(105)/LD(4)))*l11*h1*m11
    b22=l12*h0*m11
    b23=l14*h0*m13
    return -L*b0*m3 + ((LD(1)/LD(2)))*L*h1*m1 - b0*b12 + ((LD(7)/LD(2)))*b11 - (LD(3)/LD(4))*b13*h0*m5 + ((LD(17)/LD(8)))*b14 + ((LD(81)/LD(4)))*b16 + ((LD(7)/LD(4)))*b17*h2*m9 + 15*b18 - (LD(35)/LD(2))*b19*h1*m9 + b2 - (LD(833)/LD(8))*b20*h0*m9 - b21 + ((LD(21)/LD(2)))*b22 + ((LD(1155)/LD(8)))*b23 + ((LD(1)/LD(8)))*b3 + ((LD(3)/LD(2)))*b5 - b6 + b8*b9 + b8*h0*m3,((LD(5005)/LD(4)))*l16*h0*m15 - (LD(385)/LD(2))*l13*h1*m13 + ((LD(1)/LD(6)))*L*h1*m1 - (LD(1)/LD(24))*L*h3*m3 + ((LD(1)/LD(12)))*b1*h4*m5 - (LD(2)/LD(3))*b11 - (LD(11)/LD(8))*b12*h3 + ((LD(9)/LD(4)))*b13*h0*m5 - (LD(69)/LD(8))*b14 - (LD(7)/LD(6))*b15*h3*m9 - (LD(251)/LD(12))*b16 + ((LD(35)/LD(3)))*b17*h2*m9 - (LD(75)/LD(2))*b18 + ((LD(273)/LD(2)))*b19*h1*m9 - (LD(1)/LD(3))*b2 + ((LD(3059)/LD(8)))*b20*h0*m9 + ((LD(77)/LD(4)))*b20*h2*m11 - b21 - (LD(2653)/LD(4))*b22 - (LD(4389)/LD(8))*b23 - (LD(13)/LD(24))*b3 - (LD(1)/LD(3))*b5 - b6 - (LD(17)/LD(12))*b7*b9 + ((LD(1)/LD(6)))*b7*h0*m3 + ((LD(1)/LD(24)))*b7*h4*m7 - (LD(1)/LD(6))*h2*m1 + ((LD(1)/LD(24)))*h4*m3
