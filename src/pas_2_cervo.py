#### MODULES
from pyrsistent import v
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES




valeurs_photosynthese = np.array([429859.7489], dtype=float) 
def PHOTOSYNTHESE_simple () :

    eq = valeurs_photosynthese

    return eq



temp_20 = 293
gaz_p = 8.314
viscosity = 1E6
longueur_commune_entrenoeuds = np.array([2.91047327], dtype=float)
rayon_commun_entrenoeuds = np.array([0.160215922], dtype=float)/100
def RESISTANCES () :

    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    eq = constante_resistance * (longueur_commune_entrenoeuds)/(rayon_commun_entrenoeuds**4)

    eq = np.append(eq, (eq,eq))

    return eq
"PAS DE DIFFERENCE REISSTANCE"

C_0 = 360
def FLUX_2puits (Ci_m) :

    eq = np.empty(0)
    R_ = RESISTANCES ()
    denominateur = R_[0]*(R_[1]+R_[2])+R_[1]*R_[2]

    F01 = C_0*((R_[2]*(C_0-Ci_m[0]) + R_[0]*(Ci_m[1]-Ci_m[0])) / denominateur)
    F02 = C_0*((R_[1]*(C_0-Ci_m[1]) + R_[0]*(Ci_m[0]-Ci_m[1])) / denominateur)

    eq = np.append(eq, (F01,F02))

    return eq


k = np.array([405.197,405.197], dtype=float)
v = np.array([220,220], dtype=float)
def RER (Ci_m) :
    eq = np.empty(0)

    for i in range(2) :
        RERi = (v[i] * Ci_m[0+i])/(k[i] + Ci_m[0+i])
        eq = np.append(eq, RERi)

    b = np.array([eq/10], dtype=float)
    "b c'est les RER pour le puits 2"
    
    return eq
"PAS DE DIFFERENCE RER"

hauteur = 2E-3
rayon = 1.25E-3
formule = ((math.pi * rayon**2 * hauteur)/(3))
volume_fixe_bourgeon = np.array([formule], dtype=float)
volume_fixe_feuilles = np.array([1.79192E-06], dtype=float)
def VOLUME () :

    eq = np.append(volume_fixe_feuilles, volume_fixe_bourgeon)
    
    return eq
"DIFFERENCE VOLUME"

delta = np.array([1.03453E+13, 1.03453E+13], dtype=float)
def UTILISATION () :
    eq = delta * VOLUME ()

    eq = np.append(eq, eq/50)

    return eq 
"DIFFERENCE UTILISATION"


def A_RESOUDRE (Ci_m) :

    eq = np.empty(0)

    C1_m = FLUX_2puits (Ci_m)[0] - UTILISATION ()[0]
    C2_m = FLUX_2puits (Ci_m)[1] - UTILISATION ()[1]

    return [C1_m, C2_m]


conditions_ini = np.array([518.3472222, 0.069444444], dtype=float)


Ci_m = scipy.fsolve(A_RESOUDRE, x0=conditions_ini, args=(), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 