#### MODULES
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES


#### CONSTANTES
# CONSTANTES PHOTOSYNTHESE __________
# P c'est la durée d'éclairement
P = 16 * 20
"Durée d'éclairement en ??? pour 16h  de photopériode et 20°Cj"

valeurs_photosynthese = np.array([429859.7489, 641689.2869, 237700.9616], dtype=float) 
# CONSTANTES PHOTOSYNTHESE __________

# CONSTANTES RESISTANCE __________
temp_20 = 293
gaz_p = 8.314
# Viscosité selon Bancal et al. 2002 = 1M, on met en micromolaire ici
viscosity = 1E6

longueur_commune_entrenoeuds = np.array([2.91047327, 2.91047327, 2.91047327], dtype=float)
rayon_commun_entrenoeuds = np.array([0.160215922,0.160215922, 0.160215922], dtype=float)
"longueur et rayon moyen d'entrenoeud mature - COMMUNE"
# CONSTANTES RESISTANCE __________

# CONSTANTES UTILISATION ET VOLUME __________
# CHECKER VOIR CE QUI EST CALCULE VIA LA FONCTION UTILISATION (INCREMENTALE DE VOLUME) COMPARER AVEC LA FONCTION QUE J'AI
utilisation = np.array([236419.6269, 244842.0206, 186507.4465], dtype=float)
"Utilisation à tBFV - HH LH LL"

delta = np.array([1.03453E+13, 1.03453E+13, 1.03453E+13], dtype=float)
"moyenne des delta à tBFV - HH LH LL - COMMUNE"

# # On garde un volume fixe pour cette étape
volume_fixe_feuilles = np.array([1.79192E-06, 1.79192E-06, 1.79192E-06], dtype=float)
"Volume foliaire - COMMUNE"

hauteur = 2E-3
rayon = 1.25E-3
# volume en m3
formule = ((math.pi * rayon**2 * hauteur)/(3))
volume_fixe_bourgeon = np.array([formule, formule, formule], dtype=float)
"bourgeon assimilé à un cône - COMMUNE"
# CONSTANTES UTILISATION ET VOLUME __________

# CONSTANTES RER __________
RER = np.array([110, 90, 90], dtype=float)
"Valeurs de RER pour le puits 1 aka les organes en croissance"

k = np.array([405.197, 397.561, 357.624], dtype=float)
v = np.array([220, 180, 180], dtype=float)
"Paramètres k et v calculés via cm_t0_puits1 - HH LH LL"
"Résolu pour v = 2*RER"
# CONSTANTES RER __________
#### CONSTANTES



#### FONCTIONS
def PHOTOSYNTHESE_complete (alpha, PPFD, Ai_max, P) :
    """Equation de la photosynthèse, complète"""

    eq = (alpha * PPFD * Ai_max)/(alpha * PPFD + Ai_max) * P 

    return eq
def CONCENTRATION_SOURCE (surface_fixe_feuilles, volume_fixe_feuilles) :
    "Variation de la concentration dans la source"

    #  voir a nouveau pour les ajouts
    eq = PHOTOSYNTHESE_simple () * (surface_fixe_feuilles)/(volume_fixe_feuilles) 

    return eq


def PHOTOSYNTHESE_simple () :
    "valeurs de photosynthèse à tBFV"

    eq = valeurs_photosynthese

    return eq


def RESISTANCES (viscosity, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) :
    "Equation de la résistance au transport (Minchin et al. 1993)"

    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    eq = constante_resistance * (longueur_commune_entrenoeuds)/(rayon_commun_entrenoeuds**4)

    eq = np.append(eq, (eq,eq))
    "Pas correct pour la suite il faudra corriger l'output pour qu'il s'adapte aux valeurs variables entre conduits"

    return eq


def FLUX_complet (Ci_m) :
    "Equations de F01 et F02"

    R_ = RESISTANCES (viscosity, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds)

    "R_0[0] = R_0 HH"
    "R_0[1] = R_1 HH"
    "R_0[2] = R_2 HH"

    "Ci_m[0] = C0_m HH"
    "Ci_m[1] = C1_m HH"
    "Ci_m[2] = C2_m HH"

    # a c'est un array pour correctement naviguer R_ 
    a = np.array([0,3,6])
    # besoin de déclarer eq pour la suite 
    eq = np.empty(0)
    for i in range(3) :
        i = a[i]
        denominateur = R_[1+i]*(R_[2+i]+R_[2+i])+R_[2+i]*R_[2+i]

        F01 = Ci_m[0]*((R_[2+i]*(Ci_m[0]-Ci_m[1]) + R_[0+i]*(Ci_m[2]-Ci_m[1])) / denominateur)
        F02 = Ci_m[0]*((R_[1+i]*(Ci_m[0]-Ci_m[2]) + R_[0+i]*(Ci_m[1]-Ci_m[2])) / denominateur)

        eq = np.append(eq, (F01,F02))

    return eq


def RERi (v, k, Ci_m) :
    """Expression du RER"""

    eq = (v * Ci_m)/(k + Ci_m)

    return eq


def VOLUME_simple (Vi, Ci_m) :
    """Variation du volume des organes, simplifiée"""

    eq = Vi * RERi (v, k, Ci_m)
    
    return eq


def UTILISATION_simple (delta, Vi) :
    """Equation de l'utilisation, simplifiée"""

    eq = delta * VOLUME_simple (Vi)

    return eq 


def A_RESOUDRE (Ci_m, Vi) :
    """equation à résoudre """

    C0_m = (PHOTOSYNTHESE_simple() * P) / Vi[0] - FLUX_complet (Ci_m)[0] - FLUX_complet (Ci_m)[2]
    C1_m = FLUX_complet (Ci_m)[1] - UTILISATION_simple (delta, Vi)[0]
    C2_m = FLUX_complet (Ci_m)[2] - UTILISATION_simple (delta, Vi)[1]

    return [C0_m, C1_m, C2_m]
#### FONCTIONS


#### CONDITIONS INITIALES
cm_t0_puits1 = np.array([405.197, 397.561, 357.624], dtype=float)
"Concentrations mobiles en C à t0 - HH LH LL"

cm_t0_puits2 = np.array([0, 0, 0], dtype=float)
"Données à extraire de la publication Girault et al. 2008 - HH LH LL"
#### CONDITIONS INITIALES

"CONDITIONS INITIALES"
# On a besoin des concentrations et des utilisations 
Ci_m_t0 = np.array([405.197, 397.561, 357.624, 0, 0, 0], dtype=float)
"Voir comment agencer cet array pour que la solution soit adaptée"


"CONDITIONS INITIALES"







#### RESOLUTION
# La résolution va nous donner les valeurs de Cm_i dans les puits et la source à t+1
Cm_i = scipy.fsolve(A_RESOUDRE, x0=x_0, args=(param), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)
x_result = scipy.fsolve(equa, x0=x_0, args=(Vi), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)

#### RESOLUTION

































































































