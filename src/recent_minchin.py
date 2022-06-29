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
rayon_commun_entrenoeuds = np.array([0.160215922,0.160215922, 0.160215922], dtype=float)/100
"Rayson divisé par 100 car on a pris les rayons des entrenoeud et pas du phloème"
"longueur et rayon moyen d'entrenoeud mature - COMMUNE"

"Possibilité via Minchin et al. 1993 de faire un rayon différent selon l'age des organes (proto et métaphloème) "

# CONSTANTES RESISTANCE __________

# CONSTANTES UTILISATION ET VOLUME __________
# CHECKER VOIR CE QUI EST CALCULE VIA LA FONCTION UTILISATION (INCREMENTALE DE VOLUME) COMPARER AVEC LA FONCTION QUE J'AI

# # On garde un volume fixe pour cette étape
volume_fixe_feuilles = np.array([1.79192E-06, 1.79192E-06, 1.79192E-06], dtype=float)
"Volume foliaire - COMMUNE"

hauteur = 2E-3
rayon = 1.25E-3
# volume en m3
formule = ((math.pi * rayon**2 * hauteur)/(3))
volume_fixe_bourgeon = np.array([formule, formule, formule], dtype=float)
"bourgeon assimilé à un cône - COMMUNE"

U1 = np.array([236419.6269, 244842.0206, 186507.4465], dtype=float)
"Utilisation à tBFV pour le puis 1 - HH LH LL"
# CONSTANTES UTILISATION ET VOLUME __________

# CONSTANTES RER __________
k = np.array([405.197, 397.561, 357.624], dtype=float)
v = np.array([220, 180, 180], dtype=float)
"Paramètres k et v calculés via cm_t0_puits1 - HH LH LL"
"Résolu pour v = 2*RER"



delta = np.array([1.03453E+13, 1.03453E+13, 1.03453E+13, 1.03453E+13, 1.03453E+13, 1.03453E+13], dtype=float)
"moyenne des delta à tBFV - HH LH LL - COMMUNE"
"Attention le delta ici est identique pour les deux puits"
# CONSTANTES RER __________
#### CONSTANTES


#### FONCTIONS
def PHOTOSYNTHESE_complete (alpha, PPFD, Ai_max, P) :
    """Equation de la photosynthèse, complète"""

    eq = (alpha * PPFD * Ai_max)/(alpha * PPFD + Ai_max) * P 

    return eq
"OK"

def PHOTOSYNTHESE_simple () :
    "valeurs de photosynthèse à tBFV"

    eq = valeurs_photosynthese

    return eq
"OK"

def RESISTANCES (longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) :
    "Equation de la résistance au transport (Minchin et al. 1993)"

    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    eq = constante_resistance * (longueur_commune_entrenoeuds)/(rayon_commun_entrenoeuds**4)

    eq = np.append(eq, (eq,eq))
    
    return eq
"OK"
"Pas correct pour la suite, il faudra corriger l'output pour qu'il s'adapte aux valeurs variables entre conduits"
"Aussi il faudra conserver l'ordre établi de l'output ou alors changer la fonction FLUX aussi"

def FLUX_2puits (Ci_m) :
    "Equations de F01 et F02"

    R_ = RESISTANCES (longueur_commune_entrenoeuds, rayon_commun_entrenoeuds)

    "R_0[0] = R_0 HH"
    "R_0[1] = R_1 HH"
    "R_0[2] = R_2 HH"

    "Ci_m[0] = C0_m HH"
    "Ci_m[1] = C1_m HH"
    "Ci_m[2] = C2_m HH"

    # a = pour exploiter R_ et Ci_m correctement
    a = np.array([0,3,6])
    # besoin de déclarer eq pour la suite 
    eq = np.empty(0)
    for i in range(3) :
        i = a[i]
        denominateur = R_[1+i]*(R_[2+i]+R_[2+i])+R_[2+i]*R_[2+i]

        F01 = Ci_m[0+i]*((R_[2+i]*(Ci_m[0+i]-Ci_m[1+i]) + R_[0+i]*(Ci_m[2+i]-Ci_m[1+i])) / denominateur)

        eq = np.append(eq, F01)

    for i in range(3) :
        i = a[i]
        denominateur = R_[1+i]*(R_[2+i]+R_[2+i])+R_[2+i]*R_[2+i]

        F02 = Ci_m[0+i]*((R_[1+i]*(Ci_m[0+i]-Ci_m[2+i]) + R_[0+i]*(Ci_m[1+i]-Ci_m[2+i])) / denominateur)

        eq = np.append(eq, F02)

    return eq
"OK"

def FLUX (Ci_m) :
    "Expression du flux une fois la résolution effectuée par les équations cradingues"

    R_ = RESISTANCES (longueur_commune_entrenoeuds, rayon_commun_entrenoeuds)

    a = np.array([0,3,6])
        # besoin de déclarer eq pour la suite 
    eq = np.empty(0)
    for i in range(3) :
        i = a[i]

        F01 = (Ci_m[0+i](Ci_m[0+i] - Ci_m[1+i]))/ (R_[0+i] + R_[1+i])
        F02 = (Ci_m[0+i](Ci_m[1+i] - Ci_m[2+i]))/ (R_[0+i] + R_[2+i])

        eq = np.append(eq,(F01,F02))
    return eq
"PAS OK car cette formule ne marche que pour un puits"

def RER () :
    """Expression du RER"""

    eq = np.array([110, 90, 90], dtype=float)
    "Valeurs de RER pour le puits 1 aka les organes en croissance"

    return eq
"OK"

def RER_dynamique (v, k, Ci_m) :
    "Expression du RER pour le puits 1 qui est le puits en croissance"

    a = np.array([0,3,6])
    eq = np.empty(0)

    for i in range(3) :
        i = a[i]
        RERi = (v[i] * Ci_m[1+i])/(k[i] + Ci_m[1+i])

    eq = np.append(eq, RERi)

    b = np.array([eq/10], dtype=float)
    "b c'est les RER pour le puits 2"
    
    return eq
"OK"

def VOLUME () :
    "Volume fixé"

    eq = np.append(volume_fixe_feuilles, volume_fixe_bourgeon)
    
    return eq
"OK"

def VOLUME_dynamique (Vi, v, k, Ci_m) :
    """Variation du volume des organes, simplifiée"""

    eq = Vi * RER_dynamique (v, k, Ci_m)
    
    return eq
"PAS OK"

def UTILISATION () :
    "Utilisation sans variation de volume prise en compte"

    "Méthode numérique pour le calcul de l'utilisation"
    eq = delta * VOLUME ()

    eq = np.append(eq, eq/50)

    return eq 
"OK"

def UTILISATION_dynamique (delta, Vi, v, k, Ci_m) :
    "Utilisation avec dynamique de volume"

    eq = delta * VOLUME_dynamique (Vi, v, k, Ci_m)

    return eq
"PAS OK manque de formalisation"

def A_RESOUDRE (Ci_m) :
    "Equations à résoudre"

    eq = np.empty(0)

    for i in range(3) :
        # C0_m = (PHOTOSYNTHESE_simple ()[0+i]) / VOLUME ()[0+i] + (RESOLUTION_Flux (Ci_m)[0+i] + RESOLUTION_Flux (Ci_m)[3+i])
        C0_m = Ci_m_t0[0+i]
        eq = np.append(eq, C0_m)
        C1_m = FLUX_2puits (Ci_m)[0+i] - UTILISATION ()[0+i]
        eq = np.append(eq, C1_m)
        C2_m = FLUX_2puits (Ci_m)[3+i] - UTILISATION ()[3+i]
        eq = np.append(eq, C2_m)

    return eq 
#### FONCTIONS


#### CONDITIONS INITIALES
C0_m_t0_Sucrose_fit = np.array([360,300,300], dtype=float)
"Valeurs t0 sucrose fittées autour de la moyenne pour séparer les conditions"
C0_m_t0_Sucrose_moyen = np.array([327.2154541, 327.2154541, 327.2154541], dtype=float)
"Valeur moyenne du sucrose identique dans toutes les conditions"

C1_m_t0 = np.array([518.3472222, 479.0671296, 479.0671296], dtype=float)
"valeurs en hexoses totaux dans la Z1 des entrenoeuds+feuilles"

C2_m_t0 = np.array([0.069444444, 0.069444444, 0.069444444], dtype=float)
"Concentrations de Girault et al. 2008, saccharose + glucose"

Ci_m_t0 = np.append(C0_m_t0_Sucrose_moyen, (C1_m_t0, C2_m_t0))
#### CONDITIONS INITIALES


#### RESOLUTION
# La résolution va nous donner les valeurs de Cm_i dans les puits et la source à t+1
Ci_m = scipy.fsolve(A_RESOUDRE, x0=Ci_m_t0, args=(), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)
#### RESOLUTION

print(Ci_m)

#### PLOTTING
# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.plot(range_C_0, valeurs_F1, color="blue", marker="+")
ax.plot(range_C_0, valeurs_F2, color="green", marker="+")
# set x-axis label
ax.set_xlabel("Concentration C_0", fontsize=14)
# set y-axis label
ax.set_ylabel("Flux F1 et F2", color="black", fontsize=14)

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(range_C_0, F1F2, color="black", linestyle= "-.")
ax2.set_ylabel("F1/F2", color="black", fontsize=14)
ax2.set_ylim([0, 1])
plt.show()
#### PLOTTING