#### MODULES
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES


#### CONSTANTES
# CONSTANTES PHOTOSYNTHESE __________
# P c'est la durée d'éclairement
P = 20
"Durée d'éclairement de 20°Cj"

valeurs_photosynthese = np.array([429859.7489, 641689.2869, 237700.9616], dtype=float)
"valeurs de photosynthèse à tBFV" 
# CONSTANTES PHOTOSYNTHESE __________

# CONSTANTES RESISTANCE __________
temp_20 = 293
gaz_p = 8.314
viscosity = 1E6 # Viscosité selon Bancal et al. 2002

longueur_commune_entrenoeuds = np.array([2.91047327, 2.91047327, 2.91047327], dtype=float)
rayon_commun_entrenoeuds = np.array([0.160215922,0.160215922, 0.160215922], dtype=float)/100
"Unité = cm"
"Possibilité via Minchin et al. 1993 de faire un rayon différent selon l'age des organes (proto et métaphloème) "
# CONSTANTES RESISTANCE __________

# CONSTANTES UTILISATION ET VOLUME __________
volume_fixe_feuilles = np.array([1.79192E-06, 1.79192E-06, 1.79192E-06], dtype=float)
"Volume foliaire"

hauteur = 2E-3
rayon = 1.25E-3
# volume en m3
formule = ((math.pi * rayon**2 * hauteur)/(3))
volume_fixe_bourgeon = np.array([formule, formule, formule], dtype=float)
"bourgeon bourgeon"

U1 = np.array([236419.6269, 244842.0206, 186507.4465], dtype=float)
"Utilisation puits 1 - POUR VALIDER LES RESULTATS"
# CONSTANTES UTILISATION ET VOLUME __________

# CONSTANTES RER __________
k = np.array([405.197, 397.561, 357.624], dtype=float)
k = np.append(k, k)
v = np.array([220, 180, 180], dtype=float)
v = np.append(v, v*100)
"k et v du RER puits 1 et 2"
# Paramètres k et v calculés via cm_t0_puits1

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
"OK - PAS UTILISE MAINTENANT"

def PHOTOSYNTHESE_simple () :
    "valeurs de photosynthèse à tBFV"

    eq = valeurs_photosynthese

    return eq
"OK - PAS UTILISE MAINTENANT"

def RESISTANCES () :
    "Equation de la résistance au transport (Minchin et al. 1993)"

    "R_0[0] = R_0 HH"
    "R_0[1] = R_1 HH"
    "R_0[2] = R_2 HH"
    "R_0[3] = R_0 LH"
    "R_0[4] = R_1 LH"
    "R_0[5] = R_2 LH"
    "R_0[6] = R_0 LL"
    "R_0[7] = R_1 LL"
    "R_0[8] = R_2 LL"

    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    eq = constante_resistance * (longueur_commune_entrenoeuds)/(rayon_commun_entrenoeuds**4)

    eq = np.append(eq, (eq,eq))
    
    return eq
"OK - 3 RESISTANCES PAR CONDITION = 9 RESISTANCES"

def FLUX_2puits (Ci_m) :
    "Equations de F01 et F02"

    "eq[0] = F01 HH"
    "eq[1] = F01 LH"
    "eq[2] = F01 LL"
    "eq[3] = F02 HH"
    "eq[4] = F02 LH"
    "eq[5] = F02 LL"

    R_ = RESISTANCES ()

    # a = pour exploiter R_ et Ci_m correctement
    a = np.array([0,3,6])
    # besoin de déclarer eq pour la suite 
    eq = np.empty(0)
    for i in range(3) :
        i = a[i]
        denominateur = R_[0+i]*(R_[1+i]+R_[2+i]) + R_[1+i]*R_[2+i]

        F01 = Ci_m[0+i]*((R_[2+i]*(Ci_m[0+i]-Ci_m[1+i]) + R_[0+i]*(Ci_m[2+i]-Ci_m[1+i])) / denominateur)

        eq = np.append(eq, F01)

    for i in range(3) :
        i = a[i]
        denominateur = R_[0+i]*(R_[1+i]+R_[2+i]) + R_[1+i]*R_[2+i]

        F02 = Ci_m[0+i]*((R_[1+i]*(Ci_m[0+i]-Ci_m[2+i]) + R_[0+i]*(Ci_m[1+i]-Ci_m[2+i])) / denominateur)

        eq = np.append(eq, F02)

    return eq
"OK - 2 FLUX PAR CONDITION = 6 FLUX "

def FLUX_2puits_UPDATE (Ci_m, C0_m) :
    "Equations de F01 et F02"

    "eq[0] = F01 HH"
    "eq[1] = F01 LH"
    "eq[2] = F01 LL"
    "eq[3] = F02 HH"
    "eq[4] = F02 LH"
    "eq[5] = F02 LL"

    R_ = RESISTANCES ()

    # a = pour exploiter R_ et Ci_m correctement
    a = np.array([0,3,6])
    c = np.array([0,2,4])
    # besoin de déclarer eq pour la suite 
    eq = np.empty(0)
    for i in range(3) :
        b = a[i]
        d = c[i]
        denominateur = R_[0+b]*(R_[1+b]+R_[2+b]) + R_[1+b]*R_[2+b]

        F01 = C0_m*((R_[2+b]*(C0_m-Ci_m[0+d]) + R_[0+b]*(Ci_m[1+d]-Ci_m[0+d])) / denominateur)

        eq = np.append(eq, F01)

    for i in range(3) :
        i = a[i]
        b = i
        denominateur = R_[0+b]*(R_[1+b]+R_[2+b]) + R_[1+b]*R_[2+b]

        F02 = C0_m*((R_[1+b]*(C0_m-Ci_m[1+d]) + R_[0+b]*(Ci_m[0+d]-Ci_m[1+d])) / denominateur)

        eq = np.append(eq, F02)

    return eq
"OK - 2 FLUX PAR CONDITION = 6 FLUX "

def RER () :
    """Expression du RER"""

    eq = np.array([110, 90, 90], dtype=float)
    eq = np.append(eq, eq/100)


    return eq
"OK - 3 RER PAR PUITS = 6 RER - PAS UTILISE MAINTENANT"

def RER_dynamique (Ci_m) :

    "eq[0] = RER_1 HH"
    "eq[1] = RER_1 LH"
    "eq[2] = RER_1 LL"
    "eq[3] = RER_2 HH"
    "eq[4] = RER_2 LH"
    "eq[5] = RER_2 LL"

    a = np.array([0,3,6])
    eq = np.empty(0)

    for i in range(3) :
        n = a[i]
        RERi = (v[i] * Ci_m[1+n])/(k[i] + Ci_m[1+n])
        eq = np.append(eq, RERi)

    for i in range(3) :   
        n = a[i] 
        RERi = (v[i] * Ci_m[2+n])/(k[i] + Ci_m[2+n])
        eq = np.append(eq, RERi)
    
    return eq
"OK- 3 RER PAR PUITS = 6 RER"

def RER_dynamique_UPDATE (Ci_m) :

    "eq[0] = RER_1 HH"
    "eq[1] = RER_1 LH"
    "eq[2] = RER_1 LL"
    "eq[3] = RER_2 HH"
    "eq[4] = RER_2 LH"
    "eq[5] = RER_2 LL"

    a = np.array([0,3,6])
    eq = np.empty(0)

    for i in range(3) :
        n = a[i]
        RERi = (v[i] * Ci_m[0+n])/(k[i] + Ci_m[0+n])
        eq = np.append(eq, RERi)

    for i in range(3) :   
        n = a[i] 
        RERi = (v[i] * Ci_m[1+n])/(k[i] + Ci_m[1+n])
        eq = np.append(eq, RERi)
    
    return eq
"OK- 3 RER PAR PUITS = 6 RER"

def VOLUME () :
    "eq[0] = V_1 HH"
    "eq[1] = V_1 LH"
    "eq[2] = V_1 LL"
    "eq[3] = V_2 HH"
    "eq[4] = V_2 LH"
    "eq[5] = V_2 LL"
    eq = np.append(volume_fixe_feuilles, volume_fixe_bourgeon)
    
    return eq
"OK - 3 VOLUMES POUR CHAQUE PUITS = 6 VOLUMES"

def VOLUME_dynamique (Vi, Ci_m) :
    """Variation du volume des organes, simplifiée"""

    eq = Vi * RER_dynamique (Ci_m)
    
    return eq
"PAS OK - PAS UTILISE MAINTENANT"

def UTILISATION () :

    "eq[0] = U_1 HH"
    "eq[1] = U_1 LH"
    "eq[2] = U_1 LL"
    "eq[3] = U_2 HH"
    "eq[4] = U_2 LH"
    "eq[5] = U_2 LL"

    eq = delta * VOLUME ()

    return eq 
"OK - 3 UTILISATIONS PAR PUITS = 6 UTILISATIONS"

def UTILISATION_dynamique (Vi, Ci_m) :
    "Utilisation avec dynamique de volume"

    eq = delta * VOLUME_dynamique (Vi, Ci_m)

    return eq
"PAS OK - PAS UTILISE MAINTENANT"

def A_RESOUDRE (Ci_m) :
    "Equations à résoudre"

    "Ci_m[0] = C0_m HH"
    "Ci_m[1] = C1_m HH"
    "Ci_m[2] = C2_m HH"
    "Ci_m[3] = C0_m LH"
    "Ci_m[4] = C1_m LH"
    "Ci_m[5] = C2_m LH"
    "Ci_m[6] = C0_m LL"
    "Ci_m[7] = C1_m LL"
    "Ci_m[8] = C2_m LL"

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

def A_RESOUDRE_UPDATE (Ci_m, C0_m) :
    "Equations à résoudre"

    "Ci_m[0] = C1_m HH"
    "Ci_m[1] = C2_m HH"
    "Ci_m[2] = C1_m LH"
    "Ci_m[3] = C2_m LH"
    "Ci_m[4] = C1_m LL"
    "Ci_m[5] = C2_m LL"


    eq = np.empty(0)

    for i in range(3) :
        C1_m = FLUX_2puits_UPDATE (Ci_m, C0_m)[0+i] - UTILISATION ()[0+i]
        eq = np.append(eq, C1_m)
        C2_m = FLUX_2puits_UPDATE (Ci_m, C0_m)[3+i] - UTILISATION ()[3+i]
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

# Ci_m_t0 = np.append(C0_m_t0_Sucrose_fit, (C1_m_t0, C2_m_t0))
# Ci_m_t0 = np.append(C0_m_t0_Sucrose_moyen, (C1_m_t0, C2_m_t0))
Ci_m_t0 = np.append(C1_m_t0, C2_m_t0)
#### CONDITIONS INITIALES


#### RESOLUTION
# La résolution va nous donner les valeurs de Cm_i dans les puits et la source à t+1
# Ci_m = scipy.fsolve(A_RESOUDRE, x0=Ci_m_t0, args=(), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)
Ci_m = scipy.fsolve(A_RESOUDRE_UPDATE, x0=Ci_m_t0, args=(C0_m_t0_Sucrose_moyen), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)
#### RESOLUTION

print(Ci_m)

# C0_m_equilibre = np.empty(0)
# C1_m_equilibre = np.empty(0)
# C2_m_equilibre = np.empty(0)
# a = np.array([0,3,6])
# for i in range(3) :
#     i = a[i]
#     C0_m_equilibre = np.append(C0_m_equilibre, Ci_m[i])
#     C1_m_equilibre = np.append(C1_m_equilibre, Ci_m[1+i])
#     C2_m_equilibre = np.append(C2_m_equilibre, Ci_m[2+i])    

# print(C0_m_equilibre)
# print(C1_m_equilibre)
# print(C2_m_equilibre)

#### PLOTTING
# create figure and axis objects with subplots()
# fig,ax = plt.subplots(2)


# # ax[0].plot(range(3), C0_m_equilibre , color="blue", marker="+")
# ax[0].plot(range(3), C1_m_equilibre, color="red", marker="+")
# # ax[0].plot(range(3), C2_m_equilibre, color="yellow", marker="+")

# ax[0].set_title("Concentrations à l'équilibre")
# ax[0].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0].set_ylabel("Concentration en umolEqGlucose/gMS", color="black", fontsize=14)

# location = 0 # For the best location
# legend_drawn_flag = True
# ax[0].legend(["C0_m", "C1_m", "C2_m"], loc=0, frameon=legend_drawn_flag)

# ax[1,0].set_title("Nom du plot")
# ax[1,0].set_xlabel("nom axe x", color="black", fontsize=14)
# ax[1,0].set_ylabel("nom axe y", color="black", fontsize=14)

# ax[2,0].set_title("Nom du plot")
# ax[2,0].set_xlabel("nom axe x", color="black", fontsize=14)
# ax[2,0].set_ylabel("nom axe y", color="black", fontsize=14)

# plt.show()
#### PLOTTING
