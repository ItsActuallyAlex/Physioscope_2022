#### MODULES
import numpy as np
import math
import scipy.optimize as scipy
#### MODULES

#### CONSTANTES QUI VARIENT PAS DU TOUT
viscosity = 1E6
temp_20 = 293
gaz_p = 8.314
#### CONSTANTES QUI VARIENT PAS DU TOUT

#### CONSTANTES QUI VARIENT ENTRE TRAITEMENTS
photosynthese_journaliere = np.array([429859.7489, 641689.2869, 237700.9616], dtype=float)

longueur_entrenoeuds = np.array([2.91047327, 2.91047327, 2.91047327], dtype=float)
rayon_entrenoeuds = np.array([2E-4, 2E-4, 2E-4], dtype=float)
# Valeurs en cm

volumes_puits1 = np.array([1.79192E-06, 1.79192E-06, 1.79192E-06], dtype=float)

# C0_m = np.array([360,300,300], dtype=float)
# Fitté
C0_m = np.array([327.2154541, 327.2154541, 327.2154541], dtype=float)
# Identique

C1_m = np.array([518.3472222, 479.0671296, 479.0671296], dtype=float)
C2_m = np.array([0.069444444, 0.069444444, 0.069444444], dtype=float)

k = np.array([405.197, 397.561, 357.624], dtype=float)
v = np.array([220, 180, 180], dtype=float)

delta_masse_volume = np.array([1.03453E+13, 1.03453E+13, 1.03453E+13], dtype=float)
#### CONSTANTES QUI VARIENT ENTRE TRAITEMENTS

#### FONCTIONS DU MODELE
def PHOTOSYNTHESE_SOURCE () :
    return photosynthese_journaliere
"PAS UTILISE POUR L'INSTANT"

# def RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds) :
    # constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p)

    # resistances_modele = constante_resistance * (longueur_entrenoeuds)/(rayon_entrenoeuds**4)

    # resistances_modele = np.append(resistances_modele, (resistances_modele, resistances_modele))
    
    # return resistances_modele 
    
def RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds) :
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p)


    resistances_modele = constante_resistance * ((longueur_entrenoeuds)/(rayon_entrenoeuds**4))

    print("BBBB", resistances_modele)
    resistances_modele = np.append(resistances_modele, (resistances_modele, resistances_modele))
    print("BBBB", resistances_modele)
    return resistances_modele

R_ = RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds)
print("Longueur R_",len(R_),R_)

def FLUX_CARBONE (longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m) :
    R_ = RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds)

    denominateur_delta = R_[0]*(R_[1] + R_[2]) + R_[1]*R_[2]
    print("DDD", denominateur_delta)
    # print("Longueur dénominateur_delta",len(denominateur_delta), denominateur_delta)
    
    F01 = (R_[2]*(C0_m - C1_m) + R_[0]*(C2_m - C1_m)) / denominateur_delta
    F02 = (R_[1]*(C0_m - C1_m) + R_[0]*(C1_m - C2_m)) / denominateur_delta
    print("CCC", F01)



    FLUX = np.append(F01, F02)

    return FLUX

# def FLUX_CARBONE (longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m) :
#     R_ = RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds)

#     denominateur_delta = R_[0]*(R_[1] + R_[2]) + R_[1]*R_[2]

#     F01 = (R_[2]*(C0_m - C1_m) + R_[0]*(C2_m - C1_m)) / denominateur_delta
#     F02 = (R_[1]*(C0_m - C1_m) + R_[0]*(C1_m - C2_m)) / denominateur_delta

#     FLUX = np.append(F01, F02)

#     return FLUX

# def VOLUME_PUITS (volumes_puits1) :
#     return volumes_puits1 

def RER_PUITS (v, k, C1_m, C2_m) :

    RER_puits1 = (v*C1_m)/(k+C1_m)
    RER_puits2 = (v*C2_m)/(k+C2_m)

    return [RER_puits1, RER_puits2]
"PAS UTILISE POUR L'INSTANT"

def UTILISATION_CARBONE (delta_masse_volume) :

    utilisation_puits1 = delta_masse_volume * volumes_puits1
    utilisation_puits2 = utilisation_puits1/100

    UTILISATIONS = np.append(utilisation_puits1, utilisation_puits2)

    return UTILISATIONS

def A_RESOUDRE (longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m, delta_masse_volume) :

    C1_m_equilibre = C0_m*FLUX_CARBONE(longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m) - UTILISATION_CARBONE(delta_masse_volume)
    C2_m_equilibre = C0_m*FLUX_CARBONE(longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m) - UTILISATION_CARBONE(delta_masse_volume)

    EQUILIBRE_SOLUTIONS_PUITS = np.append(C1_m_equilibre, C2_m_equilibre)

    return EQUILIBRE_SOLUTIONS_PUITS
#### FONCTIONS DU MODELE   


# print("resistances : ", RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds))
# print("flux : ", FLUX_CARBONE (longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m))
# print("utilisation : " ,UTILISATION_CARBONE (delta_masse_volume))


A = np.array([longueur_entrenoeuds])
B = np.array([rayon_entrenoeuds])
C = np.array([C0_m])
X = np.array([C1_m])
Y = np.array([C2_m])
D = np.array([volumes_puits1])
E = np.array([delta_masse_volume])
    
X0 = np.append(X, Y)
args = A,B,C,D,E
    
equilibre_concentrations = scipy.fsolve(A_RESOUDRE, x0=(X, Y), args=(args))

print("FONCTION PASSEE")




# for condition in range(3) :
#     # = np.array([photosynthese_journaliere[condition]])
#     A = np.array([longueur_entrenoeuds[condition]])
#     B = np.array([rayon_entrenoeuds[condition]])
#     # C0_m_t0_Sucrose_fitted = np.array([C0_m_t0_Sucrose_fitted[condition]])
#     C = np.array([C0_m[condition]])
#     X = np.array([C1_m[condition]])
#     Y = np.array([C2_m[condition]])
#     D = np.array([volumes_puits1[condition]])
#     # = np.array([k[condition]])
#     # = np.array([v[condition]])
#     E = np.array([delta_masse_volume[condition]])
    
#     X0 = np.append(X, Y)
#     args = A,B,C,D,E
    
#     equilibre_concentrations = scipy.fsolve(A_RESOUDRE, x0=(X, Y), args=(args))

