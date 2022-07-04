#### MODULES
from pyrsistent import v
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES

#### CONSTANTES VRAIMENT CONSTANTES
temp_20 = 293
gaz_p = 8.314
viscosity = 1E6
hauteur_bourgeon = 2E-3
rayon_bourgeon = 1.25E-3
#### CONSTANTES VRAIMENT CONSTANTES

#### FORMULES
formule_cone = ((math.pi * rayon_bourgeon**2 * hauteur_bourgeon)/(3))
#### FORMULES

#### VARIABLES ENTRE CONDITIONS
BASE_valeur_photosynthese = np.array([429859.7489], dtype=float)
BASE_longueur_commune_entrenoeuds = np.array([2.91047327E-2], dtype=float)
BASE_rayon_commun_entrenoeuds = np.array([5.1E-6], dtype=float)
BASE_concentration_0 = np.array([360], dtype=float)
BASE_valeurs_RER_bourgeon = np.array([110], dtype=float)
BASE_valeurs_RER_feuilles = np.array([110], dtype=float)
BASE_volume_fixe_bourgeon = np.array([formule_cone], dtype=float)
BASE_volume_fixe_feuilles = np.array([1.79192E-06], dtype=float)
BASE_coef_delta_feuilles = np.array([1.03453E+13], dtype=float)
BASE_coef_delta_bourgeon = np.array([1.03453E+13], dtype=float)

conditions_ini_C1C2 = np.array([518.3472222, 0.069444444], dtype=float)
#### VARIABLES ENTRE CONDITIONS

#### FONCTIONS
def PHOTOSYNTHESE_simple (valeur_photosynthese) :

    eq = valeur_photosynthese

    return eq
print("Valeurs photosynthèse : ", PHOTOSYNTHESE_simple (BASE_valeur_photosynthese))

def RESISTANCES (longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) :

    # Particule qui ne varie pas
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    # calcul de la résistance
    eq = constante_resistance * (longueur_commune_entrenoeuds)/(rayon_commun_entrenoeuds**4)

    # 3 fois la même valeur R0 R1 R2
    eq = np.append(eq, (eq,eq))

    return eq
print("Valeurs resistances : ", RESISTANCES (BASE_longueur_commune_entrenoeuds, BASE_rayon_commun_entrenoeuds))

def FLUX_2puits (concentration_0, Ci_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) :

    R_ = RESISTANCES (longueur_commune_entrenoeuds, rayon_commun_entrenoeuds)
    denominateur = R_[0]*(R_[1]+R_[2])+R_[1]*R_[2]

    F01 = concentration_0*((R_[2]*(concentration_0-Ci_m[0]) + R_[0]*(Ci_m[1]-Ci_m[0])) / denominateur)
    F02 = concentration_0*((R_[1]*(concentration_0-Ci_m[1]) + R_[0]*(Ci_m[0]-Ci_m[1])) / denominateur)

    return F01, F02
print("Valeurs flux : ", FLUX_2puits (BASE_concentration_0, conditions_ini_C1C2, BASE_longueur_commune_entrenoeuds, BASE_rayon_commun_entrenoeuds))

def VOLUME (volume_fixe_feuilles, volume_fixe_bourgeon) :

    volume_feuilles = volume_fixe_feuilles
    volume_bourgeon = volume_fixe_bourgeon
    
    return volume_feuilles, volume_bourgeon
print("Valeurs volume : ", VOLUME (BASE_volume_fixe_feuilles, BASE_volume_fixe_bourgeon))

def RER (valeurs_RER_bourgeon, valeurs_RER_feuilles) :
    
    RER1 = valeurs_RER_feuilles
    RER2 = valeurs_RER_bourgeon

    return RER1, RER2
print("Valeurs RER : ", RER (BASE_valeurs_RER_bourgeon, BASE_valeurs_RER_feuilles))

def COEFF_delta (coefficient_delta_feuilles, coefficient_delta_bourgeon) :
    
    delta_feuilles = coefficient_delta_feuilles
    delta_bourgeon = coefficient_delta_bourgeon

    return delta_feuilles, delta_bourgeon
print("Valeurs delta : ", COEFF_delta (BASE_coef_delta_feuilles, BASE_coef_delta_bourgeon))

def UTILISATION (coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, valeurs_RER_bourgeon, valeurs_RER_feuilles) :
    
    U1 = COEFF_delta (coef_delta_feuilles, coef_delta_bourgeon)[0] * VOLUME (volume_fixe_feuilles, volume_fixe_bourgeon)[0] * RER (valeurs_RER_bourgeon, valeurs_RER_feuilles)[0]
    U2 = COEFF_delta (coef_delta_feuilles, coef_delta_bourgeon)[1] * VOLUME (volume_fixe_feuilles, volume_fixe_bourgeon)[1] * RER (valeurs_RER_bourgeon, valeurs_RER_feuilles)[1]

    return U1, U2
print("Valeurs utilisation : ", UTILISATION (BASE_coef_delta_feuilles, BASE_coef_delta_bourgeon, BASE_volume_fixe_feuilles, BASE_volume_fixe_bourgeon, BASE_valeurs_RER_bourgeon, BASE_valeurs_RER_feuilles))

def A_RESOUDRE (concentration_0, Ci_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds, coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, valeurs_RER_bourgeon, valeurs_RER_feuilles) :

    C1_m = FLUX_2puits (concentration_0, Ci_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds)[0] - UTILISATION (coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, valeurs_RER_bourgeon, valeurs_RER_feuilles)[0]
    C2_m = FLUX_2puits (concentration_0, Ci_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds)[1] - UTILISATION (coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, valeurs_RER_bourgeon, valeurs_RER_feuilles)[1]

    return C1_m, C2_m
#### FONCTIONS

# condition[0] = HH
# condition[1] = LH
# condition[2] = LL

for condition in range(1) : 
    longueur_commune_entrenoeuds = BASE_longueur_commune_entrenoeuds[condition]
    rayon_commun_entrenoeuds = BASE_rayon_commun_entrenoeuds[condition]
    valeur_photosynthese = BASE_valeur_photosynthese[condition]
    concentration_0 = BASE_concentration_0[condition]
    coef_delta_feuilles = BASE_coef_delta_feuilles[condition]
    coef_delta_bourgeon = BASE_coef_delta_bourgeon[condition]
    volume_fixe_feuilles = BASE_volume_fixe_feuilles[condition]
    volume_fixe_bourgeon = BASE_volume_fixe_bourgeon[condition]
    valeurs_RER_bourgeon = BASE_valeurs_RER_bourgeon[condition]
    valeurs_RER_feuilles = BASE_valeurs_RER_feuilles[condition]

    # Ci_m = scipy.fsolve(A_RESOUDRE, x0=conditions_ini, args=(concentration_0, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds, coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, valeurs_RER_bourgeon, valeurs_RER_feuilles), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
    # print(Ci_m)