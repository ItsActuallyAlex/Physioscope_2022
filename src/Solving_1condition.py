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

nom_conditions = np.array(["HH", "LH", "LL"], dtype=str)
#### CONSTANTES VRAIMENT CONSTANTES

#### FORMULES
formule_cone = ((math.pi * rayon_bourgeon**2 * hauteur_bourgeon)/(3))
#### FORMULES

#### VARIABLES ENTRE CONDITIONS
# BASE_valeur_photosynthese = np.array([429859.7489], dtype=float)
# BASE_longueur_commune_entrenoeuds = np.array([2.91047327E-2], dtype=float)
# BASE_rayon_commun_entrenoeuds = np.array([5.1E-6], dtype=float)
# # Valeur moyenne de rayon Minchin
# BASE_C_0 = np.array([360], dtype=float)
# BASE_valeurs_RER_bourgeon = np.array([110], dtype=float)
# BASE_valeurs_RER_feuilles = np.array([110], dtype=float)
# BASE_volume_fixe_bourgeon = np.array([formule_cone], dtype=float)
# BASE_volume_fixe_feuilles = np.array([1.79192E-06], dtype=float)
# BASE_coef_delta_feuilles = np.array([1.03453E+13], dtype=float)
# BASE_coef_delta_bourgeon = np.array([1.03453E+13], dtype=float)

# conditions_ini_C1C2 = np.array([518.3472222, 0.069444444], dtype=float)

BASE_valeur_photosynthese = np.array([429859.7489, 641689.2869, 237700.9616], dtype=float)
BASE_longueur_commune_entrenoeuds = np.array([0.02732, 0.02732, 0.02732], dtype=float)
BASE_rayon_commun_entrenoeuds = np.array([35.1E-6, 35.1E-6, 35.1E-6], dtype=float)
# Valeur moyenne de rayon Minchin
BASE_C_0 = np.array([360, 300, 300], dtype=float)
BASE_valeurs_RER_bourgeon = np.array([110, 90, 90], dtype=float)
BASE_valeurs_RER_feuilles = np.array([110, 90, 90], dtype=float)
BASE_volume_fixe_bourgeon = np.array([formule_cone, formule_cone, formule_cone], dtype=float)
BASE_volume_fixe_feuilles = np.array([1.79192E-06, 1.79192E-06, 1.79192E-06], dtype=float)
BASE_coef_delta_feuilles = np.array([1.03453E+13, 1.03453E+13, 1.03453E+13], dtype=float)
BASE_coef_delta_bourgeon = np.array([1.03453E+13, 1.03453E+13, 1.03453E+13], dtype=float)

BASE_k1 = np.array([405.197, 397.561, 357.624], dtype=float)
BASE_v1 = np.array([220, 180, 180], dtype=float)
BASE_k2 = np.array([405.197, 397.561, 357.624], dtype=float)
BASE_v2 = np.array([220, 180, 180], dtype=float)
# BASE_conditions_ini_C1C2_équivalentCARBONE = np.array([7.775208333, , 7.186006944, , 7.186006944], dtype=float) 

# BASE_conditions_ini_C1C2 = np.array([518.3472222, 0.069444444, 479.0671296, 0.069444444, 479.0671296, 0.069444444], dtype=float)
BASE_conc_ini_C1 = np.array([518.3472222, 479.0671296, 479.0671296], dtype=float)
BASE_conc_ini_C2 = np.array([0.069444444, 0.069444444, 0.069444444], dtype=float)
#### VARIABLES ENTRE CONDITIONS

#### FONCTIONS

def RESISTANCES (longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) :

    # Particule qui ne varie pas
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    # calcul de la résistance
    eq = constante_resistance * (longueur_commune_entrenoeuds)/(rayon_commun_entrenoeuds**4)

    # 3 fois la même valeur R0 R1 R2
    eq = np.append(eq, (eq,eq))

    return eq
# print("Valeurs resistances : ", RESISTANCES (BASE_longueur_commune_entrenoeuds, BASE_rayon_commun_entrenoeuds))

def FLUX_2puits (C1_m, C2_m, C0_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) :

    R_ = RESISTANCES (longueur_commune_entrenoeuds, rayon_commun_entrenoeuds)
    R_0 = R_[0]
    R_1 = R_[1]
    R_2 = R_[2]
    denominateur = R_0*(R_1+R_2)+R_1*R_2

    F01 = ((R_2*(C0_m-C1_m) + R_0*(C2_m-C1_m)) / denominateur)
    F02 = ((R_1*(C0_m-C2_m) + R_0*(C1_m-C2_m)) / denominateur)

    eq = np.append(F01,F02)

    return eq
# print("Valeurs flux : ", FLUX_2puits (BASE_conditions_ini_C1C2, BASE_C_0, BASE_longueur_commune_entrenoeuds, BASE_rayon_commun_entrenoeuds))

def RER (C1_m, C2_m, v_1, k_1, v_2, k_2) :
    
    RER1 = (v_1*C1_m)/(k_1+C1_m)
    RER2 = (v_2*C2_m)/(k_2+C2_m)

    eq = np.append([RER1, RER2], dtype=float)

    return eq
# print("Valeurs RER : ", RER (BASE_valeurs_RER_bourgeon, BASE_valeurs_RER_feuilles))

def UTILISATION (coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, C1_m, C2_m, v_1, k_1, v_2, k_2) :
    
    U1 = coef_delta_feuilles * volume_fixe_feuilles * RER (C1_m, C2_m, v_1, k_1, v_2, k_2)[0]
    U2 = coef_delta_bourgeon * volume_fixe_bourgeon * RER (C1_m, C2_m, v_1, k_1, v_2, k_2)[1]

    eq = np.append(U1, U2)

    return eq
# print("Valeurs utilisation : ", UTILISATION (BASE_coef_delta_feuilles, BASE_coef_delta_bourgeon, BASE_volume_fixe_feuilles, BASE_volume_fixe_bourgeon, BASE_valeurs_RER_bourgeon, BASE_valeurs_RER_feuilles))

def A_RESOUDRE (C1_m, C2_m, C0_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds, coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, v_1, k_1, v_2, k_2) :

    dCi_m_dt = C0_m*FLUX_2puits (C1_m, C2_m, C0_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) - UTILISATION (coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, C1_m, C2_m, v_1, k_1, v_2, k_2)

    return dCi_m_dt
#### FONCTIONS

# condition[0] = HH
# condition[1] = LH
# condition[2] = LL

print("DEBUT LOOP CONDITION _________________________")
for condition, nom_condition in enumerate(nom_conditions) : 
    longueur_commune_entrenoeuds = BASE_longueur_commune_entrenoeuds[condition]
    rayon_commun_entrenoeuds = BASE_rayon_commun_entrenoeuds[condition]
    valeur_photosynthese = BASE_valeur_photosynthese[condition]
    concentration_0 = BASE_C_0[condition]
    print("PRINT DE C_0 :", concentration_0)
    coef_delta_feuilles = BASE_coef_delta_feuilles[condition]
    coef_delta_bourgeon = BASE_coef_delta_bourgeon[condition]
    volume_fixe_feuilles = BASE_volume_fixe_feuilles[condition]
    volume_fixe_bourgeon = BASE_volume_fixe_bourgeon[condition]
    valeurs_RER_bourgeon = BASE_valeurs_RER_bourgeon[condition]
    valeurs_RER_feuilles = BASE_valeurs_RER_feuilles[condition]
    coef_delta_feuilles = BASE_coef_delta_feuilles[condition]
    coef_delta_bourgeon = BASE_coef_delta_bourgeon [condition]
    k_1 = BASE_k1[condition]
    v_1 = BASE_v1[condition]
    k_2 = BASE_k2[condition]
    v_2 = BASE_v2[condition]

    C1_m_t0 = BASE_conc_ini_C1[condition]
    C2_m_t0 = BASE_conc_ini_C2[condition]

    conditions_ini_C1C2 = np.array([C1_m_t0, C2_m_t0])

    # increment_C1 = np.array([0, 2, 4], dtype=int)
    # increment_C2 = np.array([1, 3, 5], dtype=int)
    # conditions_ini_C1C2 = np.append(BASE_conditions_ini_C1C2[increment_C1[condition]], BASE_conditions_ini_C1C2[increment_C2[condition]])
    # print("PRINT INCREMENT C1 :", BASE_conditions_ini_C1C2[increment_C1[condition]])
    # print("PRINT INCREMENT C2 :", BASE_conditions_ini_C1C2[increment_C2[condition]])

    dCi_m_dt = scipy.fsolve(A_RESOUDRE, x0=(conditions_ini_C1C2), args=(concentration_0, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds, coef_delta_feuilles, coef_delta_bourgeon, volume_fixe_feuilles, volume_fixe_bourgeon, v_1, k_1, v_2, k_2), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
    print("Les solutions à l'équilibre C1_m et C2_m pour la condition",nom_condition[condition], "sont :", dCi_m_dt)

print("FIN LOOP CONDITION _________________________")




# #### PLOTTING
# fig,ax = plt.subplots(2)


# # ax[0].plot(range(3), C0_m_equilibre , color="blue", marker="+")
# ax[0].plot(range(3), Ci_m[increment_C1], color="red", marker="+")

# ax[0].set_title("C1_m à l'équilibre")
# ax[0].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0].set_ylabel("Concentration en umolEqGlucose/gMS", color="black", fontsize=14)

# # location = 0
# # legend_drawn_flag = True
# # ax[0].legend(["C1_m", "C2_m"], loc=0, frameon=legend_drawn_flag)

# ax[0].plot(range(3), Ci_m[increment_C2], color="yellow", marker="+")
# ax[0].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0].set_ylabel("Concentration en umolEqGlucose/gMS", color="black", fontsize=14)
#     #### PLOTTING



# plt.show()



