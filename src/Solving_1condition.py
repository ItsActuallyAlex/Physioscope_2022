#### MODULES
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES

#### CONSTANTES VRAIMENT CONSTANTES
temp_20 = 293
gaz_p = 8.314
viscosity = 1E6

nom_conditions = np.array(["HH", "LH", "LL"], dtype=str)

delta = 1.379E+10
#### CONSTANTES VRAIMENT CONSTANTES

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

# BASE_valeur_photosynthese = np.array([429859.7489, 641689.2869, 237700.9616], dtype=float)
BASE_longueur_commune_entrenoeuds = np.array([0.02732, 0.02732, 0.02732], dtype=float)
BASE_rayon_commun_entrenoeuds = np.array([35.1E-6, 35.1E-6, 35.1E-6], dtype=float)
# Valeur moyenne de rayon Minchin

# BASE_valeurs_RER_bourgeon = np.array([110, 90, 90], dtype=float)
# BASE_valeurs_RER_feuilles = np.array([110, 90, 90], dtype=float)
BASE_volume_fixe_bourgeon = np.array([8,849E-10, 8,849E-10, 8,849E-10], dtype=float)
BASE_volume_fixe_feuilles = np.array([1.515E-07, 2.096E-07, 3.027E-07], dtype=float)


# BASE_coef_delta_feuilles = np.array([1.379E+10, 1.379E+10, 1.379E+10], dtype=float)
# BASE_coef_delta_bourgeon = np.array([1.379E+10, 1.379E+10, 1.379E+10], dtype=float)
BASE_v1 = np.array([1.993E-02, 1.884E-02, 1.777E-02], dtype=float)
BASE_k1 = np.array([1.069E+12, 1.566E+12, 9.148E+11], dtype=float)
BASE_v2 = np.array([9.617E-11, 9.617E-11, 9.617E-11], dtype=float)
BASE_k2 = np.array([2.099E+12, 2.099E+12, 2.099E+12], dtype=float)


# BASE_vi = np.array([220, 220, 180, 180, 180, 180], dtype=float)
# BASE_ki = np.array([405.197, 405.197,397.561, 397.561, 357.624,  357.624], dtype=float)
# BASE_conditions_ini_C1C2_équivalentCARBONE = np.array([7.775208333, , 7.186006944, , 7.186006944], dtype=float) 

BASE_C_0 = np.array([6.619E+10, 4.723E+11, 5.972E+10], dtype=float)
# BASE_conditions_ini_C1C2 = np.array([518.3472222, 0.069444444, 479.0671296, 0.069444444, 479.0671296, 0.069444444], dtype=float)
BASE_conc_ini_C1 = np.array([5.344E+11, 7.829E+11, 4.574E+11], dtype=float)
BASE_conc_ini_C2 = np.array([1.04937E+12, 1.04937E+12, 1.04937E+12], dtype=float)
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

def FLUX_2puits (Ci_m, C0_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) :

    R_ = RESISTANCES (longueur_commune_entrenoeuds, rayon_commun_entrenoeuds)
    R_0 = R_[0]
    R_1 = R_[1]
    R_2 = R_[2]
    denominateur = R_0*(R_1+R_2)+R_1*R_2

    F01 = ((R_2*(C0_m-Ci_m[0]) + R_0*(Ci_m[1]-Ci_m[0])) / denominateur)
    F02 = ((R_1*(C0_m-Ci_m[1]) + R_0*(Ci_m[0]-Ci_m[1])) / denominateur)

    eq = np.append(F01,F02)

    return eq
# print("Valeurs flux : ", FLUX_2puits (BASE_conditions_ini_C1C2, BASE_C_0, BASE_longueur_commune_entrenoeuds, BASE_rayon_commun_entrenoeuds))

def RER (Ci_m, v_i, k_i) :
    
    eq = (v_i * Ci_m)/(k_i + Ci_m)

    return eq
# print("Valeurs RER : ", RER (BASE_valeurs_RER_bourgeon, BASE_valeurs_RER_feuilles))

def UTILISATION (Ci_m, volume_i, v_i, k_i) :
    
    eq = delta * volume_i * RER (Ci_m, v_i, k_i)

    return eq
# print("Valeurs utilisation : ", UTILISATION (BASE_coef_delta_feuilles, BASE_coef_delta_bourgeon, BASE_volume_fixe_feuilles, BASE_volume_fixe_bourgeon, BASE_valeurs_RER_bourgeon, BASE_valeurs_RER_feuilles))

def A_RESOUDRE (Ci_m, C0_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds, volume_i, v_i, k_i) :

    dCi_m_dt = C0_m*FLUX_2puits (Ci_m, C0_m, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds) - UTILISATION (Ci_m, volume_i, v_i, k_i)

    return dCi_m_dt
#### FONCTIONS

# condition[0] = HH
# condition[1] = LH
# condition[2] = LL

C1_m = np.empty(0)
C2_m = np.empty(0)



print("DEBUT LOOP CONDITION _________________________")
for condition, nom_condition in enumerate(nom_conditions) : 

    longueur_commune_entrenoeuds = BASE_longueur_commune_entrenoeuds[condition]
    # print("longueur_commune_entrenoeuds",longueur_commune_entrenoeuds)
    rayon_commun_entrenoeuds = BASE_rayon_commun_entrenoeuds[condition]
    # print("rayon_commun_entrenoeuds",rayon_commun_entrenoeuds)
    concentration_0 = BASE_C_0[condition]
    # print("concentration_0", concentration_0)
    # coef_delta_feuilles = BASE_coef_delta_feuilles[condition]
    # coef_delta_bourgeon = BASE_coef_delta_bourgeon[condition]
    # coef_delta = np.append(BASE_coef_delta_feuilles[condition], BASE_coef_delta_bourgeon[condition])
    # print("coef_delta",coef_delta)
    # volume_fixe_feuilles = BASE_volume_fixe_feuilles[condition]
    # volume_fixe_bourgeon = BASE_volume_fixe_bourgeon[condition]
    volume_i = np.append(BASE_volume_fixe_bourgeon[condition],BASE_volume_fixe_feuilles[condition])
    # print("volume_i", volume_i)
    v_i = np.append(BASE_v1[condition], BASE_v2[condition])
    # print("v_i",v_i)
    k_i = np.append(BASE_k1[condition], BASE_k2[condition])
    # print("k_i",k_i)
    # valeurs_RER_bourgeon = BASE_valeurs_RER_bourgeon[condition]
    # valeurs_RER_feuilles = BASE_valeurs_RER_feuilles[condition]
    C1_m_t0 = BASE_conc_ini_C1[condition]
    C2_m_t0 = BASE_conc_ini_C2[condition]

    conditions_ini_C1C2 = np.array([C1_m_t0, C2_m_t0])

    # print("AAAAA",conditions_ini_C1C2)
    # increment_C1 = np.array([0, 2, 4], dtype=int)
    # increment_C2 = np.array([1, 3, 5], dtype=int)
    # conditions_ini_C1C2 = np.append(BASE_conditions_ini_C1C2[increment_C1[condition]], BASE_conditions_ini_C1C2[increment_C2[condition]])
    # print("PRINT INCREMENT C1 :", BASE_conditions_ini_C1C2[increment_C1[condition]])
    # print("PRINT INCREMENT C2 :", BASE_conditions_ini_C1C2[increment_C2[condition]])
    # print("AAAAA",conditions_ini_C1C2)


    print("FLUX", FLUX_2puits (conditions_ini_C1C2, concentration_0, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds))
    print("UTILISATION", UTILISATION(conditions_ini_C1C2, volume_i, v_i, k_i))
    print("RER", RER (conditions_ini_C1C2, v_i, k_i))

    dCi_m_dt = scipy.fsolve(A_RESOUDRE, x0=(conditions_ini_C1C2), args=(concentration_0, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds, volume_i, v_i, k_i), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
    print("Les solutions à l'équilibre C1_m et C2_m pour la condition", nom_conditions[condition], "sont :", dCi_m_dt)

    print("FLUX EQUILIBRE", FLUX_2puits (dCi_m_dt, concentration_0, longueur_commune_entrenoeuds, rayon_commun_entrenoeuds))
    print("UTILISATION EQUILIBRE", UTILISATION(dCi_m_dt, volume_i, v_i, k_i))
    print("RER EQUILIBRE", RER (dCi_m_dt, v_i, k_i))
    print("______________")

    C1_m = np.append(C1_m, dCi_m_dt[0])
    C2_m = np.append(C2_m, dCi_m_dt[1])

print("FIN LOOP CONDITION _________________________")




#### PLOTTING
fig,ax = plt.subplots(2)

# ax[0].plot(range(3), BASE_conc_ini_C1, color="black", marker="+")
# ax[1].plot(range(3), BASE_conc_ini_C2, color="black", marker="+")

# ax[0].plot(range(3), BASE_C_0 , color="black", marker="+")
# ax[1].plot(range(3), BASE_C_0 , color="black", marker="+")
ax[0].plot(range(3), C1_m, color="red", marker="+")
ax[1].plot(range(3), C2_m, color="blue", marker="+")

ax[0].set_title("C1_m et C2_m à l'équilibre")
ax[0].set_xlabel("HH LH LL", color="black", fontsize=14)
ax[0].set_ylabel("Concentration en umolC/m3", color="black", fontsize=14)


#### PLOTTING



plt.show()



