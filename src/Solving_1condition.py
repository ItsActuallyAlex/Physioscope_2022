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

C1_m = np.empty(0)
C2_m = np.empty(0)

RER1 = np.empty(0)
RER2 = np.empty(0)

U1 = np.empty(0)
U2 = np.empty(0)

F01 = np.empty(0)
F02 = np.empty(0)
#### CONSTANTES VRAIMENT CONSTANTES

#### VARIABLES ENTRE CONDITIONS
BASE_photosynthese = np.array([6.978E+2, 6.978E+2, 6.978E+2], dtype=float)

BASE_longueur_entrenoeuds = np.array([4.339E-2, 4.339E-2, 4.339E-2], dtype=float)
BASE_rayon_entrenoeuds = np.array([35.1E-6, 35.1E-6, 35.1E-6], dtype=float)

BASE_volume_ini_bourgeon = np.array([8,849E-10, 8,849E-10, 8,849E-10], dtype=float)
BASE_volume_ini_feuilles = np.array([1.063E-06, 1.063E-06, 1.063E-06], dtype=float)

BASE_v1 = np.array([1.779E-02, 1.012E-02, 1.440E-02], dtype=float)
BASE_k1 = np.array([1.183E+12, 1.183E+12, 1.183E+12], dtype=float)
BASE_v2 = np.array([6.510E-2, 4.510E-2, 4.510E-2], dtype=float)
BASE_k2 = np.array([2.099E+12, 2.099E+12, 2.099E+12], dtype=float)

BASE_conc_ini_C0 = np.array([6.619E+10, 4.723E+11, 5.972E+10], dtype=float)
BASE_conc_ini_C1 = np.array([5.344E+11, 7.829E+11, 4.574E+11], dtype=float)
BASE_conc_ini_C2 = np.array([1.049E+12, 1.049E+12, 1.049E+12], dtype=float)
#### VARIABLES ENTRE CONDITIONS

#### FONCTIONS
def RESISTANCES (longueur_entrenoeuds, rayon_entrenoeuds) :

    # Particule qui ne varie pas
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    # calcul de la résistance
    eq = constante_resistance * (longueur_entrenoeuds)/(rayon_entrenoeuds**4)

    # 3 fois la même valeur R0 R1 R2
    eq = np.append(eq, (eq,eq))

    return eq

def FLUX_2puits (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds) :

    R_ = RESISTANCES (longueur_entrenoeuds, rayon_entrenoeuds)
    R_0 = R_[0]
    R_1 = R_[1]
    R_2 = R_[2]
    denominateur = R_0*(R_1+R_2)+R_1*R_2

    F01 = ((R_2*(C0_m-Ci_m[0]) + R_0*(Ci_m[1]-Ci_m[0])) / denominateur)
    F02 = ((R_1*(C0_m-Ci_m[1]) + R_0*(Ci_m[0]-Ci_m[1])) / denominateur)

    eq = np.append(F01,F02)

    return eq

def RER (Ci_m, v_i, k_i) :
    
    eq = (v_i * Ci_m)/(k_i + Ci_m)

    return eq

def UTILISATION (Ci_m, volume_i, v_i, k_i) :
    
    eq = delta * volume_i * RER (Ci_m, v_i, k_i)

    return eq

def A_RESOUDRE (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds, volume_i, v_i, k_i) :

    dCi_m_dt = C0_m*FLUX_2puits (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds) - UTILISATION (Ci_m, volume_i, v_i, k_i)

    return dCi_m_dt
#### FONCTIONS

# #### PRINT VALEURS INI
# print("UTILISATION", UTILISATION (Ci_m, volume_i, v_i, k_i))
# print("FLUX", FLUX_2puits (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds))
# print("RER", RER (Ci_m, v_i, k_i))
# print("VOLUME", )
# #### PRINT VALEURS INI



# print("DEBUT LOOP CONDITION _________________________")
for condition, nom_condition in enumerate(nom_conditions) : 

    longueur_entrenoeuds = BASE_longueur_entrenoeuds[condition]
    rayon_entrenoeuds = BASE_rayon_entrenoeuds[condition]

    volume_ini = np.append(BASE_volume_ini_bourgeon[condition],BASE_volume_ini_feuilles[condition])

    v_i = np.append(BASE_v1[condition], BASE_v2[condition])
    k_i = np.append(BASE_k1[condition], BASE_k2[condition])

    C0_m_t0 = BASE_conc_ini_C0[condition]
    C1_m_t0 = BASE_conc_ini_C1[condition]
    C2_m_t0 = BASE_conc_ini_C2[condition]

    conditions_ini_C1C2 = np.array([C1_m_t0, C2_m_t0])

    # print("__________________________")
    # print("FLUX", FLUX_2puits (conditions_ini_C1C2, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds))
    # print("UTILISATION", UTILISATION(conditions_ini_C1C2, volume_ini, v_i, k_i))
    # print("RER", RER (conditions_ini_C1C2, v_i, k_i))


    dCi_m_dt = scipy.fsolve(A_RESOUDRE, x0=(conditions_ini_C1C2), args=(C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds, volume_ini, v_i, k_i), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
    print("Les solutions à l'équilibre C1_m et C2_m pour la condition", nom_conditions[condition], "sont :", dCi_m_dt)

    # print("FLUX EQUILIBRE", FLUX_2puits (dCi_m_dt, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds))
    # print("UTILISATION EQUILIBRE", UTILISATION(dCi_m_dt, volume_ini, v_i, k_i))
    # print("RER EQUILIBRE", RER (dCi_m_dt, v_i, k_i))
    # print("______________")


    C1_m = np.append(C1_m, dCi_m_dt[0])
    C2_m = np.append(C2_m, dCi_m_dt[1])

    RER1 = np.append(RER1, RER (dCi_m_dt, v_i, k_i)[0])
    RER2 = np.append(RER2, RER (dCi_m_dt, v_i, k_i)[1])

    U1 = np.append(U1,UTILISATION (dCi_m_dt, volume_ini[0], v_i, k_i))
    U2 = np.append(U2, UTILISATION (dCi_m_dt, volume_ini[1], v_i, k_i))

    U1 = np.append(F01, FLUX_2puits (dCi_m_dt, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[0])
    U2 = np.append(F02, FLUX_2puits (dCi_m_dt, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[1])

    print("RER1",RER1)
    print("RER2",RER2)

    print("__________________________")
print("FIN LOOP CONDITION _________________________")




#### PLOTTING
fig,ax = plt.subplots(2,3)

# PLOTTING CONCENTRATIONS
ax[0,0].plot(range(3), C1_m, color="red", marker="+")
ax[1,0].plot(range(3), C2_m, color="blue", marker="+")

ax[0,0].set_title("C1_m et C2_m à l'équilibre")
ax[0,0].set_xlabel("HH LH LL", color="black", fontsize=14)
ax[0,0].set_ylabel("Concentration en umolC/m3", color="black", fontsize=14)
# PLOTTING CONCENTRATIONS

# PLOTTING RER
ax[0,1].plot(range(3), RER1, color="red", marker="+")
ax[1,1].plot(range(3), RER2, color="blue", marker="+")

ax[0,1].set_title("RER1 et RER2 à l'équilibre")
ax[0,1].set_xlabel("HH LH LL", color="black", fontsize=14)
ax[0,1].set_ylabel("RER /°Cj", color="black", fontsize=14)
# PLOTTING RER

# PLOTTING RER
# ax[0,1].plot(range(3), RER1, color="red", marker="+")
# ax[1,1].plot(range(3), RER2, color="blue", marker="+")

# ax[0,2].set_title("U1 et U2 à l'équilibre")
# ax[0,2].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,2].set_ylabel("RER /°Cj", color="black", fontsize=14)
# PLOTTING RER



#### PLOTTING



plt.show()



