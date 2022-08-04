#### MODULES
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES

#### CONSTANTES
temp_20 = 293
gaz_p = 8.314
viscosity = 1E6

nom_conditions = np.array(["HH", "LH", "LL"], dtype=str)

delta = 1.694E+10
#### CONSTANTES

#### VARIABLES ENTRE CONDITIONS
BASE_photosynthese = np.array([6.978E+2, 6.978E+2, 6.978E+2], dtype=float)

BASE_longueur_entrenoeuds = np.array([4.339E-2, 4.339E-2, 4.339E-2], dtype=float)
BASE_rayon_entrenoeuds = np.array([35.1E-6, 35.1E-6, 35.1E-6], dtype=float)

BASE_volume_ini_bourgeon = np.array([8.849E-10, 8.849E-10, 8.849E-10], dtype=float)
BASE_volume_ini_feuilles = np.array([4.422E-07, 4.422E-07, 4.422E-07], dtype=float)

# BASE_v1 = np.array([1.779E-02, 1.012E-02, 1.012E-02], dtype=float)
# BASE_k1 = np.array([1.183E+12, 1.183E+12, 1.183E+12], dtype=float)
# BASE_v2 = np.array([3.510E-2, 3.510E-2, 3.510E-2], dtype=float)
# BASE_k2 = np.array([2.099E+12, 2.099E+12, 2.099E+12], dtype=float)

BASE_v1 = np.array([1.779E-0, 1.012E-1, 1.012E-1], dtype=float)
BASE_k1 = np.array([1.183E+12, 1.183E+12, 1.183E+12], dtype=float)
BASE_v2 = np.array([3.510E-2, 3.510E-2, 3.510E-2], dtype=float)
BASE_k2 = np.array([2.099E+12, 2.099E+12, 2.099E+12], dtype=float)

# BASE_conc_ini_C0 = np.array([5.173E+10, 5.173E+10, 5.173E+10], dtype=float)
# BASE_conc_ini_C1 = np.array([5.916E+11, 5.916E+11, 5.916E+11], dtype=float)
# BASE_conc_ini_C2 = np.array([1.049E+12, 1.049E+12, 1.049E+12], dtype=float)

BASE_conc_ini_C0 = np.array([5.173E+12, 5.173E+12, 5.173E+12], dtype=float)
BASE_conc_ini_C1 = np.array([5.916E+11, 5.916E+11, 5.916E+11], dtype=float)
BASE_conc_ini_C2 = np.array([1.049E+12, 1.049E+12, 1.049E+12], dtype=float)
#### VARIABLES ENTRE CONDITIONS

#### ARRAYS INIT
C1_m = np.empty(0)
C2_m = np.empty(0)

RER1 = np.empty(0)
RER2 = np.empty(0)

VOLUME1 = np.empty(0)
VOLUME2 = np.empty(0)

U1 = np.empty(0)
U2 = np.empty(0)
#### ARRAYS INIT

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

def VOLUME(Ci_m, v_i, k_i, V_t) :

    eq = V_t + (RER (Ci_m, v_i, k_i) * V_t)

    return eq

def UTILISATION (Ci_m, v_i, k_i, V_t) :
    
    eq = delta * V_t * RER (Ci_m, v_i, k_i)

    return eq

def A_RESOUDRE (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds, v_i, k_i, V_t) :

    eq = C0_m*FLUX_2puits (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds) - UTILISATION (Ci_m, v_i, k_i, V_t)

    return eq

def FLUX_puits1 (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds) :

    R_ = RESISTANCES (longueur_entrenoeuds, rayon_entrenoeuds)
    R_0 = R_[0]
    R_1 = R_[1]
    R_2 = R_[2]
    denominateur = R_0*(R_1+R_2)+R_1*R_2

    F01 = ((R_2*(C0_m-Ci_m[0]) + R_0*(Ci_m[1]-Ci_m[0])) / denominateur)

    return F01

def FLUX_puits2 (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds) :

    R_ = RESISTANCES (longueur_entrenoeuds, rayon_entrenoeuds)
    R_0 = R_[0]
    R_1 = R_[1]
    R_2 = R_[2]
    denominateur = R_0*(R_1+R_2)+R_1*R_2

    F02 = ((R_1*(C0_m-Ci_m[1]) + R_0*(Ci_m[0]-Ci_m[1])) / denominateur)

    return F02
#### FONCTIONS

#### A BOUGER 
C1_m_var = np.empty(0)
C2_m_var = np.empty(0)

FLUX01_var = np.empty(0)
FLUX02_var = np.empty(0)

U1_var = np.empty(0)
U2_var = np.empty(0)

RER1_var = np.empty(0)
RER2_var = np.empty(0)
#### A BOUGER 


#### SOLVING
print("DEBUT LOOP CONDITION _________________________")
for condition, nom_condition in enumerate(nom_conditions) : 

    v_i = np.empty(0)
    k_i = np.empty(0)
    volume_ini = np.empty(0)

    longueur_entrenoeuds = BASE_longueur_entrenoeuds[condition]
    rayon_entrenoeuds = BASE_rayon_entrenoeuds[condition]

    v_i = np.append(BASE_v1[condition], BASE_v2[condition])
    k_i = np.append(BASE_k1[condition], BASE_k2[condition])

    C0_m_t0 = BASE_conc_ini_C0[condition]
    C1_m_t0 = BASE_conc_ini_C1[condition]
    C2_m_t0 = BASE_conc_ini_C2[condition]
    conditions_ini_C1C2 = np.array([C1_m_t0, C2_m_t0])

    volume_ini = np.append(BASE_volume_ini_bourgeon[condition],BASE_volume_ini_feuilles[condition])

    V_t = volume_ini
    Ci_m = conditions_ini_C1C2

    print("PRINT FONCTION RER", RER (conditions_ini_C1C2, v_i, k_i))

    print("RER 1 a t0", RER (C1_m_t0, BASE_v1[0], BASE_k1[0]))
    a = (5.9E11*1.779E-1)/(1.183E12+5.9E11)
    b = (4.83E10*1.779E-1)/(1.183E12+4.83E10)
    print("RER1 t0 ", a)
    print("RER1 t1 ", b)
    
    ## BOUCLE UPDATE VOLUME
    for t in range(0, 300, 20) :
        V_t = VOLUME(Ci_m, v_i, k_i, V_t)

        Ci_m = scipy.fsolve(A_RESOUDRE, x0=(Ci_m), args=(C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds, v_i, k_i, V_t), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
        print("Les solutions à l'équilibre C1_m et C2_m pour la condition", nom_conditions[condition], "sont :", Ci_m)
    
        C1_m_var = np.append(C1_m_var, Ci_m[0])
        C2_m_var = np.append(C2_m_var, Ci_m[1])

        FLUX01_var = np.append(FLUX01_var, FLUX_2puits (Ci_m, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[0])
        FLUX02_var = np.append(FLUX02_var, FLUX_2puits (Ci_m, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[1])

        U1_var = np.append(U1_var, UTILISATION (Ci_m, v_i, k_i, V_t)[0])
        U2_var = np.append(U2_var, UTILISATION (Ci_m, v_i, k_i, V_t)[1])

        RER1_var = np.append(RER1_var, RER (Ci_m, v_i, k_i)[0])
        RER2_var = np.append(RER2_var, RER (Ci_m, v_i, k_i)[1])

        print("aaaaa",C0_m_t0)

        ## BOUCLE UPDATE VOLUME

    fig,ax = plt.subplots(2,2)

    # PLOTTING FLUX VAR
    ax[0,0].plot(range(0, 300, 20), FLUX01_var, color="red", marker="+")
    ax[0,0].plot(range(0, 300, 20), FLUX02_var, color="green", marker="+")
    ax[0,0].set_title("FLUX01 et FLUX02 var")
    ax[0,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
     # PLOTTING FLUX VAR

    # PLOTTING UTILISATION VAR
    ax[1,0].plot(range(0, 300, 20), U1_var, color="red", marker="+")
    ax[1,0].plot(range(0, 300, 20), U2_var, color="green", marker="+")
    ax[1,0].set_title("U1 et U2 var")
    ax[1,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
    # PLOTTING UTILISATION VAR

    # PLOTTING RER VAR
    ax[1,1].plot(range(0, 300, 20), RER1_var, color="red", marker="+")
    ax[1,1].plot(range(0, 300, 20), RER2_var, color="green", marker="+")
    ax[1,1].set_title("RER1 et RER2 var")
    ax[1,1].set_ylabel("/°Cj", color="black", fontsize=14)
    # PLOTTING RER VAR

    # PLOTTING CI_M VAR
    ax[0,1].plot(range(0, 300, 20), C1_m_var, color="red", marker="+")
    ax[0,1].plot(range(0, 300, 20), C2_m_var, color="green", marker="+")
    ax[0,1].set_title("C1_m et C2_m var")
    ax[0,1].set_ylabel("umol/m3", color="black", fontsize=14)
    plt.show()
    # PLOTTING CI_M VAR

    # RESET DES VAR
    FLUX01_var = np.empty(0)
    FLUX02_var = np.empty(0)

    U1_var = np.empty(0)
    U2_var = np.empty(0)

    RER1_var = np.empty(0)
    RER2_var = np.empty(0)
    
    C1_m_var = np.empty(0)
    C2_m_var = np.empty(0)
    # RESET DES VAR





    # DONNEES A L'EQUILIBRE
    # C1_m = np.append(C1_m, Ci_m[0])
    # C2_m = np.append(C2_m, Ci_m[1])

    # RER1 = np.append(RER1, RER(C1_m[condition], BASE_v1[condition], BASE_k1[condition])) 
    # RER2 = np.append(RER2, RER(C2_m[condition], BASE_v2[condition], BASE_k2[condition])) 

    # VOLUME1 = np.append(VOLUME1, V_t[0])
    # VOLUME2 = np.append(VOLUME2, V_t[1])

    # U1 = np.append(U1, UTILISATION(C1_m[condition], BASE_v1[condition], BASE_k1[condition], V_t[0]))
    # U2 = np.append(U2, UTILISATION(C2_m[condition], BASE_v2[condition], BASE_k2[condition], V_t[1])) 
    # DONNEES A L'EQUILIBRE


    print("_________________________")
print("FIN LOOP CONDITION _________________________")


# #### PLOTTING
# fig,ax = plt.subplots(2,4)

# # PLOTTING CONCENTRATIONS INI
# ax[0,0].plot(range(3), BASE_conc_ini_C1, color="red", marker="+")
# ax[1,0].plot(range(3), BASE_conc_ini_C2, color="red", marker="+")

# ax[0,0].set_title("C1_m et C2_m")
# ax[0,0].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,0].set_ylabel("Concentration en umolC/m3", color="black", fontsize=14)
# ax[1,0].set_ylabel("Concentration en umolC/m3", color="black", fontsize=14)
# # PLOTTING CONCENTRATIONS INI

# # PLOTTING RER INI
# ax[0,1].plot(range(3), RER(BASE_conc_ini_C1, BASE_v1, BASE_k1), color="red", marker="+")
# ax[1,1].plot(range(3), RER(BASE_conc_ini_C2, BASE_v2, BASE_k2), color="red", marker="+")

# ax[0,1].set_title("RER1 et RER2")
# ax[0,1].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,1].set_ylabel("/°Cj", color="black", fontsize=14)
# ax[1,1].set_ylabel("/°Cj", color="black", fontsize=14)
# # PLOTTING RER INI

# # PLOTTING VOLUME INI
# ax[0,2].plot(range(3), BASE_volume_ini_feuilles, color="red", marker="+")
# ax[1,2].plot(range(3), BASE_volume_ini_bourgeon, color="red", marker="+")

# ax[0,2].set_title("VOLUME 1 et VOLUME 2 ")
# ax[0,2].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,2].set_ylabel("m3", color="black", fontsize=14)
# ax[1,2].set_ylabel("m3", color="black", fontsize=14)
# # PLOTTING VOLUME INI

# # PLOTTING UTILISATION INI
# ax[0,3].plot(range(3), UTILISATION(BASE_conc_ini_C1, BASE_v1, BASE_k1, BASE_volume_ini_feuilles), color="red", marker="+")
# ax[1,3].plot(range(3), UTILISATION(BASE_conc_ini_C2, BASE_v2, BASE_k2, BASE_volume_ini_bourgeon) , color="red", marker="+")

# ax[0,3].set_title("UTLISATION 1 et UTILISATION 2 ")
# ax[0,3].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,3].set_ylabel("umol/°Cj", color="black", fontsize=14)
# ax[1,3].set_ylabel("umol/°Cj", color="black", fontsize=14)
# # PLOTTING UTILISATION INI
# plt.show()






# fig,ax = plt.subplots(2,4)

# # PLOTTING CONCENTRATIONS EQ
# ax[0,0].plot(range(3), C1_m, color="green", marker="+")
# ax[1,0].plot(range(3), C2_m, color="green", marker="+")

# ax[0,0].set_title("C1_m et C2_m")
# ax[0,0].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,0].set_ylabel("Concentration en umolC/m3", color="black", fontsize=14)
# ax[1,0].set_ylabel("Concentration en umolC/m3", color="black", fontsize=14)
# # PLOTTING CONCENTRATIONS EQ

# # PLOTTING RER EQ
# ax[0,1].plot(range(3), RER1, color="green", marker="+")
# ax[1,1].plot(range(3), RER2, color="green", marker="+")

# ax[0,1].set_title("RER1 et RER2")
# ax[0,1].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,1].set_ylabel("/°Cj", color="black", fontsize=14)
# ax[1,1].set_ylabel("/°Cj", color="black", fontsize=14)
# # PLOTTING RER EQ

# # PLOTTING VOLUME EQ
# ax[0,2].plot(range(3), VOLUME1, color="green", marker="+")
# ax[1,2].plot(range(3), VOLUME2, color="green", marker="+")

# ax[0,2].set_title("VOLUME 1 et VOLUME 2 ")
# ax[0,2].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,2].set_ylabel("m3", color="black", fontsize=14)
# ax[1,2].set_ylabel("m3", color="black", fontsize=14)
# # PLOTTING VOLUME EQ

# # PLOTTING UTILISATION EQ
# ax[0,3].plot(range(3), U1, color="green", marker="+")
# ax[1,3].plot(range(3), U2, color="green", marker="+")

# ax[0,3].set_title("UTLISATION 1 et UTILISATION 2 ")
# ax[0,3].set_xlabel("HH LH LL", color="black", fontsize=14)
# ax[0,3].set_ylabel("umol/°Cj", color="black", fontsize=14)
# ax[1,3].set_ylabel("umol/°Cj", color="black", fontsize=14)
# # PLOTTING UTILISATION EQ
# plt.show()

#### PLOTTING





