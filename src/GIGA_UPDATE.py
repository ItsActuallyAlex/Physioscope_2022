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

degre_jour = 20
#### CONSTANTES

#### VARIABLES ENTRE CONDITIONS
BASE_longueur_entrenoeuds = np.array([4.339E-2, 4.339E-2, 4.339E-2], dtype=float)
BASE_rayon_entrenoeuds = np.array([35.1E-6, 35.1E-6, 35.1E-6], dtype=float)

BASE_volume_ini_bourgeon = np.array([8.849E-10, 8.849E-10, 8.849E-10], dtype=float)
BASE_volume_ini_feuilles = np.array([4.422E-07, 4.422E-07, 4.422E-07], dtype=float)

# BASES POUR PUITS SIMILAIRES
BASE_v1 = np.array([1.186E-2, 6.747E-3, 6.747E-3], dtype=float)
BASE_k1 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)
BASE_v2 = np.array([2.100E-2, 2.100E-2, 2.100E-2], dtype=float)
BASE_k2 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)

BASE_conc_ini_C0 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)
BASE_conc_ini_C1 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)
BASE_conc_ini_C2 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)
# BASES POUR PUITS SIMILAIRES

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

    F01 = (R_2*(C0_m-Ci_m[0]) + R_0*(Ci_m[1]-Ci_m[0])) / (denominateur)
    F02 = (R_1*(C0_m-Ci_m[1]) + R_0*(Ci_m[0]-Ci_m[1])) / (denominateur)

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
#### FONCTIONS

## SETUP PLOTS
fig,ax = plt.subplots(3,2)
## SETUP PLOTS

#### SOLVING
print("DEBUT LOOP CONDITION _________________________")
for condition, nom_condition in enumerate(nom_conditions) : 
    print("_________________________")
    # RESET SUR LES APPEND ENTRE CONDITIONS
    v_i = np.empty(0)
    k_i = np.empty(0)
    V_t0 = np.empty(0)

    ARRAY_C0_m = np.empty(0)
    ARRAY_C1_m = np.empty(0)
    ARRAY_C2_m = np.empty(0)
    ARRAY_F01 = np.empty(0)
    ARRAY_F02 = np.empty(0)
    ARRAY_RER1 = np.empty(0)
    ARRAY_RER2 = np.empty(0)
    ARRAY_U1 = np.empty(0)
    ARRAY_U2 = np.empty(0)
    ARRAY_volume_1 = np.empty(0)
    ARRAY_volume_2 = np.empty(0)
    # RESET SUR LES APPEND ENTRE CONDITIONS


    longueur_entrenoeuds = BASE_longueur_entrenoeuds[condition]
    rayon_entrenoeuds = BASE_rayon_entrenoeuds[condition]

    v_i = np.append(BASE_v1[condition], BASE_v2[condition])
    k_i = np.append(BASE_k1[condition], BASE_k2[condition])

    C0_m_t0 = BASE_conc_ini_C0[condition]
    C1_m_t0 = BASE_conc_ini_C1[condition]
    C2_m_t0 = BASE_conc_ini_C2[condition]
    C1C2_ini = np.array([C1_m_t0, C2_m_t0])
    Ci_m = C1C2_ini

    V_t0 = np.append(BASE_volume_ini_feuilles[condition], BASE_volume_ini_bourgeon[condition])
    V_t = V_t0

    VARIABLE_VOLUME = V_t0
    VARIABLE_C0_m = C0_m_t0

    ## INCREMENT C0_m
    if nom_condition == "LL"  :
        increment_C0_m = 0.0066 
        
    else :
        increment_C0_m = 0.05
    ## INCREMENT C0_m    

## BOUCLE RESOLUTION TEMPS ____________________________________________________________________________________________________________
    for t in range(0, 300, degre_jour) :

        # RESOLUTION
        Ci_m = scipy.fsolve(A_RESOUDRE, x0=(Ci_m), args=(C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds, v_i, k_i, VARIABLE_VOLUME), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
        print("Les solutions à l'équilibre C1_m et C2_m pour la condition", nom_conditions[condition], "sont :", Ci_m)

        ## APPEND DES ARRAYS A TRACER
        if t == 0 :

            ARRAY_C0_m = np.append(ARRAY_C0_m, C0_m_t0)
            ARRAY_C1_m = np.append(ARRAY_C1_m, C1_m_t0)
            ARRAY_C2_m = np.append(ARRAY_C2_m, C2_m_t0)

            ARRAY_F01 = np.append(ARRAY_F01, FLUX_2puits (C1C2_ini, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[0])
            ARRAY_F02 = np.append(ARRAY_F02, FLUX_2puits (C1C2_ini, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[1])

            ARRAY_RER1 = np.append(ARRAY_RER1, RER (C1C2_ini, v_i, k_i)[0])
            ARRAY_RER2 = np.append(ARRAY_RER2, RER (C1C2_ini, v_i, k_i)[1])

            ARRAY_U1 = np.append(ARRAY_U1, UTILISATION (C1C2_ini, v_i, k_i, V_t)[0])
            ARRAY_U2 = np.append(ARRAY_U2, UTILISATION (C1C2_ini, v_i, k_i, V_t)[1])

            ARRAY_volume_1 = np.append(ARRAY_volume_1, VARIABLE_VOLUME[0])
            ARRAY_volume_2 = np.append(ARRAY_volume_2, VARIABLE_VOLUME[1])

        else :
            ## ACTUALISATION VALEURS
            # INCREMENT CO_m
            VARIABLE_C0_m = VARIABLE_C0_m + VARIABLE_C0_m * increment_C0_m

            # ACTUALISATION VOLUME
            VARIABLE_VOLUME = VOLUME(Ci_m, v_i, k_i, V_t)
            ## ACTUALISATION VALEURS

            ARRAY_C0_m = np.append(ARRAY_C0_m, VARIABLE_C0_m)
            ARRAY_C1_m = np.append(ARRAY_C1_m, Ci_m[0])
            ARRAY_C2_m = np.append(ARRAY_C2_m, Ci_m[1])

            ARRAY_F01 = np.append(ARRAY_F01, FLUX_2puits (Ci_m, VARIABLE_C0_m, longueur_entrenoeuds, rayon_entrenoeuds)[0])
            ARRAY_F02 = np.append(ARRAY_F02, FLUX_2puits (Ci_m, VARIABLE_C0_m, longueur_entrenoeuds, rayon_entrenoeuds)[1])

            ARRAY_RER1 = np.append(ARRAY_RER1, RER (Ci_m, v_i, k_i)[0])
            ARRAY_RER2 = np.append(ARRAY_RER2, RER (Ci_m, v_i, k_i)[1])

            ARRAY_U1 = np.append(ARRAY_U1, UTILISATION (Ci_m, v_i, k_i, V_t)[0])

            ARRAY_U2 = np.append(ARRAY_U2, UTILISATION (Ci_m, v_i, k_i, V_t)[1])

            ARRAY_volume_1 = np.append(ARRAY_volume_1, VARIABLE_VOLUME[0])
            ARRAY_volume_2 = np.append(ARRAY_volume_2, VARIABLE_VOLUME[1])
        ## APPEND DES ARRAYS A TRACER

#### MEGA PLOTTING CRADINGUE ____________________________________________________________________________________________________________
    # IF ON EST EN HH
    if nom_condition == "HH" :
        # PLOTTING FLUX VAR
        FLUX01_HH, = ax[0,0].plot(range(0, 300, degre_jour), ARRAY_F01, color="red", marker="+", label = "puits_1 HH", linestyle="dashdot")
        FLUX02_HH, = ax[0,0].plot(range(0, 300, degre_jour), ARRAY_F02, color="black", marker="+", label = "puits_2 HH", linestyle="dashdot")
        # PLOTTING FLUX VAR

        # PLOTTING UTILISATION VAR
        U1_HH, = ax[1,0].plot(range(0, 300, degre_jour), ARRAY_U1, color="red", marker="+", linestyle="dashdot")
        U2_HH, = ax[1,0].plot(range(0, 300, degre_jour), ARRAY_U2, color="black", marker="+", linestyle="dashdot")
        # PLOTTING UTILISATION VAR

        # PLOTTING RER VAR
        RER1_HH, = ax[1,1].plot(range(0, 300, degre_jour), ARRAY_RER1, color="red", marker="+", linestyle="dashdot")
        RER2_HH, = ax[1,1].plot(range(0, 300, degre_jour), ARRAY_RER2, color="black", marker="+", linestyle="dashdot")
        # PLOTTING RER VAR

        # PLOTTING CI_M VAR
        C0_m_HH, = ax[2,1].plot(range(0, 300, degre_jour), ARRAY_C0_m, color="pink", marker="+", linestyle="dashdot")
        C1_m_HH, = ax[0,1].plot(range(0, 300, degre_jour), ARRAY_C1_m, color="red", marker="+", linestyle="dashdot")
        C2_m_HH, = ax[0,1].plot(range(0, 300, degre_jour), ARRAY_C2_m, color="black", marker="+", linestyle="dashdot")
        # PLOTTING CI_M VAR

        # PLOTTING VOLUME VAR
        V1_HH, = ax[2,0].plot(range(0, 300, degre_jour), ARRAY_volume_1, color="red", marker="+", linestyle="dashdot")
        V2_HH, = ax[2,0].plot(range(0, 300, degre_jour), ARRAY_volume_2, color="black", marker="+", linestyle="dashdot")
        # PLOTTING VOLUME VAR

    # IF ON EST EN LH
    if nom_condition == "LH" :
        
        # PLOTTING FLUX VAR
        FLUX01_LH, = ax[0,0].plot(range(0, 300, degre_jour), ARRAY_F01, color="orange", marker="+", label = "puits_1 LH")
        FLUX02_LH, = ax[0,0].plot(range(0, 300, degre_jour), ARRAY_F02, color="blue", marker="+", label = "puits_2 LH")
        # PLOTTING FLUX VAR

        # PLOTTING UTILISATION VAR
        U1_LH, = ax[1,0].plot(range(0, 300, degre_jour), ARRAY_U1, color="orange", marker="+")
        U2_LH, = ax[1,0].plot(range(0, 300, degre_jour), ARRAY_U2, color="blue", marker="+")
        # PLOTTING UTILISATION VAR

        # PLOTTING RER VAR
        RER1_LH, = ax[1,1].plot(range(0, 300, degre_jour), ARRAY_RER1, color="orange", marker="+")
        RER2_LH, = ax[1,1].plot(range(0, 300, degre_jour), ARRAY_RER2, color="blue", marker="+")
        # PLOTTING RER VAR

        # PLOTTING CI_M VAR
        C0_m_LH, = ax[2,1].plot(range(0, 300, degre_jour), ARRAY_C0_m, color="pink", marker="+", linestyle="dashdot")
        C1_m_LH, = ax[0,1].plot(range(0, 300, degre_jour), ARRAY_C1_m, color="orange", marker="+")
        C2_m_LH, = ax[0,1].plot(range(0, 300, degre_jour), ARRAY_C2_m, color="blue", marker="+")
        # PLOTTING CI_M VAR

        # PLOTTING VOLUME VAR
        V1_LH, = ax[2,0].plot(range(0, 300, degre_jour), ARRAY_volume_1, color="orange", marker="+")
        V2_LH, = ax[2,0].plot(range(0, 300, degre_jour), ARRAY_volume_2, color="blue", marker="+")
        # PLOTTING VOLUME VAR

    # IF ON EST EN LL
    if nom_condition == "LL" :
        # PLOTTING FLUX VAR
        FLUX01_LL, = ax[0,0].plot(range(0, 300, degre_jour), ARRAY_F01, color="yellow", marker="+", label = "puits_1 LL", linestyle="dashed")
        FLUX02_LL, = ax[0,0].plot(range(0, 300, degre_jour), ARRAY_F02, color="cyan", marker="+", label = "puits_2 LL", linestyle="dashed")
        # PLOTTING FLUX VAR

        # PLOTTING UTILISATION VAR
        U1_LL, = ax[1,0].plot(range(0, 300, degre_jour), ARRAY_U1, color="yellow", marker="+", linestyle="dashed")
        U2_LL, = ax[1,0].plot(range(0, 300, degre_jour), ARRAY_U2, color="cyan", marker="+", linestyle="dashed")
        # PLOTTING UTILISATION VAR

        # PLOTTING RER VAR
        RER1_LL, = ax[1,1].plot(range(0, 300, degre_jour), ARRAY_RER1, color="yellow", marker="+", linestyle="dashed")
        RER2_LL, = ax[1,1].plot(range(0, 300, degre_jour), ARRAY_RER2, color="cyan", marker="+", linestyle="dashed")
        # PLOTTING RER VAR

        # PLOTTING CI_M VAR
        C0_m_LL, = ax[2,1].plot(range(0, 300, degre_jour), ARRAY_C0_m, color="pink", marker="+", linestyle="dashdot")
        C1_m_LL, = ax[0,1].plot(range(0, 300, degre_jour), ARRAY_C1_m, color="yellow", marker="+", linestyle="dashed")
        C2_m_LL, = ax[0,1].plot(range(0, 300, degre_jour), ARRAY_C2_m, color="cyan", marker="+", linestyle="dashed")
        # PLOTTING CI_M VAR

        # PLOTTING VOLUME VAR
        V1_LL, = ax[2,0].plot(range(0, 300, degre_jour), ARRAY_volume_1, color="yellow", marker="+", linestyle="dashed")
        V2_LL, = ax[2,0].plot(range(0, 300, degre_jour), ARRAY_volume_2, color="cyan", marker="+", linestyle="dashed")
        # PLOTTING VOLUME VAR 
    # ELSE ERROR SI ON EST NI HH NI LH NI LL
    else : 
        print("ERROR")
#### MEGA PLOTTING CRADINGUE ____________________________________________________________________________________________________________

    print("_________________________")
print("FIN LOOP CONDITION _________________________")

#### PLT SHOW
ax[0,0].set_title("FLUX01 et FLUX02 var")
ax[0,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)

ax[1,0].set_title("U1 et U2 var")
ax[1,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)

ax[1,1].set_title("RER1 et RER2 var")
ax[1,1].set_ylabel("/°Cj", color="black", fontsize=14)

ax[0,1].set_title("C1_m et C2_m var")
ax[0,1].set_ylabel("umol/m3", color="black", fontsize=14)

ax[2,0].set_title("V1 et V2 var")
ax[2,0].set_ylabel("m3", color="black", fontsize=14)

ax[2,1].set_title("C0_m var")
ax[2,1].set_ylabel("umolC/m3", color="black", fontsize=14)

plt.legend(handles=[FLUX01_HH, FLUX01_LH, FLUX01_LL, FLUX02_HH, FLUX02_LH, FLUX02_LL])
plt.show()
#### PLT SHOW

    