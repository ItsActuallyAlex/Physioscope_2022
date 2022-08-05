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


# BASES AJUSTEES FONCTIONNAIT
# BASE_v1 = np.array([1.779E0, 1.012E0, 1.012E0], dtype=float)
# BASE_k1 = np.array([1.183E+12, 1.183E+12, 1.183E+12], dtype=float)
# BASE_v2 = np.array([3.510E0, 3.510E0, 3.510E0], dtype=float)
# BASE_k2 = np.array([2.099E+12, 2.099E+12, 2.099E+12], dtype=float)

# BASE_conc_ini_C0 = np.array([5.173E+12, 5.173E+12, 5.173E+12], dtype=float)
# BASE_conc_ini_C1 = np.array([5.916E+11, 5.916E+11, 5.916E+11], dtype=float)
# BASE_conc_ini_C2 = np.array([1.049E+11, 1.049E+11, 1.049E+11], dtype=float)
# BASES AJUSTEES FONCTIONNAIT

# BASES POUR PUITS SIMILAIRES
BASE_v1 = np.array([1.186E-2, 6.747E-3, 6.747E-3], dtype=float)
BASE_k1 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)
BASE_v2 = np.array([2.100E-2, 2.100E-2, 2.100E-2], dtype=float)
BASE_k2 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)

BASE_conc_ini_C0 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)
BASE_conc_ini_C1 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)
BASE_conc_ini_C2 = np.array([5.642E+11, 5.642E+11, 5.642E+11], dtype=float)
# BASES POUR PUITS SIMILAIRES

# PARAM DE BASE
# BASE_v1 = np.array([1.779E-02, 1.012E-02, 1.012E-02], dtype=float)
# BASE_k1 = np.array([1.183E+12, 1.183E+12, 1.183E+12], dtype=float)
# BASE_v2 = np.array([3.510E-2, 3.510E-2, 3.510E-2], dtype=float)
# BASE_k2 = np.array([2.099E+12, 2.099E+12, 2.099E+12], dtype=float)

# BASE_conc_ini_C0 = np.array([5.173E+10, 5.173E+10, 5.173E+10], dtype=float)
# BASE_conc_ini_C1 = np.array([5.916E+11, 5.916E+11, 5.916E+11], dtype=float)
# BASE_conc_ini_C2 = np.array([1.049E+12, 1.049E+12, 1.049E+12], dtype=float)
# PARAM DE BASE
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

VOLUME1_var = np.empty(0)
VOLUME2_var = np.empty(0)
#### A BOUGER 

## SETUP PLOTS
fig,ax = plt.subplots(3,2)
## SETUP PLOTS

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
    Ci_m = conditions_ini_C1C2

    V_ini = np.append(BASE_volume_ini_bourgeon[condition],BASE_volume_ini_feuilles[condition])
    V_t = V_ini

## BOUCLE RESOLUTION TEMPS ____________________________________________________________________________________________________________
    for t in range(0, 300, degre_jour) :
        V_t = VOLUME(Ci_m, v_i, k_i, V_t)

    #     if condition == "HH" or condition == "LH" :
    #         C0_m = C0_m_HHLH[t] + C0_m_HHLH[t]*0.05


    #     else :
    #         C0_m = C0_m_LL[t]
    # # EVOLUTION C0_m


        Ci_m = scipy.fsolve(A_RESOUDRE, x0=(Ci_m), args=(C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds, v_i, k_i, V_t), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
        print("Les solutions à l'équilibre C1_m et C2_m pour la condition", nom_conditions[condition], "sont :", Ci_m)
        
        # CHECK VALEURS INI INTEGREES
        if t == 0 :
            print("VALEUR RER1 T0", RER (Ci_m, v_i, k_i)[0])
            print("VALEUR RER2 T0", RER (Ci_m, v_i, k_i)[1])
            print("VALEURS C1_m à T0",Ci_m[0])
            print("VALEURS C1_m à T0", Ci_m[1])
            print("VALEURS V1 à T0", V_t[0])
            print("VALEURS V2 à T0", V_t[1])
        # CHECK VALEURS INI INTEGREES



        C1_m_var = np.append(C1_m_var, Ci_m[0])
        C2_m_var = np.append(C2_m_var, Ci_m[1])

        FLUX01_var = np.append(FLUX01_var, FLUX_2puits (Ci_m, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[0])
        FLUX02_var = np.append(FLUX02_var, FLUX_2puits (Ci_m, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[1])

        U1_var = np.append(U1_var, UTILISATION (Ci_m, v_i, k_i, V_t)[0])
        U2_var = np.append(U2_var, UTILISATION (Ci_m, v_i, k_i, V_t)[1])

        RER1_var = np.append(RER1_var, RER (Ci_m, v_i, k_i)[0])
        RER2_var = np.append(RER2_var, RER (Ci_m, v_i, k_i)[1])

        VOLUME1_var = np.append(VOLUME1_var, V_t[0])
        VOLUME2_var = np.append(VOLUME2_var, V_t[1])
## BOUCLE RESOLUTION TEMPS ____________________________________________________________________________________________________________

#### MEGA PLOTTING CRADINGUE ____________________________________________________________________________________________________________
    if nom_condition == "HH" :
        # PLOTTING FLUX VAR
        FLUX01_HH, = ax[0,0].plot(range(0, 300, degre_jour), FLUX01_var, color="red", marker="+", label = "puits_1 HH", linestyle="dashdot")
        FLUX02_HH, = ax[0,0].plot(range(0, 300, degre_jour), FLUX02_var, color="black", marker="+", label = "puits_2 HH", linestyle="dashdot")
        ax[0,0].set_title("FLUX01 et FLUX02 var")
        ax[0,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING FLUX VAR

        # PLOTTING UTILISATION VAR
        U1_HH, = ax[1,0].plot(range(0, 300, degre_jour), U1_var, color="red", marker="+", linestyle="dashdot")
        U2_HH, = ax[1,0].plot(range(0, 300, degre_jour), U2_var, color="black", marker="+", linestyle="dashdot")
        ax[1,0].set_title("U1 et U2 var")
        ax[1,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING UTILISATION VAR

        # PLOTTING RER VAR
        RER1_HH, = ax[1,1].plot(range(0, 300, degre_jour), RER1_var, color="red", marker="+", linestyle="dashdot")
        RER2_HH, = ax[1,1].plot(range(0, 300, degre_jour), RER2_var, color="black", marker="+", linestyle="dashdot")
        ax[1,1].set_title("RER1 et RER2 var")
        ax[1,1].set_ylabel("/°Cj", color="black", fontsize=14)
        # PLOTTING RER VAR

        # PLOTTING CI_M VAR
        C1_m_HH, = ax[0,1].plot(range(0, 300, degre_jour), C1_m_var, color="red", marker="+", linestyle="dashdot")
        C2_m_HH, = ax[0,1].plot(range(0, 300, degre_jour), C2_m_var, color="black", marker="+", linestyle="dashdot")
        ax[0,1].set_title("C1_m et C2_m var")
        ax[0,1].set_ylabel("umol/m3", color="black", fontsize=14)
        # PLOTTING CI_M VAR

        # PLOTTING VOLUME VAR
        V1_HH, = ax[2,0].plot(range(0, 300, degre_jour), VOLUME1_var, color="red", marker="+", linestyle="dashdot")
        V2_HH, = ax[2,0].plot(range(0, 300, degre_jour), VOLUME2_var, color="black", marker="+", linestyle="dashdot")
        ax[2,0].set_title("V1 et V2 var")
        ax[2,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING VOLUME VAR

    if nom_condition == "LH" :
        # PLOTTING FLUX VAR
        FLUX01_LH, = ax[0,0].plot(range(0, 300, degre_jour), FLUX01_var, color="orange", marker="+", label = "puits_1 LH")
        FLUX02_LH, = ax[0,0].plot(range(0, 300, degre_jour), FLUX02_var, color="blue", marker="+", label = "puits_2 LH")
        ax[0,0].set_title("FLUX01 et FLUX02 var")
        ax[0,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING FLUX VAR

        # PLOTTING UTILISATION VAR
        U1_LH, = ax[1,0].plot(range(0, 300, degre_jour), U1_var, color="orange", marker="+")
        U2_LH, = ax[1,0].plot(range(0, 300, degre_jour), U2_var, color="blue", marker="+")
        ax[1,0].set_title("U1 et U2 var")
        ax[1,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING UTILISATION VAR

        # PLOTTING RER VAR
        RER1_LH, = ax[1,1].plot(range(0, 300, degre_jour), RER1_var, color="orange", marker="+")
        RER2_LH, = ax[1,1].plot(range(0, 300, degre_jour), RER2_var, color="blue", marker="+")
        ax[1,1].set_title("RER1 et RER2 var")
        ax[1,1].set_ylabel("/°Cj", color="black", fontsize=14)
        # PLOTTING RER VAR

        # PLOTTING CI_M VAR
        C1_m_LH, = ax[0,1].plot(range(0, 300, degre_jour), C1_m_var, color="orange", marker="+")
        C2_m_LH, = ax[0,1].plot(range(0, 300, degre_jour), C2_m_var, color="blue", marker="+")
        ax[0,1].set_title("C1_m et C2_m var")
        ax[0,1].set_ylabel("umol/m3", color="black", fontsize=14)
        # PLOTTING CI_M VAR

        # PLOTTING VOLUME VAR
        V1_LH, = ax[2,0].plot(range(0, 300, degre_jour), VOLUME1_var, color="orange", marker="+")
        V2_LH, = ax[2,0].plot(range(0, 300, degre_jour), VOLUME2_var, color="blue", marker="+")
        ax[2,0].set_title("V1 et V2 var")
        ax[2,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING VOLUME VAR

    if nom_condition == "LL" :
        # PLOTTING FLUX VAR
        FLUX01_LL, = ax[0,0].plot(range(0, 300, degre_jour), FLUX01_var, color="yellow", marker="+", label = "puits_1 LL", linestyle="dashed")
        FLUX02_LL, = ax[0,0].plot(range(0, 300, degre_jour), FLUX02_var, color="cyan", marker="+", label = "puits_2 LL", linestyle="dashed")
        ax[0,0].set_title("FLUX01 et FLUX02 var")
        ax[0,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING FLUX VAR

        # PLOTTING UTILISATION VAR
        U1_LL, = ax[1,0].plot(range(0, 300, degre_jour), U1_var, color="yellow", marker="+", linestyle="dashed")
        U2_LL, = ax[1,0].plot(range(0, 300, degre_jour), U2_var, color="cyan", marker="+", linestyle="dashed")
        ax[1,0].set_title("U1 et U2 var")
        ax[1,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING UTILISATION VAR

        # PLOTTING RER VAR
        RER1_LL, = ax[1,1].plot(range(0, 300, degre_jour), RER1_var, color="yellow", marker="+", linestyle="dashed")
        RER2_LL, = ax[1,1].plot(range(0, 300, degre_jour), RER2_var, color="cyan", marker="+", linestyle="dashed")
        ax[1,1].set_title("RER1 et RER2 var")
        ax[1,1].set_ylabel("/°Cj", color="black", fontsize=14)
        # PLOTTING RER VAR

        # PLOTTING CI_M VAR
        C1_m_LL, = ax[0,1].plot(range(0, 300, degre_jour), C1_m_var, color="yellow", marker="+", linestyle="dashed")
        C2_m_LL, = ax[0,1].plot(range(0, 300, degre_jour), C2_m_var, color="cyan", marker="+", linestyle="dashed")
        ax[0,1].set_title("C1_m et C2_m var")
        ax[0,1].set_ylabel("umol/m3", color="black", fontsize=14)
        # PLOTTING CI_M VAR

        # PLOTTING VOLUME VAR
        V1_LL, = ax[2,0].plot(range(0, 300, degre_jour), VOLUME1_var, color="yellow", marker="+", linestyle="dashed")
        V2_LL, = ax[2,0].plot(range(0, 300, degre_jour), VOLUME2_var, color="cyan", marker="+", linestyle="dashed")
        ax[2,0].set_title("V1 et V2 var")
        ax[2,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)
        # PLOTTING VOLUME VAR

    else : 
        print("ERROR")

    # RESET DES VAR
    FLUX01_var = np.empty(0)
    FLUX02_var = np.empty(0)

    U1_var = np.empty(0)
    U2_var = np.empty(0)

    RER1_var = np.empty(0)
    RER2_var = np.empty(0)

    C1_m_var = np.empty(0)
    C2_m_var = np.empty(0)

    VOLUME1_var = np.empty(0)
    VOLUME2_var = np.empty(0)
    # RESET DES VAR
#### MEGA PLOTTING CRADINGUE ____________________________________________________________________________________________________________

    print("_________________________")
print("FIN LOOP CONDITION _________________________")

#### PLT SHOW
plt.legend(handles=[FLUX01_HH, FLUX02_HH, FLUX01_LH, FLUX02_LH, FLUX01_LL, FLUX02_LL])
plt.show()
#### PLT SHOW




#### PLOTTING





