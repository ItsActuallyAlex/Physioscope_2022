#### MODULES
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES

#### DONNEES EXPERIMETALES
# EXP_C1m_HH = np.array([2.773E+08, np.nan, 1.232E+06], dtype=float)
# EXP_C1m_LH = np.array([2.773E+08, 1.011E+06, 1.080E+06], dtype=float)
# EXP_C1m_LL = np.array([2.773E+08, 3.954E+05, 5.763E+05], dtype=float)

# EXP_C2m_HH = np.array([5.358E+03, 2.548E+00, np.nan], dtype=float)
# EXP_C2m_LH = np.array([5.358E+03, 5.247E+00, np.nan], dtype=float)
# EXP_C2m_LL = np.array([5.358E+03, 8.631E-01, np.nan], dtype=float)

EXP_RER1_HH = np.array([6.644E-03, 6.281E-03, 5.923E-03, 5.570E-03, 5.223E-03, 4.882E-03, 4.548E-03, 4.221E-03, 3.903E-03, 3.594E-03, 3.296E-03, 3.009E-03, 2.735E-03, 2.474E-03, 2.229E-03, 1.998E-03, 1.784E-03, 1.586E-03, 1.404E-03, 1.238E-03, 1.087E-03, 9.519E-04, 8.307E-04, 7.230E-04, 6.275E-04, 5.434E-04, 4.696E-04, 4.051E-04, 3.488E-04, 3.000E-04], dtype=float)
EXP_RER1_LH = np.array([3.803E-03, 3.576E-03, 3.353E-03, 3.133E-03, 2.916E-03, 2.705E-03, 2.498E-03, 2.297E-03, 2.104E-03, 1.917E-03, 1.739E-03, 1.571E-03, 1.412E-03, 1.263E-03, 1.125E-03, 9.983E-04, 8.820E-04, 7.763E-04, 6.808E-04, 5.951E-04, 5.187E-04, 4.508E-04, 3.908E-04, 3.381E-04, 2.919E-04, 2.516E-04, 2.165E-04, 1.861E-04, 1.597E-04, 1.370E-04], dtype=float)
EXP_RER1_LL = np.array([4.056E-03, 3.825E-03, 3.596E-03, 3.370E-03, 3.148E-03, 2.929E-03, 2.715E-03, 2.506E-03, 2.303E-03, 2.107E-03, 1.919E-03, 1.739E-03, 1.569E-03, 1.409E-03, 1.259E-03, 1.120E-03, 9.927E-04, 8.761E-04, 7.702E-04, 6.748E-04, 5.893E-04, 5.131E-04, 4.455E-04, 3.859E-04, 3.336E-04, 2.879E-04, 2.480E-04, 2.133E-04, 1.832E-04, 1.572E-04], dtype=float)

EXP_U1_HH = np.array([2.556E+02, 2.573E+02, 2.574E+02, 2.559E+02, 2.526E+02, 2.477E+02, 2.413E+02, 2.333E+02, 2.240E+02, 2.136E+02, 2.021E+02, 1.899E+02, 1.772E+02, 1.641E+02, 1.509E+02, 1.379E+02, 1.252E+02, 1.129E+02, 1.013E+02, 9.036E+01, 8.019E+01, 7.083E+01, 6.230E+01, 5.458E+01, 4.766E+01, 4.148E+01, 3.600E+01, 3.117E+01, 2.693E+01, 2.322E+00], dtype=float)
EXP_U1_LH = np.array([1.409E+02, 1.442E+02, 1.468E+02, 1.488E+02, 1.502E+02, 1.508E+02, 1.508E+02, 1.502E+02, 1.490E+02, 1.472E+02, 1.450E+02, 1.423E+02, 1.394E+02, 1.362E+02, 1.328E+02, 1.294E+02, 1.260E+02, 1.227E+02, 1.194E+02, 1.163E+02, 1.134E+02, 1.106E+02, 1.081E+02, 1.038E+02, 1.805E+01, 1.561E+01, 1.347E+01, 1.161E+01, 9.985E+00, 8.578E+00], dtype=float)
EXP_U1_LL = np.array([1.213E+02, 1.203E+02, 1.186E+02, 1.163E+02, 1.134E+02, 1.099E+02, 1.058E+02, 1.012E+02, 9.611E+01, 9.068E+01, 8.496E+01, 7.905E+01, 7.304E+01, 6.702E+01, 6.110E+01, 5.533E+01, 4.981E+01, 4.458E+01, 3.968E+01, 3.515E+01, 3.099E+01, 2.721E+01, 2.380E+01, 2.075E+01, 1.804E+01, 1.564E+01, 1.353E+01, 1.167E+01, 1.006E+01, 8.654E+00], dtype=float)

EXP_V1_HH = np.array([1.28654E-06, 1.37811E-06, 1.47047E-06, 1.56305E-06, 1.65526E-06, 1.74648E-06, 1.83612E-06, 1.9236E-06, 2.00837E-06, 2.08994E-06, 2.16785E-06, 2.24173E-06, 2.31127E-06, 2.37626E-06, 2.43655E-06, 2.4921E-06, 2.54291E-06, 2.5891E-06, 2.63082E-06, 2.66827E-06, 2.70171E-06, 2.7314E-06, 2.75765E-06, 2.78075E-06, 2.801E-06, 2.81869E-06, 2.83409E-06, 2.84746E-06, 2.85904E-06, 2.86905E-06])
EXP_V1_LH = np.array([1.78383E-06, 1.85435E-06, 1.92312E-06, 1.98983E-06, 2.05418E-06, 2.11589E-06, 2.17471E-06, 2.23043E-06, 2.28288E-06, 2.33193E-06, 2.37752E-06, 2.4196E-06, 2.45821E-06, 2.49341E-06, 2.52531E-06, 2.55405E-06, 2.5798E-06, 2.60276E-06, 2.62312E-06, 2.6411E-06, 2.65691E-06, 2.67077E-06, 2.68286E-06, 2.69339E-06, 2.70253E-06, 2.71044E-06, 2.71728E-06, 2.72317E-06, 2.72825E-06, 2.73261E-06])
EXP_V1_LL = np.array([2.30537E-06, 2.40283E-06, 2.49838E-06, 2.59157E-06, 2.68196E-06, 2.76913E-06, 2.85269E-06, 2.9323E-06, 3.00768E-06, 3.07858E-06, 3.14484E-06, 3.20636E-06, 3.2631E-06, 3.31511E-06, 3.36248E-06, 3.40535E-06, 3.44394E-06, 3.47847E-06, 3.50922E-06, 3.53646E-06, 3.56048E-06, 3.58159E-06, 3.60006E-06, 3.61617E-06, 3.63018E-06, 3.64233E-06, 3.65284E-06, 3.66192E-06, 3.66975E-06, 3.67649E-06])
DATES = np.array([0, 130, 240])




# ax[1,1].set_xticklabels(range(0, 300, 10), fontsize= 7)
# ax[2,1].set_xticklabels(range(0, 300, 10), fontsize= 7)
# ax[3,1].set_xticklabels(range(0, 300, 10), fontsize= 7)
# ax[4,1].set_xticklabels(range(0, 300, 10), fontsize= 7)
# ax[5,1].set_xticklabels(range(0, 300, 10), fontsize= 7)

# ax[1,1].set_yticklabels(range(0, 300, 10), fontsize= 7)
# ax[2,1].set_yticklabels(range(0, 300, 10), fontsize= 7)
# ax[3,1].set_yticklabels(range(0, 300, 10), fontsize= 7)
# ax[4,1].set_yticklabels(range(0, 300, 10), fontsize= 7)
# ax[5,1].set_yticklabels(range(0, 300, 10), fontsize= 7)



# ax[0, 1].axis('off')
# ax[1, 1].axis('off')
# ax[2, 1].axis('off')
#### DONNEES EXPERIMETALES

#### CONSTANTES
temp_20 = 293
gaz_p = 8.314
viscosity = 1E6
nom_conditions = np.array(["HH", "LH", "LL"], dtype=str)
delta = 1.694E+10
degre_jour = 20

HH_puits1 = "#39f559"
LH_puits1 = "#218E34"
LL_puits1 = "#0C3514"

HH_puits2 = "#39e7f1"
LH_puits2 = "#249096"
LL_puits2 = "#114346"

HH_source0 = "#F7ed35"
LH_source0 = "#A49E23"
LL_source0 = "#4C4910"
#### CONSTANTES

#### VARIABLES ENTRE CONDITIONS
BASE_longueur_entrenoeuds = np.array([4.339E-2, 4.339E-2, 4.339E-2], dtype=float)
BASE_rayon_entrenoeuds = np.array([35.1E-6, 35.1E-6, 35.1E-6], dtype=float)

BASE_volume_ini_bourgeon = np.array([8.849E-10, 8.849E-10, 8.849E-10], dtype=float)
BASE_volume_ini_feuilles = np.array([2.12668E-07, 2.12668E-07, 2.12668E-07], dtype=float)

# Ci_m_t0 identiques
BASE_v1 = np.array([1.779E-02, 1.329E-02, 1.329E-02], dtype=float)
BASE_k1 = np.array([1.387E+08, 1.387E+08, 1.387E+08], dtype=float)
BASE_v2 = np.array([2.100E-02, 2.100E-02, 2.100E-02], dtype=float)
BASE_k2 = np.array([1.387E+08, 1.387E+08, 1.387E+08], dtype=float)

# BASE_v1 = np.array([1.329E-02, 1.329E-02, 1.329E-02], dtype=float)
# BASE_k1 = np.array([1.387E+08, 1.387E+08, 1.387E+08], dtype=float)
# BASE_v2 = np.array([1.329E-02, 1.329E-02, 1.329E-02], dtype=float)
# BASE_k2 = np.array([1.387E+08, 1.387E+08, 1.387E+08], dtype=float)

BASE_conc_ini_C0 = np.array([1.525E+09, 1.525E+09, 1.525E+09], dtype=float)
BASE_conc_ini_C1 = np.array([2.773E+08, 2.773E+08, 2.773E+08], dtype=float)
BASE_conc_ini_C2 = np.array([5.358E+03, 5.358E+03, 5.358E+03], dtype=float)

# BASE_conc_ini_C0 = np.array([1.525E+09, 1.525E+09, 1.525E+09], dtype=float)
# BASE_conc_ini_C1 = np.array([2.773E+08, 2.773E+08, 2.773E+08], dtype=float)
# BASE_conc_ini_C2 = np.array([2.773E+08, 2.773E+08, 2.773E+08], dtype=float)

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
def PLOT_C0m (données_C0m) : 
    if nom_condition == "HH" :
        ax[0,0].plot(np.linspace(0,300, 16), données_C0m, color=HH_source0, marker="+", label = "source_0_HH")
    else :
        if nom_condition == "LH" :
            ax[0,0].plot(np.linspace(0,300, 16), données_C0m, color=LH_source0, marker="+", label = "source_0_LH")
        else :
            ax[0,0].plot(np.linspace(0,300, 16), données_C0m, color=LL_source0, marker="+", label = "source_0_LL")

    # ax[0,0].set_xticklabels(np.linspace(0,300, 16), fontsize= 7)
    # ax[0,0].set_yticklabels(np.linspace(0,300, 16), fontsize= 7)
    # ax[0,0].legend(fontsize= 7)

def PLOT_FLUX (données_F01, données_F02) : 
    if nom_condition == "HH" :
        ax[1,0].plot(np.linspace(0,300, 16), données_F01, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[1,1].plot(np.linspace(0,300, 16), données_F02, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[1,0].plot(np.linspace(0,300, 16), données_F01, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[1,1].plot(np.linspace(0,300, 16), données_F02, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[1,0].plot(np.linspace(0,300, 16), données_F01, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[1,1].plot(np.linspace(0,300, 16), données_F02, color=LL_puits2, marker="+", label = "puits_2_LL")
            # print("a")

    # ax[1,0].set_xticklabels(np.linspace(0,300, 16), fontsize= 7)
    # ax[1,0].set_yticklabels(np.linspace(0,300, 16), fontsize= 7)

def PLOT_Cim (données_C1m, données_C2m) : 
    if nom_condition == "HH" :
        ax[2,0].plot(np.linspace(0,300, 16), données_C1m, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[2,1].plot(np.linspace(0,300, 16), données_C2m, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[2,0].plot(np.linspace(0,300, 16), données_C1m, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[2,1].plot(np.linspace(0,300, 16), données_C2m, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[2,0].plot(np.linspace(0,300, 16), données_C1m, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[2,1].plot(np.linspace(0,300, 16), données_C2m, color=LL_puits2, marker="+", label = "puits_2_LL")

    # ax[2,0].set_xticklabels(np.linspace(0,300, 16), fontsize= 7)
    # ax[2,0].set_yticklabels(np.linspace(0,300, 16), fontsize= 7)


def PLOT_RERi (données_RER1, données_RER2) : 
    if nom_condition == "HH" :
        ax[0,0].plot(np.linspace(0,300, 16), données_RER1, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[0,1].plot(np.linspace(0,300, 16), données_RER2, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[0,0].plot(np.linspace(0,300, 16), données_RER1, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[0,1].plot(np.linspace(0,300, 16), données_RER2, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[0,0].plot(np.linspace(0,300, 16), données_RER1, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[0,1].plot(np.linspace(0,300, 16), données_RER2, color=LL_puits2, marker="+", label = "puits_2_LL")

    # ax[3,0].set_xticklabels(np.linspace(0,300, 16), fontsize= 7)
    # ax[3,0].set_yticklabels(np.linspace(0,300, 16), fontsize= 7)

def PLOT_Ui (données_U1, données_U2) : 
    if nom_condition == "HH" :
        ax[1,0].plot(np.linspace(0,300, 16), données_U1, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[1,1].plot(np.linspace(0,300, 16), données_U2, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[1,0].plot(np.linspace(0,300, 16), données_U1, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[1,1].plot(np.linspace(0,300, 16), données_U2, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[1,0].plot(np.linspace(0,300, 16), données_U1, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[1,1].plot(np.linspace(0,300, 16), données_U2, color=LL_puits2, marker="+", label = "puits_2_LL")

    # ax[4,0].set_xticklabels(np.linspace(0,300, 16), fontsize= 7)
    # ax[4,0].set_yticklabels(np.linspace(0,300, 16), fontsize= 7)

def PLOT_Vi (données_V1, données_V2) : 
    if nom_condition == "HH" :
        ax[2,0].plot(np.linspace(0,300, 16), données_V1, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[2,1].plot(np.linspace(0,300, 16), données_V2, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[2,0].plot(np.linspace(0,300, 16), données_V1, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[2,1].plot(np.linspace(0,300, 16), données_V2, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[2,0].plot(np.linspace(0,300, 16), données_V1, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[2,1].plot(np.linspace(0,300, 16), données_V2, color=LL_puits2, marker="+", label = "puits_2_LL")

    # ax[5,0].set_xticklabels(np.linspace(0,300, 16), fontsize= 7)
    # ax[5,0].set_yticklabels(np.linspace(0,300, 16), fontsize= 7)
    # ax[5,0].legend(title = "Traitement et puits", loc="upper right", fontsize= 7, bbox_to_anchor=(1.85, 6.4))

def RESISTANCES (longueur_entrenoeuds, rayon_entrenoeuds) :

    # Particule qui ne varie pas
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p)

    # calcul de la résistance
    eq = constante_resistance * (longueur_entrenoeuds)/(rayon_entrenoeuds**4)

    # 3 fois la même valeur R0 R1 R2
    E = 10
    eq = eq*E**(-6)

    eq = np.append(eq, (eq,eq))

    print("RESISAN", eq)

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
    
    eq = ((v_i * Ci_m)/(k_i + Ci_m))


    return eq

def VOLUME(Ci_m, v_i, k_i, VARIABLE_VOLUME) :

    eq = VARIABLE_VOLUME + (RER (Ci_m, v_i, k_i) * VARIABLE_VOLUME)

    return eq

def UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME) :
    
    eq = delta * VARIABLE_VOLUME * RER (Ci_m, v_i, k_i)

    return eq

def A_RESOUDRE (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds, v_i, k_i, VARIABLE_VOLUME) :

    eq = C0_m*FLUX_2puits (Ci_m, C0_m, longueur_entrenoeuds, rayon_entrenoeuds) - UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME)

    return eq
#### FONCTIONS

#### SOLVING
print("DEBUT LOOP CONDITION _________________________")
for condition, nom_condition in enumerate(nom_conditions) : 
    print("_________________________")

    ## INCREMENT C0_m
    if nom_condition == "LL"  :
        increment_C0_m = 0.0066 
        
    else :
        increment_C0_m = 0.05
    ## INCREMENT C0_m   

    v_i = np.array([BASE_v1[condition], BASE_v2[condition]])
    k_i = np.array([BASE_k1[condition], BASE_k2[condition]])

    longueur_entrenoeuds = BASE_longueur_entrenoeuds[condition]
    rayon_entrenoeuds = BASE_rayon_entrenoeuds[condition]

    C0_m_t0 = BASE_conc_ini_C0[condition]
    C1_m_t0 = BASE_conc_ini_C1[condition]
    C2_m_t0 = BASE_conc_ini_C2[condition]
    C1C2_ini = np.array([C1_m_t0, C2_m_t0])
    Ci_m = C1C2_ini

    VARIABLE_VOLUME = np.array([BASE_volume_ini_feuilles[condition], BASE_volume_ini_bourgeon[condition]])
    VARIABLE_C0_m = C0_m_t0

    # RESET SUR LES APPEND ENTRE CONDITIONS
    ARRAY_C0_m = np.array([C0_m_t0])
    ARRAY_C1_m = np.array([C1_m_t0])
    ARRAY_C2_m = np.array([C2_m_t0])
    ARRAY_F01 = np.array([FLUX_2puits (C1C2_ini, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[0]])
    ARRAY_F02 = np.array([FLUX_2puits (C1C2_ini, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[1]])
    ARRAY_RER1 = np.array([RER (C1C2_ini, v_i, k_i)[0]])
    ARRAY_RER2 = np.array([RER (C1C2_ini, v_i, k_i)[1]])
    ARRAY_U1 = np.array([UTILISATION (C1C2_ini, v_i, k_i, VARIABLE_VOLUME)[0]])
    ARRAY_U2 = np.array([UTILISATION (C1C2_ini, v_i, k_i,  VARIABLE_VOLUME)[1]])
    ARRAY_volume_1 = np.array([VARIABLE_VOLUME[0]])
    ARRAY_volume_2 = np.array([VARIABLE_VOLUME[1]])
    # RESET SUR LES APPEND ENTRE CONDITIONS
 

## BOUCLE RESOLUTION TEMPS ____________________________________________________________________________________________________________
    for t in range(0, 300, degre_jour) :

        # RESOLUTION
        Ci_m = scipy.fsolve(A_RESOUDRE, x0=(Ci_m), args=(VARIABLE_C0_m, longueur_entrenoeuds, rayon_entrenoeuds, v_i, k_i, VARIABLE_VOLUME), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
        print("Les solutions à l'équilibre C1_m et C2_m pour la condition", nom_conditions[condition], "sont :", Ci_m)

        ## ACTUALISATION C0_m
        VARIABLE_C0_m = VARIABLE_C0_m + VARIABLE_C0_m * increment_C0_m
        ## ACTUALISATION C0_m

        
        ## ACTUALISATION VOLUME ET DILUTION
        Mi_m = Ci_m * VARIABLE_VOLUME
        VARIABLE_VOLUME = VOLUME(Ci_m, v_i, k_i, VARIABLE_VOLUME)
        Ci_m = Mi_m / VARIABLE_VOLUME
        ## ACTUALISATION VOLUME ET DILUTION

        ## APPEND DES ARRAYS A TRACER
        ARRAY_C0_m = np.append(ARRAY_C0_m, VARIABLE_C0_m)
        ARRAY_C1_m = np.append(ARRAY_C1_m, Ci_m[0])
        ARRAY_C2_m = np.append(ARRAY_C2_m, Ci_m[1])

        ARRAY_F01 = np.append(ARRAY_F01, FLUX_2puits (Ci_m, VARIABLE_C0_m, longueur_entrenoeuds, rayon_entrenoeuds)[0])
        ARRAY_F02 = np.append(ARRAY_F02, FLUX_2puits (Ci_m, VARIABLE_C0_m, longueur_entrenoeuds, rayon_entrenoeuds)[1])

        ARRAY_RER1 = np.append(ARRAY_RER1, RER (Ci_m, v_i, k_i)[0])
        ARRAY_RER2 = np.append(ARRAY_RER2, RER (Ci_m, v_i, k_i)[1])

        ARRAY_U1 = np.append(ARRAY_U1, UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME)[0])
        ARRAY_U2 = np.append(ARRAY_U2, UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME)[1])

        ARRAY_volume_1 = np.append(ARRAY_volume_1, VARIABLE_VOLUME[0])
        ARRAY_volume_2 = np.append(ARRAY_volume_2, VARIABLE_VOLUME[1])

        if nom_condition == "HH" :
            ARRAY_C0_m_HH = ARRAY_C0_m
            ARRAY_C1_m_HH = ARRAY_C1_m
            ARRAY_C2_m_HH = ARRAY_C2_m

            ARRAY_F01_HH = ARRAY_F01
            ARRAY_F02_HH = ARRAY_F02

            ARRAY_RER1_HH = ARRAY_RER1
            ARRAY_RER2_HH = ARRAY_RER2

            ARRAY_U1_HH = ARRAY_U1
            ARRAY_U2_HH = ARRAY_U2

            ARRAY_volume_1_HH = ARRAY_volume_1
            ARRAY_volume_2_HH = ARRAY_volume_2
        else :
            if nom_condition == "LH" :
                ARRAY_C0_m_LH = ARRAY_C0_m
                ARRAY_C1_m_LH = ARRAY_C1_m
                ARRAY_C2_m_LH = ARRAY_C2_m

                ARRAY_F01_LH = ARRAY_F01
                ARRAY_F02_LH = ARRAY_F02

                ARRAY_RER1_LH = ARRAY_RER1
                ARRAY_RER2_LH = ARRAY_RER2

                ARRAY_U1_LH = ARRAY_U1
                ARRAY_U2_LH = ARRAY_U2

                ARRAY_volume_1_LH = ARRAY_volume_1
                ARRAY_volume_2_LH = ARRAY_volume_2
            else :
                ARRAY_C0_m_LL = ARRAY_C0_m
                ARRAY_C1_m_LL = ARRAY_C1_m
                ARRAY_C2_m_LL = ARRAY_C2_m

                ARRAY_F01_LL = ARRAY_F01
                ARRAY_F02_LL = ARRAY_F02

                ARRAY_RER1_LL = ARRAY_RER1
                ARRAY_RER2_LL = ARRAY_RER2

                ARRAY_U1_LL = ARRAY_U1
                ARRAY_U2_LL = ARRAY_U2

                ARRAY_volume_1_LL = ARRAY_volume_1
                ARRAY_volume_2_LL = ARRAY_volume_2
    print("_________________________")
print("FIN LOOP CONDITION _________________________")

print("C1M", ARRAY_C1_m_HH)
print("C2M", ARRAY_C2_m_HH)


fig, ax = plt.subplots(3,2)
nom_condition = "HH"
PLOT_C0m(ARRAY_C0_m_HH)
PLOT_FLUX(ARRAY_F01_HH, ARRAY_F02_HH)
PLOT_Cim(ARRAY_C1_m_HH, ARRAY_C2_m_HH)

nom_condition = "LH"
PLOT_C0m(ARRAY_C0_m_LH)
PLOT_FLUX(ARRAY_F01_LH, ARRAY_F02_LH)
PLOT_Cim(ARRAY_C1_m_LH, ARRAY_C2_m_LH)

nom_condition = "LL"
PLOT_C0m(ARRAY_C0_m_LL)
PLOT_FLUX(ARRAY_F01_LL, ARRAY_F02_LL)
PLOT_Cim(ARRAY_C1_m_LL, ARRAY_C2_m_LL)

ax[0,0].set_title("C0m simulé", fontsize=7)
ax[0,0].set_ylabel("umolC/°Cj", color="black", fontsize=7)
ax[0,0].legend(fontsize= 7)

ax[1,0].set_title("FLUX01 simulé", fontsize=7)
ax[1,0].set_ylabel("umolC/°Cj", color="black", fontsize=7)
ax[1,1].set_title("FLUX02 simulé", fontsize=7)
ax[1,1].set_ylabel("umolC/°Cj", color="black", fontsize=7)

ax[2,0].set_title("C1_m simulé", fontsize=7)
ax[2,0].set_ylabel("umol/m3", color="black", fontsize=7)
ax[2,1].set_title("C2_m simulé", fontsize=7)
ax[2,1].set_ylabel("umol/m3", color="black", fontsize=7)

ax[1,0].legend(fontsize= 7)
ax[1,1].legend(fontsize= 7)

ax[0, 1].axis('off')



fig, ax = plt.subplots(3,2)
nom_condition = "HH"
PLOT_RERi(ARRAY_RER1_HH, ARRAY_RER2_HH)
PLOT_Ui(ARRAY_U1_HH, ARRAY_U2_HH)
PLOT_Vi(ARRAY_volume_1_HH, ARRAY_volume_2_HH)

nom_condition = "LH"
PLOT_RERi(ARRAY_RER1_LH, ARRAY_RER2_LH)
PLOT_Ui(ARRAY_U1_LH, ARRAY_U2_LH)
PLOT_Vi(ARRAY_volume_1_LH, ARRAY_volume_2_LH)
nom_condition = "LL"
PLOT_RERi(ARRAY_RER1_LL, ARRAY_RER2_LL)
PLOT_Ui(ARRAY_U1_LL, ARRAY_U2_LL)
PLOT_Vi(ARRAY_volume_1_LL, ARRAY_volume_2_LL)

ax[0,0].set_title("RER1 simulé", fontsize=7)
ax[0,1].set_title("RER2 simulé", fontsize=7)
ax[0,0].set_ylabel("/°Cj", color="black", fontsize=7)
ax[0,1].set_ylabel("/°Cj", color="black", fontsize=7)

ax[1,0].set_title("U1 simulé", fontsize=7)
ax[1,1].set_title("U2 simulé", fontsize=7)
ax[1,0].set_ylabel("umolC/°Cj", color="black", fontsize=7)
ax[1,1].set_ylabel("umolC/°Cj", color="black", fontsize=7)

ax[2,0].set_title("V1 simulé", fontsize=7)
ax[2,1].set_title("V2 simulé", fontsize=7)
ax[2,0].set_ylabel("m3", color="black", fontsize=7)
ax[2,1].set_ylabel("m3", color="black", fontsize=7)


fig, ax = plt.subplots(3,1)
ax[0].plot(range(0, 300, 10), EXP_RER1_HH, color=HH_puits1, marker="+")
ax[0].plot(range(0, 300, 10), EXP_RER1_LH, color=LH_puits1, marker="+")
ax[0].plot(range(0, 300, 10), EXP_RER1_LL, color=LL_puits1, marker="+")

ax[1].plot(range(0, 300, 10), EXP_U1_HH, color=HH_puits1, marker="+")
ax[1].plot(range(0, 300, 10), EXP_U1_LH, color=LH_puits1, marker="+")
ax[1].plot(range(0, 300, 10), EXP_U1_LL, color=LL_puits1, marker="+")

ax[2].plot(range(0, 300, 10), EXP_V1_HH, color=HH_puits1, marker="+")
ax[2].plot(range(0, 300, 10), EXP_V1_LH, color=LH_puits1, marker="+")
ax[2].plot(range(0, 300, 10), EXP_V1_LL, color=LL_puits1, marker="+")

ax[0].set_title("RER1 observé", fontsize=7)
ax[0].set_ylabel("/°Cj", color="black", fontsize=7)

ax[1].set_title("U1 observé", fontsize=7)
ax[1].set_ylabel("umolC/°Cj", color="black", fontsize=7)

ax[2].set_title("V1 observé", fontsize=7)
ax[2].set_ylabel("m3", color="black", fontsize=7)
#### PLT SHOW

plt.show()


    