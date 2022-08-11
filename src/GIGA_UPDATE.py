#### MODULES
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES

## SETUP PLOTS
fig,ax = plt.subplots(6,2)
## SETUP PLOTS

#### DONNEES EXPERIMETALES
EXP_C1m_HH = np.array([2.773E+08, np.nan, 1.232E+06], dtype=float)
EXP_C1m_LH = np.array([2.773E+08, 1.011E+06, 1.080E+06], dtype=float)
EXP_C1m_LL = np.array([2.773E+08, 3.954E+05, 5.763E+05], dtype=float)

EXP_C2m_HH = np.array([5.358E+03, 2.548E+00, np.nan], dtype=float)
EXP_C2m_LH = np.array([5.358E+03, 5.247E+00, np.nan], dtype=float)
EXP_C2m_LL = np.array([5.358E+03, 8.631E-01, np.nan], dtype=float)

EXP_RER1_HH = np.array([6.644E-03, 6.281E-03, 5.923E-03, 5.570E-03, 5.223E-03, 4.882E-03, 4.548E-03, 4.221E-03, 3.903E-03, 3.594E-03, 3.296E-03, 3.009E-03, 2.735E-03, 2.474E-03, 2.229E-03, 1.998E-03, 1.784E-03, 1.586E-03, 1.404E-03, 1.238E-03, 1.087E-03, 9.519E-04, 8.307E-04, 7.230E-04, 6.275E-04, 5.434E-04, 4.696E-04, 4.051E-04, 3.488E-04, 3.000E-04], dtype=float)
EXP_RER1_LH = np.array([3.803E-03, 3.576E-03, 3.353E-03, 3.133E-03, 2.916E-03, 2.705E-03, 2.498E-03, 2.297E-03, 2.104E-03, 1.917E-03, 1.739E-03, 1.571E-03, 1.412E-03, 1.263E-03, 1.125E-03, 9.983E-04, 8.820E-04, 7.763E-04, 6.808E-04, 5.951E-04, 5.187E-04, 4.508E-04, 3.908E-04, 3.381E-04, 2.919E-04, 2.516E-04, 2.165E-04, 1.861E-04, 1.597E-04, 1.370E-04], dtype=float)
EXP_RER1_LL = np.array([4.056E-03, 3.825E-03, 3.596E-03, 3.370E-03, 3.148E-03, 2.929E-03, 2.715E-03, 2.506E-03, 2.303E-03, 2.107E-03, 1.919E-03, 1.739E-03, 1.569E-03, 1.409E-03, 1.259E-03, 1.120E-03, 9.927E-04, 8.761E-04, 7.702E-04, 6.748E-04, 5.893E-04, 5.131E-04, 4.455E-04, 3.859E-04, 3.336E-04, 2.879E-04, 2.480E-04, 2.133E-04, 1.832E-04, 1.572E-04], dtype=float)

EXP_U1_HH = np.array([2.556E+02, 2.573E+02, 2.574E+02, 2.559E+02, 2.526E+02, 2.477E+02, 2.413E+02, 2.333E+02, 2.240E+02, 2.136E+02, 2.021E+02, 1.899E+02, 1.772E+02, 1.641E+02, 1.509E+02, 1.379E+02, 1.252E+02, 1.129E+02, 1.013E+02, 9.036E+01, 8.019E+01, 7.083E+01, 6.230E+01, 5.458E+01, 4.766E+01, 4.148E+01, 3.600E+01, 3.117E+01, 2.693E+01, 2.322E+00], dtype=float)
EXP_U1_LH = np.array([1.409E+02, 1.442E+02, 1.468E+02, 1.488E+02, 1.502E+02, 1.508E+02, 1.508E+02, 1.502E+02, 1.490E+02, 1.472E+02, 1.450E+02, 1.423E+02, 1.394E+02, 1.362E+02, 1.328E+02, 1.294E+02, 1.260E+02, 1.227E+02, 1.194E+02, 1.163E+02, 1.134E+02, 1.106E+02, 1.081E+02, 1.038E+02, 1.805E+01, 1.561E+01, 1.347E+01, 1.161E+01, 9.985E+00, 8.578E+00], dtype=float)
EXP_U1_LL = np.array([1.213E+02, 1.203E+02, 1.186E+02, 1.163E+02, 1.134E+02, 1.099E+02, 1.058E+02, 1.012E+02, 9.611E+01, 9.068E+01, 8.496E+01, 7.905E+01, 7.304E+01, 6.702E+01, 6.110E+01, 5.533E+01, 4.981E+01, 4.458E+01, 3.968E+01, 3.515E+01, 3.099E+01, 2.721E+01, 2.380E+01, 2.075E+01, 1.804E+01, 1.564E+01, 1.353E+01, 1.167E+01, 1.006E+01, 8.654E+00], dtype=float)

EXP_V1_HH = np.array([2.12668E-07, 2.40575E-07, 2.7086E-07, 3.03418E-07, 3.38072E-07, 3.74571E-07, 4.12594E-07, 4.51764E-07, 4.91654E-07, 5.31814E-07, 5.71782E-07, 6.11114E-07, 6.49394E-07, 6.86257E-07, 7.21398E-07, 7.54581E-07, 7.85636E-07, 8.14466E-07, 8.41032E-07, 8.65352E-07, 8.87488E-07, 9.07538E-07, 9.25625E-07, 9.4189E-07, 9.56481E-07, 9.69551E-07, 9.8125E-07, 9.91723E-07, 1.00111E-06, 1.00953E-06], dtype=float)
EXP_V1_LH = np.array([3.15765E-07, 3.49798E-07, 3.85339E-07, 4.22075E-07, 4.5965E-07, 4.97681E-07, 5.35775E-07, 5.7354E-07, 6.10601E-07, 6.46615E-07, 6.8128E-07, 7.1434E-07, 7.45596E-07, 7.74901E-07, 8.02162E-07, 8.27339E-07, 8.50433E-07, 8.71488E-07, 8.90577E-07, 9.07799E-07, 9.23271E-07, 9.37122E-07, 9.49485E-07, 9.60496E-07, 9.70289E-07, 9.78992E-07, 9.86725E-07, 9.936E-07, 9.9972E-07, 1.00518E-06], dtype=float)
EXP_V1_LL = np.array([4.3972E-07, 4.87299E-07, 5.37054E-07, 5.88559E-07, 6.41327E-07, 6.94832E-07, 7.48526E-07, 8.01859E-07, 8.54304E-07, 9.05368E-07, 9.54611E-07, 1.00166E-06, 1.04621E-06, 1.08803E-06, 1.12697E-06, 1.16295E-06, 1.19596E-06, 1.22604E-06, 1.25328E-06, 1.27782E-06, 1.29982E-06, 1.31945E-06, 1.33691E-06, 1.35238E-06, 1.36606E-06, 1.37814E-06, 1.38879E-06, 1.39818E-06, 1.40645E-06, 1.41375E-06], dtype=float)

DATES = np.array([0, 130, 240])

HH_puits1 = "#39f559"
LH_puits1 = "#218E34"
LL_puits1 = "#0C3514"

HH_puits2 = "#39e7f1"
LH_puits2 = "#249096"
LL_puits2 = "#114346"

HH_source0 = "#F7ed35"
LH_source0 = "#A49E23"
LL_source0 = "#4C4910"

# Ci_m
ax[2,1].scatter(DATES, EXP_C1m_HH, color=HH_puits1, marker="+")
ax[2,1].scatter(DATES, EXP_C1m_LH, color=LH_puits1, marker="+")
ax[2,1].scatter(DATES, EXP_C1m_LL, color=LL_puits1, marker="+")
ax[2,1].scatter(DATES, EXP_C2m_HH, color=HH_puits2, marker="+")
ax[2,1].scatter(DATES, EXP_C2m_LH, color=LH_puits2, marker="+")
ax[2,1].scatter(DATES, EXP_C2m_LL, color=LL_puits2, marker="+")

ax[3,1].plot(range(0, 300, 10), EXP_RER1_HH, color=HH_puits1, marker="+")
ax[3,1].plot(range(0, 300, 10), EXP_RER1_LH, color=LH_puits1, marker="+")
ax[3,1].plot(range(0, 300, 10), EXP_RER1_LL, color=LL_puits1, marker="+")

ax[4,1].plot(range(0, 300, 10), EXP_U1_HH, color=HH_puits1, marker="+")
ax[4,1].plot(range(0, 300, 10), EXP_U1_LH, color=LH_puits1, marker="+")
ax[4,1].plot(range(0, 300, 10), EXP_U1_LL, color=LL_puits1, marker="+")

ax[5,1].plot(range(0, 300, 10), EXP_V1_HH, color=HH_puits1, marker="+")
ax[5,1].plot(range(0, 300, 10), EXP_V1_LH, color=LH_puits1, marker="+")
ax[5,1].plot(range(0, 300, 10), EXP_V1_LL, color=LL_puits1, marker="+")
#### DONNEES EXPERIMETALES

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

# BASE_volume_ini_bourgeon = np.array([8.849E-10, 8.849E-10, 8.849E-10], dtype=float)
# BASE_volume_ini_feuilles = np.array([4.422E-07, 4.422E-07, 4.422E-07], dtype=float)

# VALEURS ACTUALISEES CAR JE SAIS PAS FAIRE DE CALCULS
BASE_volume_ini_bourgeon = np.array([8.849E-10, 8.849E-10, 8.849E-10], dtype=float)
BASE_volume_ini_feuilles = np.array([2.12668E-07, 3.15765E-07, 4.3972E-07], dtype=float)


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
def PLOT_C0m (données_C0m) : 
    if nom_condition == "HH" :
        ax[0,0].plot(np.linspace(0,300, 16), données_C0m, color=HH_source0, marker="+", label = "source_0_HH")
    else :
        if nom_condition == "LH" :
            ax[0,0].plot(np.linspace(0,300, 16), données_C0m, color=LH_source0, marker="+", label = "source_0_LH")
        else :
            ax[0,0].plot(np.linspace(0,300, 16), données_C0m, color=LL_source0, marker="+", label = "source_0_LL")

    ax[0,0].legend()

def PLOT_FLUX (données_F01, données_F02) : 
    if nom_condition == "HH" :
        ax[1,0].plot(np.linspace(0,300, 16), données_F01, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[1,0].plot(np.linspace(0,300, 16), données_F02, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[1,0].plot(np.linspace(0,300, 16), données_F01, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[1,0].plot(np.linspace(0,300, 16), données_F02, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[1,0].plot(np.linspace(0,300, 16), données_F01, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[1,0].plot(np.linspace(0,300, 16), données_F02, color=LL_puits2, marker="+", label = "puits_2_LL")

def PLOT_Cim (données_C1m, données_C2m) : 
    if nom_condition == "HH" :
        ax[2,0].plot(np.linspace(0,300, 16), données_C1m, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[2,0].plot(np.linspace(0,300, 16), données_C2m, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[2,0].plot(np.linspace(0,300, 16), données_C1m, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[2,0].plot(np.linspace(0,300, 16), données_C2m, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[2,0].plot(np.linspace(0,300, 16), données_C1m, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[2,0].plot(np.linspace(0,300, 16), données_C2m, color=LL_puits2, marker="+", label = "puits_2_LL")

def PLOT_RERi (données_RER1, données_RER2) : 
    if nom_condition == "HH" :
        ax[3,0].plot(np.linspace(0,300, 16), données_RER1, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[3,0].plot(np.linspace(0,300, 16), données_RER2, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[3,0].plot(np.linspace(0,300, 16), données_RER1, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[3,0].plot(np.linspace(0,300, 16), données_RER2, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[3,0].plot(np.linspace(0,300, 16), données_RER1, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[3,0].plot(np.linspace(0,300, 16), données_RER2, color=LL_puits2, marker="+", label = "puits_2_LL")

def PLOT_Ui (données_U1, données_U2) : 
    if nom_condition == "HH" :
        ax[4,0].plot(np.linspace(0,300, 16), données_U1, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[4,0].plot(np.linspace(0,300, 16), données_U2, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[4,0].plot(np.linspace(0,300, 16), données_U1, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[4,0].plot(np.linspace(0,300, 16), données_U2, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[4,0].plot(np.linspace(0,300, 16), données_U1, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[4,0].plot(np.linspace(0,300, 16), données_U2, color=LL_puits2, marker="+", label = "puits_2_LL")

def PLOT_Vi (données_V1, données_V2) : 
    if nom_condition == "HH" :
        ax[5,0].plot(np.linspace(0,300, 16), données_V1, color=HH_puits1, marker="+", label = "puits_1_HH")
        ax[5,0].plot(np.linspace(0,300, 16), données_V2, color=HH_puits2, marker="+", label = "puits_2_HH")
    else :
        if nom_condition == "LH" :
            ax[5,0].plot(np.linspace(0,300, 16), données_V1, color=LH_puits1, marker="+", label = "puits_1_LH")
            ax[5,0].plot(np.linspace(0,300, 16), données_V2, color=LH_puits2, marker="+", label = "puits_2_LH")
        else :
            ax[5,0].plot(np.linspace(0,300, 16), données_V1, color=LL_puits1, marker="+", label = "puits_1_LL")
            ax[5,0].plot(np.linspace(0,300, 16), données_V2, color=LL_puits2, marker="+", label = "puits_2_LL")

    ax[5,0].legend(title = "Traitement et puits", loc="upper right")


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

#### SOLVING
print("DEBUT LOOP CONDITION _________________________")
for condition, nom_condition in enumerate(nom_conditions) : 
    print("_________________________")

    v_i = np.array([BASE_v1[condition], BASE_v2[condition]])
    k_i = np.array([BASE_k1[condition], BASE_k2[condition]])

    longueur_entrenoeuds = BASE_longueur_entrenoeuds[condition]
    rayon_entrenoeuds = BASE_rayon_entrenoeuds[condition]

    C0_m_t0 = BASE_conc_ini_C0[condition]
    C1_m_t0 = BASE_conc_ini_C1[condition]
    C2_m_t0 = BASE_conc_ini_C2[condition]
    C1C2_ini = np.array([C1_m_t0, C2_m_t0])
    Ci_m = C1C2_ini

    V_t0 = np.array([BASE_volume_ini_feuilles[condition], BASE_volume_ini_bourgeon[condition]])
    V_t = V_t0

    VARIABLE_VOLUME = V_t0
    VARIABLE_C0_m = C0_m_t0

    # RESET SUR LES APPEND ENTRE CONDITIONS
    ARRAY_C0_m = np.array([C0_m_t0])
    ARRAY_C1_m = np.array([C1_m_t0])
    ARRAY_C2_m = np.array([C2_m_t0])
    ARRAY_F01 = np.array([FLUX_2puits (C1C2_ini, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[0]])
    ARRAY_F02 = np.array([FLUX_2puits (C1C2_ini, C0_m_t0, longueur_entrenoeuds, rayon_entrenoeuds)[1]])
    ARRAY_RER1 = np.array([RER (C1C2_ini, v_i, k_i)[0]])
    ARRAY_RER2 = np.array([RER (C1C2_ini, v_i, k_i)[1]])
    ARRAY_U1 = np.array([UTILISATION (C1C2_ini, v_i, k_i, V_t)[0]])
    ARRAY_U2 = np.array([UTILISATION (C1C2_ini, v_i, k_i, V_t)[1]])
    ARRAY_volume_1 = np.array([VARIABLE_VOLUME[0]])
    ARRAY_volume_2 = np.array([VARIABLE_VOLUME[1]])
    # RESET SUR LES APPEND ENTRE CONDITIONS
    
    ## INCREMENT C0_m
    if nom_condition == "LL"  :
        increment_C0_m = 0.0066 
        
    else :
        increment_C0_m = 0.05
    ## INCREMENT C0_m    

## BOUCLE RESOLUTION TEMPS ____________________________________________________________________________________________________________
    for t in range(0, 300, degre_jour) :

        # RESOLUTION
        Ci_m = scipy.fsolve(A_RESOUDRE, x0=(Ci_m), args=(VARIABLE_C0_m, longueur_entrenoeuds, rayon_entrenoeuds, v_i, k_i, VARIABLE_VOLUME), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
        print("Les solutions à l'équilibre C1_m et C2_m pour la condition", nom_conditions[condition], "sont :", Ci_m)

        ## ACTUALISATION C0_m
        VARIABLE_C0_m = VARIABLE_C0_m + VARIABLE_C0_m * increment_C0_m
        ## ACTUALISATION C0_m

        ## ACTUALISATION VOLUME
        VARIABLE_VOLUME = VOLUME(Ci_m, v_i, k_i, VARIABLE_VOLUME)
        ## ACTUALISATION VOLUME

        ## APPEND DES ARRAYS A TRACER
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

    PLOT_C0m(ARRAY_C0_m)
    PLOT_FLUX(ARRAY_F01, ARRAY_F02)
    PLOT_Cim(ARRAY_C1_m, ARRAY_C2_m)
    PLOT_RERi(ARRAY_RER1, ARRAY_RER2)
    PLOT_Ui(ARRAY_U1, ARRAY_U2)
    PLOT_Vi(ARRAY_volume_1, ARRAY_volume_2)

    print(ARRAY_F01)

    print("_________________________")
print("FIN LOOP CONDITION _________________________")

#### PLT SHOW
ax[0,0].set_title("C0m")
ax[0,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)

ax[1,0].set_title("FLUX01 et FLUX02")
ax[1,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)

ax[2,0].set_title("C1_m et C2_m var")
ax[2,0].set_ylabel("umol/m3", color="black", fontsize=14)

ax[3,0].set_title("RER1 et RER2 var")
ax[3,0].set_ylabel("/°Cj", color="black", fontsize=14)

ax[4,0].set_title("U1 et U2 var")
ax[4,0].set_ylabel("umolC/°Cj", color="black", fontsize=14)

ax[5,0].set_title("V1 et V2")
ax[5,0].set_ylabel("m3", color="black", fontsize=14)

plt.show()
#### PLT SHOW

    