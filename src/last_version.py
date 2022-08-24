#### MODULES ___________________________________________________________________________________________________________________________________________________________
import scipy.optimize as scipy 
import matplotlib.pyplot as plt
import numpy as np
import math
#### MODULES ___________________________________________________________________________________________________________________________________________________________


#### DONNEES EXPERIMETALES _____________________________________________________________________________________________________________________________________________
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
#### DONNEES EXPERIMETALES _____________________________________________________________________________________________________________________________________________


#### CONSTANTES ________________________________________________________________________________________________________________________________________________________
temp_20 = 293
gaz_p = 8.314
viscosity = 1E6
nom_conditions = np.array(["HH", "LH", "LL"], dtype=str)
delta = 1.694E+10
degre_jour = 20
correction_unite_resistance = 8.46E10
#### CONSTANTES ________________________________________________________________________________________________________________________________________________________


#### VARIABLES INTERCONDITIONS _________________________________________________________________________________________________________________________________________
BASE_longueur_entrenoeuds = np.array([4.339E-2, 4.339E-2, 4.339E-2], dtype=float)
# BASE_longueur_entrenoeuds = np.array([20E-2, 20E-2, 20E-2], dtype=float)
BASE_rayon_entrenoeuds = np.array([35.1E-6, 35.1E-6, 35.1E-6], dtype=float)

BASE_volume_ini_bourgeon = np.array([8.849E-10, 8.849E-10, 8.849E-10], dtype=float)
BASE_volume_ini_feuilles = np.array([1.792E-06, 1.792E-06, 1.792E-06], dtype=float)

BASE_v1 = np.array([1.779E-02, 1.012E-02, 1.012E-02], dtype=float)
BASE_k1 = np.array([1.183E+12, 1.183E+12, 1.183E+12], dtype=float)
BASE_v2 = np.array([3.150E-02, 3.150E-02, 3.150E-02], dtype=float)
BASE_k2 = np.array([2.099E+12, 2.099E+12, 2.099E+12], dtype=float)

BASE_conc_ini_C0 = np.array([1.525E+09, 1.525E+09, 1.525E+09], dtype=float)
# BASE_conc_ini_C0 = np.array([1E+08, 1E+08, 1E+08], dtype=float)
BASE_conc_ini_C1 = np.array([2.773E+08, 2.773E+08, 2.773E+08], dtype=float)
BASE_conc_ini_C2 = np.array([5.358E+03, 5.358E+03, 5.358E+03], dtype=float)
#### VARIABLES INTERCONDITIONS _________________________________________________________________________________________________________________________________________


# ### VARIABLES ENTRE CONDITIONS
# BASE_longueur_entrenoeuds = np.array([4.339E-2, 4.339E-2, 4.339E-2], dtype=float)
# BASE_rayon_entrenoeuds = np.array([35.1E-6, 35.1E-6, 35.1E-6], dtype=float)

# BASE_volume_ini_bourgeon = np.array([8.849E-10, 8.849E-10, 8.849E-10], dtype=float)
# BASE_volume_ini_feuilles = np.array([2.12668E-07, 2.12668E-07, 2.12668E-07], dtype=float)

# # Ci_m_t0 identiques
# # BASE_v1 = np.array([1.779E-02, 1.329E-02, 1.329E-02], dtype=float)
# # BASE_k1 = np.array([1.387E+08, 1.387E+08, 1.387E+08], dtype=float)
# # BASE_v2 = np.array([2.100E-02, 2.100E-02, 2.100E-02], dtype=float)
# # BASE_k2 = np.array([1.387E+08, 1.387E+08, 1.387E+08], dtype=float)

# BASE_v1 = np.array([1.779E-02, 1.329E-02, 1.329E-02], dtype=float)
# BASE_k1 = np.array([20*1.387E+08, 20*1.387E+08, 20*1.387E+08], dtype=float)
# BASE_v2 = np.array([2.100E-02, 2.100E-02, 2.100E-02], dtype=float)
# BASE_k2 = np.array([20*1.387E+08, 20*1.387E+08, 20*1.387E+08], dtype=float)



# # BASE_conc_ini_C0 = np.array([1.525E+09, 1.525E+09, 1.525E+09], dtype=float)
# # BASE_conc_ini_C1 = np.array([1.387E+08, 1.387E+08, 1.387E+08], dtype=float)
# # BASE_conc_ini_C2 = np.array([5.358E+03, 5.358E+03, 5.358E+03], dtype=float)

# BASE_conc_ini_C0 = np.array([150*1.525E+09, 150*1.525E+09, 150*1.525E+09], dtype=float)
# BASE_conc_ini_C1 = np.array([20*1.387E+08, 20*1.387E+08, 20*1.387E+08], dtype=float)
# BASE_conc_ini_C2 = np.array([20*1.387E+08, 20*1.387E+08, 20*1.387E+08], dtype=float)
# #### VARIABLES ENTRE CONDITIONS


#### FONCTIONS _________________________________________________________________________________________________________________________________________________________
def RESISTANCES_T0 (longueur_entrenoeuds, rayon_entrenoeuds) :
    # eq = (((8*longueur_entrenoeuds*viscosity) /(rayon_entrenoeuds**4*math.pi)) / (gaz_p*temp_20))*correction_unite_resistance
    eq = np.array([5E18, 5E18, 5E18], dtype=float)
    eq[2] = eq[2] * 100

    return eq


def RESISTANCES_VARIATION (Ci_m, v_i, k_i, VARIABLE_VOLUME, longueur_entrenoeuds, rayon_entrenoeuds, resistances_temps_t) :
        # SI T=0
    if t == 0 :
        eq = RESISTANCES_T0 (longueur_entrenoeuds, rayon_entrenoeuds)

    else : 
        eq = resistances_temps_t
    
    # si R2 baisse trop
    if eq[2] <= eq[0] :
        eq[2] = eq[0]

    # Baisse classique R2
    else :
        eq[2] = eq[2] - eq[2] * UTILISATION(Ci_m, v_i, k_i, VARIABLE_VOLUME)[1]
        # eq[2] = 1
    
    return eq


def SIGNAL_PUITS2 () :
    if nom_condition == "LH" :
        eq = 3
    else :
        if nom_condition == "HH" :
            eq = 1.5
        else : 
            eq = 0.25
    return eq


def RER (Ci_m, v_i, k_i) :
    signal_puits2 = SIGNAL_PUITS2 ()

    RER_puits1 = ((v_i[0] * Ci_m[0])/(k_i[0] + Ci_m[0]))

    RER_puits2 = ((v_i[1] * Ci_m[1])/(k_i[1] + Ci_m[1]))*signal_puits2

    eq = np.append(RER_puits1, RER_puits2)   

    return eq


def UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME) :
    eq = delta * VARIABLE_VOLUME * RER (Ci_m, v_i, k_i)

    return eq


def CROISSANCE_VOLUME (Ci_m, v_i, k_i, VARIABLE_VOLUME) :
    eq = VARIABLE_VOLUME + (RER (Ci_m, v_i, k_i) * VARIABLE_VOLUME)

    return eq


def C0M_VARIATION (C0_m) :
    ## INCREMENT C0_m
    if nom_condition == "LL" :
        increment_C0_m = 0.012
    else :
        increment_C0_m = 0.048
    ## INCREMENT C0_m   

    eq = C0_m + increment_C0_m*BASE_conc_ini_C0[condition]

    return eq


def FLUX_2puits (Ci_m, C0_m) :
    R_ = resistances_temps_t
    R_0 = R_[0]
    R_1 = R_[1]
    R_2 = R_[2]
    denominateur = R_0*(R_1+R_2)+R_1*R_2

    F01 = C0_m*((R_2*(C0_m-Ci_m[0]) + R_0*(Ci_m[1]-Ci_m[0])) / (denominateur))
    F02 = C0_m*((R_1*(C0_m-Ci_m[1]) + R_0*(Ci_m[0]-Ci_m[1])) / (denominateur))

    eq = np.append(F01,F02)

    return eq


def A_RESOUDRE (Ci_m, C0_m, v_i, k_i, VARIABLE_VOLUME) :

    eq = FLUX_2puits (Ci_m, C0_m) - UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME)

    return eq
#### FONCTIONS _________________________________________________________________________________________________________________________________________________________


#### SOLVING ___________________________________________________________________________________________________________________________________________________________
print("DEBUT DES CALCULS _________________________")
for condition, nom_condition in enumerate(nom_conditions) : 
    # SETUP DU T0
    # deux values
    v_i = np.array([BASE_v1[condition], BASE_v2[condition]])
    k_i = np.array([BASE_k1[condition], BASE_k2[condition]])
    VARIABLE_VOLUME = np.array([BASE_volume_ini_feuilles[condition], BASE_volume_ini_bourgeon[condition]])

    # une value
    longueur_entrenoeuds = BASE_longueur_entrenoeuds
    rayon_entrenoeuds = BASE_rayon_entrenoeuds
    VARIABLE_C0_m = BASE_conc_ini_C0[condition]
    VARIABLE_C1_m = BASE_conc_ini_C1[condition]
    VARIABLE_C2_m = BASE_conc_ini_C2[condition]

    # deux values
    Ci_m = np.array([VARIABLE_C1_m, VARIABLE_C2_m])
    resistances_temps_t = RESISTANCES_T0 (longueur_entrenoeuds, rayon_entrenoeuds)
    # SETUP DU T0

    ## RESET ARRAYS DATA
    ARRAY_C0_m = np.array([VARIABLE_C0_m])
    ARRAY_C1_m = np.array([VARIABLE_C1_m])
    ARRAY_C2_m = np.array([VARIABLE_C2_m])
    ARRAY_R0 = np.array([resistances_temps_t[0]])
    ARRAY_R1 = np.array([resistances_temps_t[1]])
    ARRAY_R2 = np.array([resistances_temps_t[2]])
    ARRAY_V1 = np.array(VARIABLE_VOLUME[0])
    ARRAY_V2 = np.array(VARIABLE_VOLUME[1])
    ARRAY_F01 = np.array(FLUX_2puits (Ci_m, VARIABLE_C0_m)[0])
    ARRAY_F02 = np.array(FLUX_2puits (Ci_m, VARIABLE_C0_m)[1])
    ARRAY_RER1 = np.array(RER (Ci_m, v_i, k_i)[0])
    ARRAY_RER2 = np.array(RER (Ci_m, v_i, k_i)[1])
    ARRAY_U1 = np.array(UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME)[0])
    ARRAY_U2 = np.array(UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME)[1])
    ## RESET ARRAYS DATA

    #### BOUCLE TEMPS ______________________________________________________________________________________________________________________________________________________
    for t in range(0, 300, degre_jour) :

            # print (RESISTANCES_VARIATION (Ci_m, v_i, k_i, VARIABLE_VOLUME, longueur_entrenoeuds, rayon_entrenoeuds, resistances_temps_t))

            # RESOLUTION CI_m
            Ci_m = scipy.fsolve(A_RESOUDRE, x0=(Ci_m), args=(VARIABLE_C0_m, v_i, k_i, VARIABLE_VOLUME), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None) 
            # print("Les solutions à l'équilibre C1_m et C2_m pour la condition", nom_conditions[condition], "sont :", Ci_m)
            # RESOLUTION CI_m

            # ANTI NEGATIF SUR C1_m et C2_m
            if Ci_m[0] <= 0 :
                Ci_m[0] = 0
            if Ci_m[1] <= 0 :
                Ci_m[1] = 0
            # ANTI NEGATIF SUR C1_m et C2_m

            # ACTUALISATION C0_m
            VARIABLE_C0_m = C0M_VARIATION (VARIABLE_C0_m)
            # ACTUALISATION C0_m
        
            # ACTUALISATION VOLUME ET DILUTION Ci_m
            Mi_m = Ci_m * VARIABLE_VOLUME
            VARIABLE_VOLUME = CROISSANCE_VOLUME(Ci_m, v_i, k_i, VARIABLE_VOLUME)
            Ci_m = Mi_m / VARIABLE_VOLUME
            # ACTUALISATION VOLUME ET DILUTION Ci_m

            # ACTUALISATION RESISTANCES
            resistances_temps_t = RESISTANCES_VARIATION (Ci_m, v_i, k_i, VARIABLE_VOLUME, longueur_entrenoeuds, rayon_entrenoeuds, resistances_temps_t)
            # ACTUALISATION RESISTANCES

            # ARRAYS POUR EXTRACTION DATA 
            ARRAY_C0_m = np.append(ARRAY_C0_m, VARIABLE_C0_m)
            ARRAY_C1_m = np.append(ARRAY_C1_m, Ci_m[0])
            ARRAY_C2_m = np.append(ARRAY_C2_m, Ci_m[1])
            ARRAY_R0 = np.append(ARRAY_R0, resistances_temps_t[0])
            ARRAY_R1 = np.append(ARRAY_R1, resistances_temps_t[1])
            ARRAY_R2 = np.append(ARRAY_R2, resistances_temps_t[2])
            ARRAY_V1 = np.append(ARRAY_V1, VARIABLE_VOLUME[0])
            ARRAY_V2 = np.append(ARRAY_V2, VARIABLE_VOLUME[1])
            ARRAY_F01 = np.append(ARRAY_F01, FLUX_2puits (Ci_m, VARIABLE_C0_m)[0])
            ARRAY_F02 = np.append(ARRAY_F02, FLUX_2puits (Ci_m, VARIABLE_C0_m)[1])
            ARRAY_RER1 = np.append(ARRAY_RER1, RER (Ci_m, v_i, k_i)[0])
            ARRAY_RER2 = np.append(ARRAY_RER2, RER (Ci_m, v_i, k_i)[1])
            ARRAY_U1 = np.append(ARRAY_U1, UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME)[0])
            ARRAY_U2 = np.append(ARRAY_U2, UTILISATION (Ci_m, v_i, k_i, VARIABLE_VOLUME)[1])
            # ARRAYS POUR EXTRACTION DATA 

    #### BOUCLE TEMPS ______________________________________________________________________________________________________________________________________________________

    # print("R2", ARRAY_R2)
    # print("R2", ARRAY_R1)
    # print("VOLUME 1", ARRAY_V1)
    # print("VOLUME 2", ARRAY_V2)

    # SETUP POUR EXTRACTION DATA 
    if nom_condition == "HH" :
        HH_ARRAY_C0_m = ARRAY_C0_m
        HH_ARRAY_C1_m = ARRAY_C1_m
        HH_ARRAY_C2_m = ARRAY_C2_m
        HH_ARRAY_R0 = ARRAY_R0
        HH_ARRAY_R1 = ARRAY_R1
        HH_ARRAY_R2 = ARRAY_R2
        HH_ARRAY_V1 = ARRAY_V1
        HH_ARRAY_V2 = ARRAY_V2
        HH_ARRAY_F01 = ARRAY_F01
        HH_ARRAY_F02 = ARRAY_F02
        HH_ARRAY_RER1 = ARRAY_RER1
        HH_ARRAY_RER2 = ARRAY_RER2
        HH_ARRAY_U1 = ARRAY_U1
        HH_ARRAY_U2 = ARRAY_U2
    else :
        if nom_condition == "LH" :
            LH_ARRAY_C0_m = ARRAY_C0_m
            LH_ARRAY_C1_m = ARRAY_C1_m
            LH_ARRAY_C2_m = ARRAY_C2_m
            LH_ARRAY_R0 = ARRAY_R0
            LH_ARRAY_R1 = ARRAY_R1
            LH_ARRAY_R2 = ARRAY_R2
            LH_ARRAY_V1 = ARRAY_V1
            LH_ARRAY_V2 = ARRAY_V2
            LH_ARRAY_F01 = ARRAY_F01
            LH_ARRAY_F02 = ARRAY_F02
            LH_ARRAY_RER1 = ARRAY_RER1
            LH_ARRAY_RER2 = ARRAY_RER2
            LH_ARRAY_U1 = ARRAY_U1
            LH_ARRAY_U2 = ARRAY_U2  
        else :
            LL_ARRAY_C0_m = ARRAY_C0_m
            LL_ARRAY_C1_m = ARRAY_C1_m
            LL_ARRAY_C2_m = ARRAY_C2_m
            LL_ARRAY_R0 = ARRAY_R0
            LL_ARRAY_R1 = ARRAY_R1
            LL_ARRAY_R2 = ARRAY_R2
            LL_ARRAY_V1 = ARRAY_V1
            LL_ARRAY_V2 = ARRAY_V2
            LL_ARRAY_F01 = ARRAY_F01
            LL_ARRAY_F02 = ARRAY_F02
            LL_ARRAY_RER1 = ARRAY_RER1
            LL_ARRAY_RER2 = ARRAY_RER2
            LL_ARRAY_U1 = ARRAY_U1
            LL_ARRAY_U2 = ARRAY_U2 
    # SETUP POUR EXTRACTION DATA 
    print("_________________________")
print("FIN DES CALCULS _________________________")
#### SOLVING ___________________________________________________________________________________________________________________________________________________________


#### PLOTTING __________________________________________________________________________________________________________________________________________________________
axe_x = np.linspace(0, 300, 16)

HH_puits1 = "#39f559"
LH_puits1 = "#218E34"
LL_puits1 = "#0C3514"

HH_puits2 = "#39e7f1"
LH_puits2 = "#249096"
LL_puits2 = "#114346"

HH_source0 = "#F7ed35"
LH_source0 = "#A49E23"
LL_source0 = "#4C4910"


# RER OBS
plt.plot(range(0,300,10), EXP_RER1_HH, color=HH_puits1, marker="+", label = "RER1_HH")
plt.plot(range(0,300,10), EXP_RER1_LH, color=LH_puits1, marker="+", label = "RER1_LH")
plt.plot(range(0,300,10), EXP_RER1_LL, color=LL_puits1, marker="+", label = "RER1_LL")
plt.title("Evolution du RER du puits 1 d'après les données experimentales", fontsize=12)
plt.ylabel("/°Cj", color="black", fontsize=9)
plt.xlabel("°Cj", color="black", fontsize=9)
plt.legend(fontsize=9)
plt.show()

# UTILISATION OBS
plt.plot(range(0,300,10), EXP_U1_HH, color=HH_puits1, marker="+", label = "U1_HH")
plt.plot(range(0,300,10), EXP_U1_LH, color=LH_puits1, marker="+", label = "U1_LH")
plt.plot(range(0,300,10), EXP_U1_LL, color=LL_puits1, marker="+", label = "U1_LL")
plt.title("Evolution de l'utilisation en C du puits 1 d'après les données experimentales", fontsize=12)
plt.ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
plt.xlabel("°Cj", color="black", fontsize=9)
plt.legend(fontsize=9)
plt.show()

# VOLUME OBS
plt.plot(range(0,300,10), EXP_V1_HH, color=HH_puits1, marker="+", label = "V1_HH")
plt.plot(range(0,300,10), EXP_V1_LH, color=LH_puits1, marker="+", label = "V1_LH")
plt.plot(range(0,300,10), EXP_V1_LL, color=LL_puits1, marker="+", label = "V1_LL")
plt.title("Evolution du volume de puits 1 d'après les données experimentales", fontsize=12)
plt.ylabel("m\u00B3", color="black", fontsize=9)
plt.xlabel("°Cj", color="black", fontsize=9)
plt.legend(fontsize=9)
plt.show()

# RER1 OBS VS RER1 SIM
fig, ax = plt.subplots(1, 2)
ax[0].plot(range(0,300,10), EXP_RER1_HH, color=HH_puits1, marker="+", label = "RER1_HH observé")
ax[0].plot(range(0,300,10), EXP_RER1_LH, color=LH_puits1, marker="+", label = "RER1_LH observé")
ax[0].plot(range(0,300,10), EXP_RER1_LL, color=LL_puits1, marker="+", label = "RER1_LL observé")
ax[1].plot(axe_x, HH_ARRAY_RER1, color=HH_puits2, marker="+", label = "RER1_HH simulé")
ax[1].plot(axe_x, LH_ARRAY_RER1, color=LH_puits2, marker="+", label = "RER1_LH simulé")
ax[1].plot(axe_x, LL_ARRAY_RER1, color=LL_puits2, marker="+", label = "RER1_LL simulé")
ax[0].set_title("Evolution du RER observé dans l'organe puits 1", fontsize=12)
ax[0].set_ylabel("\u03BCmolC/m\u00B3", color="black", fontsize=9)
ax[0].set_xlabel("°Cj", color="black", fontsize=9)
ax[0].legend(fontsize=9)
ax[1].set_title("Evolution du RER simulé dans l'organe puits 1", fontsize=12)
ax[1].set_ylabel("\u03BCmolC/m\u00B3", color="black", fontsize=9)
ax[1].set_xlabel("°Cj", color="black", fontsize=9)
ax[1].legend(fontsize=9)
plt.show()

# U1 OBS VS U1 SIM
fig, ax = plt.subplots(1, 2)
ax[0].plot(range(0,300,10), EXP_U1_HH, color=HH_puits1, marker="+", label = "U1_HH observé")
ax[0].plot(range(0,300,10), EXP_U1_LH, color=LH_puits1, marker="+", label = "U1_LH observé")
ax[0].plot(range(0,300,10), EXP_U1_LL, color=LL_puits1, marker="+", label = "U1_LL observé")
ax[1].plot(axe_x, HH_ARRAY_U1, color=HH_puits2, marker="+", label = "U1_HH simulé")
ax[1].plot(axe_x, LH_ARRAY_U1, color=LH_puits2, marker="+", label = "U1_LH simulé")
ax[1].plot(axe_x, LL_ARRAY_U1, color=LL_puits2, marker="+", label = "U1_LL simulé")
ax[0].set_title("Evolution de l'utilisation observée dans l'organe puits 1", fontsize=12)
ax[0].set_ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
ax[0].set_xlabel("°Cj", color="black", fontsize=9)
ax[0].legend(fontsize=9)
ax[1].set_title("Evolution de l'utilisation simulée dans l'organe puits 1", fontsize=12)
ax[1].set_ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
ax[1].set_xlabel("°Cj", color="black", fontsize=9)
ax[1].legend(fontsize=9)
plt.show()

# 11 OBS VS V1 SIM
fig, ax = plt.subplots(1, 2)
ax[0].plot(range(0,300,10), EXP_V1_HH, color=HH_puits1, marker="+", label = "V1_HH observé")
ax[0].plot(range(0,300,10), EXP_V1_LH, color=LH_puits1, marker="+", label = "V1_LH observé")
ax[0].plot(range(0,300,10), EXP_V1_LL, color=LL_puits1, marker="+", label = "V1_LL observé")
ax[1].plot(axe_x, HH_ARRAY_V1, color=HH_puits2, marker="+", label = "V1_HH simulé")
ax[1].plot(axe_x, LH_ARRAY_V1, color=LH_puits2, marker="+", label = "V1_LH simulé")
ax[1].plot(axe_x, LL_ARRAY_V1, color=LL_puits2, marker="+", label = "V1_LL simulé")
ax[0].set_title("Evolution du volume observée dans l'organe puits 1", fontsize=12)
ax[0].set_ylabel("m\u00B3", color="black", fontsize=9)
ax[0].set_xlabel("°Cj", color="black", fontsize=9)
ax[0].legend(fontsize=9)
ax[1].set_title("Evolution du volume simulée dans l'organe puits 1", fontsize=12)
ax[1].set_ylabel("m\u00B3", color="black", fontsize=9)
ax[1].set_xlabel("°Cj", color="black", fontsize=9)
ax[1].legend(fontsize=9)
plt.show()


# C0_m
plt.plot(axe_x, HH_ARRAY_C0_m, color=HH_source0, marker="+", label = "C0m_HH")
plt.plot(axe_x, LH_ARRAY_C0_m, color=LH_source0, marker="+", linestyle = "dashed", label = "C0m_LH")
plt.plot(axe_x, LL_ARRAY_C0_m, color=LL_source0, marker="+", label = "C0m_LL")
plt.title("C0m simulé", fontsize=12)
plt.ylabel("\u03BCmolC/m\u00B3", color="black", fontsize=9)
plt.xlabel("°Cj", color="black", fontsize=9)
plt.legend(fontsize=9)
# plt.ylim(0, )
# plt.xlim(0, 300)
plt.show()

# R2
plt.plot(axe_x, HH_ARRAY_R2, color=HH_source0, marker="+", label = "R2_HH")
plt.plot(axe_x, LH_ARRAY_R2, color=LH_source0, marker="+", label = "R2_LH")
plt.plot(axe_x, LL_ARRAY_R2, color=LL_source0, marker="+", label = "R2_LL")
plt.title("Evolution de R2", fontsize=12)
plt.ylabel("UNITE A RAJOUTER", color="black", fontsize=9)
plt.xlabel("°Cj", color="black", fontsize=9)
plt.legend(fontsize=9)
plt.show()

# # C1_m et C2_m groupées
# plt.plot(axe_x, HH_ARRAY_C1_m, color=HH_puits1, marker="+", label = "C1m_HH")
# plt.plot(axe_x, LH_ARRAY_C1_m, color=LH_puits1, marker="+", label = "C1m_LH")
# plt.plot(axe_x, LL_ARRAY_C1_m, color=LL_puits1, marker="+", label = "C1m_LL")
# plt.plot(axe_x, HH_ARRAY_C2_m, color=HH_puits2, marker="+", label = "C2m_HH")
# plt.plot(axe_x, LH_ARRAY_C2_m, color=LH_puits2, marker="+", label = "C2m_LH")
# plt.plot(axe_x, LL_ARRAY_C2_m, color=LL_puits2, marker="+", label = "C2m_LL")
# plt.title("Evolution de la concentration en C dans les organes puits 1 et 2", fontsize=12)
# plt.ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
# plt.xlabel("°Cj", color="black", fontsize=9)
# plt.legend(fontsize=9)
# plt.show()

# C1_m et C2_m séparées
fig, ax = plt.subplots(1,2)
ax[0].plot(axe_x, HH_ARRAY_C1_m, color=HH_puits1, marker="+", label = "C1m_HH")
ax[0].plot(axe_x, LH_ARRAY_C1_m, color=LH_puits1, marker="+", label = "C1m_LH")
ax[0].plot(axe_x, LL_ARRAY_C1_m, color=LL_puits1, marker="+", label = "C1m_LL")
ax[1].plot(axe_x, HH_ARRAY_C2_m, color=HH_puits2, marker="+", label = "C2m_HH")
ax[1].plot(axe_x, LH_ARRAY_C2_m, color=LH_puits2, marker="+", label = "C2m_LH")
ax[1].plot(axe_x, LL_ARRAY_C2_m, color=LL_puits2, marker="+", label = "C2m_LL")
ax[0].set_title("Evolution de la concentration en C de l'organe puits 1", fontsize=12)
ax[0].set_ylabel("\u03BCmolC/m\u00B3", color="black", fontsize=9)
ax[0].set_xlabel("°Cj", color="black", fontsize=9)
ax[0].legend(fontsize=9)
ax[1].set_title("Evolution de la concentration en C de l'organe puits 2", fontsize=12)
ax[1].set_ylabel("\u03BCmolC/m\u00B3", color="black", fontsize=9)
ax[1].set_xlabel("°Cj", color="black", fontsize=9)
ax[1].legend(fontsize=9)
plt.show()

# # F01 et F02 groupés
# plt.plot(axe_x, HH_ARRAY_F01, color=HH_puits1, marker="+", label = "F01_HH")
# plt.plot(axe_x, LH_ARRAY_F01, color=LH_puits1, marker="+", label = "F01_LH")
# plt.plot(axe_x, LL_ARRAY_F01, color=LL_puits1, marker="+", label = "F01_LL")
# plt.plot(axe_x, HH_ARRAY_F02, color=HH_puits2, marker="+", label = "F02_HH")
# plt.plot(axe_x, LH_ARRAY_F02, color=LH_puits2, marker="+", label = "F02_LH")
# plt.plot(axe_x, LL_ARRAY_F02, color=LL_puits2, marker="+", label = "F02_LL")
# plt.title("Evolution des flux de C vers les organes puits 1 et 2", fontsize=12)
# plt.ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
# plt.xlabel("°Cj", color="black", fontsize=9)
# plt.legend(fontsize=9)
# plt.show()

# F01 et F02 séparés
fig, ax = plt.subplots(1, 2)
ax[0].plot(axe_x, HH_ARRAY_F01, color=HH_puits1, marker="+", label = "F01_HH")
ax[0].plot(axe_x, LH_ARRAY_F01, color=LH_puits1, marker="+", label = "F01_LH")
ax[0].plot(axe_x, LL_ARRAY_F01, color=LL_puits1, marker="+", label = "F01_LL")
ax[1].plot(axe_x, HH_ARRAY_F02, color=HH_puits2, marker="+", label = "F02_HH")
ax[1].plot(axe_x, LH_ARRAY_F02, color=LH_puits2, marker="+", label = "F02_LH")
ax[1].plot(axe_x, LL_ARRAY_F02, color=LL_puits2, marker="+", label = "F02_LL")
ax[0].set_title("Evolution du flux vers le puits 1", fontsize=12)
ax[0].set_ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
ax[0].set_xlabel("°Cj", color="black", fontsize=9)
ax[0].legend(fontsize=9)
ax[1].set_title("Evolution du flux vers le puits 2", fontsize=12)
ax[1].set_ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
ax[1].set_xlabel("°Cj", color="black", fontsize=9)
ax[1].legend(fontsize=9)
plt.show()

# # RER1 et RER2 groupés
# plt.plot(axe_x, HH_ARRAY_RER1, color=HH_puits1, marker="+", label = "RER1_HH")
# plt.plot(axe_x, LH_ARRAY_RER1, color=LH_puits1, marker="+", label = "RER1_LH")
# plt.plot(axe_x, LL_ARRAY_RER1, color=LL_puits1, marker="+", label = "RER1_LL")
# plt.plot(axe_x, HH_ARRAY_RER2, color=HH_puits2, marker="+", label = "RER2_HH")
# plt.plot(axe_x, LH_ARRAY_RER2, color=LH_puits2, marker="+", label = "RER2_LH")
# plt.plot(axe_x, LL_ARRAY_RER2, color=LL_puits2, marker="+", label = "RER2_LL")
# plt.title("Evolution du RER des puits 1 et 2", fontsize=12)
# plt.ylabel("/°Cj", color="black", fontsize=9)
# plt.legend(fontsize=9)
# plt.show()

# RER1 et RER2 séparés
fig, ax = plt.subplots(1, 2)
ax[0].plot(axe_x, HH_ARRAY_RER1, color=HH_puits1, marker="+", label = "RER1_HH")
ax[0].plot(axe_x, LH_ARRAY_RER1, color=LH_puits1, marker="+", label = "RER1_LH")
ax[0].plot(axe_x, LL_ARRAY_RER1, color=LL_puits1, marker="+", label = "RER1_LL")
ax[1].plot(axe_x, HH_ARRAY_RER2, color=HH_puits2, marker="+", label = "RER2_HH")
ax[1].plot(axe_x, LH_ARRAY_RER2, color=LH_puits2, marker="+", label = "RER2_LH")
ax[1].plot(axe_x, LL_ARRAY_RER2, color=LL_puits2, marker="+", label = "RER2_LL")
ax[0].set_title("Evolution du RER du puits 1", fontsize=12)
ax[0].set_ylabel("/°Cj", color="black", fontsize=9)
ax[0].set_xlabel("°Cj", color="black", fontsize=9)
ax[0].legend(fontsize=9)
ax[1].set_title("Evolution du RER du puits 2", fontsize=12)
ax[1].set_ylabel("/°Cj", color="black", fontsize=9)
ax[1].set_xlabel("°Cj", color="black", fontsize=9)
ax[1].legend(fontsize=9)
plt.show()

# # U1 et U2 groupés
# plt.plot(axe_x, HH_ARRAY_U1, color=HH_puits1, marker="+", label = "U1_HH")
# plt.plot(axe_x, LH_ARRAY_U1, color=LH_puits1, marker="+", label = "U1_LH")
# plt.plot(axe_x, LL_ARRAY_U1, color=LL_puits1, marker="+", label = "U1_LL")
# plt.plot(axe_x, HH_ARRAY_U2, color=HH_puits2, marker="+", label = "U2_HH")
# plt.plot(axe_x, LH_ARRAY_U2, color=LH_puits2, marker="+", label = "U2_LH")
# plt.plot(axe_x, LL_ARRAY_U2, color=LL_puits2, marker="+", label = "U2_LL")
# plt.title("Evolution de l'utilisation en C des puits 1 et 2", fontsize=12)
# plt.ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
# plt.legend(fontsize=9)
# plt.show()

# U1 et U2 séparés
fig, ax = plt.subplots(1, 2)
ax[0].plot(axe_x, HH_ARRAY_U1, color=HH_puits1, marker="+", label = "U1_HH")
ax[0].plot(axe_x, LH_ARRAY_U1, color=LH_puits1, marker="+", label = "U1_LH")
ax[0].plot(axe_x, LL_ARRAY_U1, color=LL_puits1, marker="+", label = "U1_LL")
ax[1].plot(axe_x, HH_ARRAY_U2, color=HH_puits2, marker="+", label = "U2_HH")
ax[1].plot(axe_x, LH_ARRAY_U2, color=LH_puits2, marker="+", label = "U2_LH")
ax[1].plot(axe_x, LL_ARRAY_U2, color=LL_puits2, marker="+", label = "U2_LL")
ax[0].set_title("Evolution de l'utilisation en C du puits 1", fontsize=12)
ax[0].set_ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
ax[0].set_xlabel("°Cj", color="black", fontsize=9)
ax[0].legend(fontsize=9)
ax[1].set_title("Evolution de l'utilisation en C du puits 2", fontsize=12)
ax[1].set_ylabel("\u03BCmolC/°Cj", color="black", fontsize=9)
ax[1].set_xlabel("°Cj", color="black", fontsize=9)
ax[1].legend(fontsize=9)
plt.show()

# # V1 et V2 groupés
# plt.plot(axe_x, HH_ARRAY_V1, color=HH_puits1, marker="+", label = "V1_HH")
# plt.plot(axe_x, LH_ARRAY_V1, color=LH_puits1, marker="+", label = "V1_LH")
# plt.plot(axe_x, LL_ARRAY_V1, color=LL_puits1, marker="+", label = "V1_LL")
# plt.plot(axe_x, HH_ARRAY_V2, color=HH_puits2, marker="+", label = "V2_HH")
# plt.plot(axe_x, LH_ARRAY_V2, color=LH_puits2, marker="+", label = "V2_LH")
# plt.plot(axe_x, LL_ARRAY_V2, color=LL_puits2, marker="+", label = "V2_LL")
# plt.title("Evolution du volume des organes puits 1 et 2", fontsize=12)
# plt.ylabel("m\u00B3j", color="black", fontsize=9)
# plt.legend(fontsize=9)
# plt.show()

# V1 et V2 séparés
fig, ax = plt.subplots(1, 2)
ax[0].plot(axe_x, HH_ARRAY_V1, color=HH_puits1, marker="+", label = "V1_HH")
ax[0].plot(axe_x, LH_ARRAY_V1, color=LH_puits1, marker="+", label = "V1_LH")
ax[0].plot(axe_x, LL_ARRAY_V1, color=LL_puits1, marker="+", label = "V1_LL")
ax[1].plot(axe_x, HH_ARRAY_V2, color=HH_puits2, marker="+", label = "V2_HH")
ax[1].plot(axe_x, LH_ARRAY_V2, color=LH_puits2, marker="+", label = "V2_LH")
ax[1].plot(axe_x, LL_ARRAY_V2, color=LL_puits2, marker="+", label = "V2_LL")
ax[0].set_title("Evolution du volume du puits 1", fontsize=12)
ax[0].set_ylabel("m\u00B3", color="black", fontsize=9)
ax[0].set_xlabel("°Cj", color="black", fontsize=9)
ax[0].legend(fontsize=9)
ax[1].set_title("Evolution du volume du puits 2", fontsize=12)
ax[1].set_ylabel("m\u00B3", color="black", fontsize=9)
ax[1].set_xlabel("°Cj", color="black", fontsize=9)
ax[1].legend(fontsize=9)
plt.show()
# #### PLOTTING __________________________________________________________________________________________________________________________________________________________



    