#### MODULES
import numpy as np
import math
import scipy.optimize as scipy
#### MODULES

#### CONSTANTES QUI VARIENT PAS DU TOUT
viscosity = 1E6
temp_20 = 293
gaz_p = 8.314
#### CONSTANTES QUI VARIENT PAS DU TOUT

#### CONSTANTES QUI VARIENT ENTRE TRAITEMENTS
photosynthese_journaliere = np.array([429859.7489, 641689.2869, 237700.9616], dtype=float)

longueur_entrenoeuds = np.array([2.91047327, 2.91047327, 2.91047327], dtype=float)
rayon_entrenoeuds = np.array([0.160215922,0.160215922, 0.160215922], dtype=float)/100

volumes_puits1 = np.array([1.79192E-06, 1.79192E-06, 1.79192E-06], dtype=float)

C0_m_t0_Sucrose_fitted = np.array([360,300,300], dtype=float)
C0_m_t0_Sucrose_moyen = np.array([327.2154541, 327.2154541, 327.2154541], dtype=float)
C1_m_t0 = np.array([518.3472222, 479.0671296, 479.0671296], dtype=float)
C2_m_t0 = np.array([0.069444444, 0.069444444, 0.069444444], dtype=float)

k = np.array([405.197, 397.561, 357.624], dtype=float)
v = np.array([220, 180, 180], dtype=float)

delta_masse_volume = np.array([1.03453E+13, 1.03453E+13, 1.03453E+13], dtype=float)
#### CONSTANTES QUI VARIENT ENTRE TRAITEMENTS

#### FONCTIONS DU MODELE
def PHOTOSYNTHESE_SOURCE () :
    return photosynthese_journaliere


def RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds) :
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    resistances_modele = constante_resistance * (longueur_entrenoeuds)/(rayon_entrenoeuds**4)

    resistances_modele = np.append(resistances_modele, (resistances_modele, resistances_modele))
    
    return resistances_modele 


def FLUX_CARBONE (longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m) :
    R_ = RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds)

    denominateur_delta = R_[0]*(R_[1] + R_[2]) + R_[1]* R_[2]

    F01 = (R_[2]*(C0_m - C1_m) + R_[0]*(C2_m - C1_m)) / denominateur_delta
    F02 = (R_[1]*(C0_m - C1_m) + R_[0]*(C1_m - C2_m)) / denominateur_delta

    return [F01, F02]
    

def VOLUME_PUITS (volumes_puits1) :
    return volumes_puits1 


def RER_PUITS (v, k, C1_m, C2_m) :

    RER_puits1 = (v*C1_m)/(k+C1_m)
    RER_puits2 = (v*C2_m)/(k+C2_m)

    return [RER_puits1, RER_puits2]


def UTILISATION_CARBONE (volumes_puits1, delta_masse_volume) :

    utilisation_puits1 = delta_masse_volume * VOLUME_PUITS(volumes_puits1)
    utilisation_puits2 = utilisation_puits1/100

    return [utilisation_puits1, utilisation_puits2]


def A_RESOUDRE (longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m, delta_masse_volume, volumes_puits1) :

    C1_m_equilibre = C0_m*FLUX_CARBONE(longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m) - UTILISATION_CARBONE(delta_masse_volume, volumes_puits1)
    C2_m_equilibre = C0_m*FLUX_CARBONE(longueur_entrenoeuds, rayon_entrenoeuds, C0_m, C1_m, C2_m) - UTILISATION_CARBONE(delta_masse_volume, volumes_puits1)
    return [C1_m_equilibre, C2_m_equilibre]
    

# print(UTILISATION_CARBONE(delta_masse_volume, volumes_puits1))
    
# print(A_RESOUDRE (longueur_entrenoeuds, rayon_entrenoeuds, C0_m_t0_Sucrose_fitted, C1_m_t0 , C2_m_t0, delta_masse_volume, volumes_puits1))

for condition in range(3) :
    # = np.array([photosynthese_journaliere[condition]])
    A = np.array([longueur_entrenoeuds[condition]])
    B = np.array([rayon_entrenoeuds[condition]])
    # C0_m_t0_Sucrose_fitted = np.array([C0_m_t0_Sucrose_fitted[condition]])
    C = np.array([C0_m_t0_Sucrose_moyen[condition]])
    X = np.array([C1_m_t0[condition]])
    Y = np.array([C2_m_t0[condition]])
    D = np.array([volumes_puits1[condition]])
    # = np.array([k[condition]])
    # = np.array([v[condition]])
    E = np.array([delta_masse_volume[condition]])
    
    # X0 = np.append(X, Y)
    args = A,B,C,D,E
    
    equilibre_concentrations = scipy.fsolve(A_RESOUDRE, x0=(X, Y), args=(args))






    # photosynthese_journaliere = np.array([429859.7489, 641689.2869, 237700.9616], dtype=float)
    longueur_entrenoeuds = np.array([2.91047327, 2.91047327, 2.91047327], dtype=float)
    rayon_entrenoeuds = np.array([0.160215922,0.160215922, 0.160215922], dtype=float)/100
    volumes_puits1 = np.array([1.79192E-06, 1.79192E-06, 1.79192E-06], dtype=float)
    C0_m_t0_Sucrose_fitted = np.array([360,300,300], dtype=float)
    C0_m_t0_Sucrose_moyen = np.array([327.2154541, 327.2154541, 327.2154541], dtype=float)
    C1_m_t0 = np.array([518.3472222, 479.0671296, 479.0671296], dtype=float)
    C2_m_t0 = np.array([0.069444444, 0.069444444, 0.069444444], dtype=float)
    # k = np.array([405.197, 397.561, 357.624], dtype=float)
    # v = np.array([220, 180, 180], dtype=float)
    delta_masse_volume = np.array([1.03453E+13, 1.03453E+13, 1.03453E+13], dtype=float) 

#### FONCTIONS DU MODELE