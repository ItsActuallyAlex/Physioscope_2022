#### MODULES
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math

#### MODULES


## CONSTANTES
R_0 = 0.5E13
#C_0 = 100

ratio = 0.33

R_1 = 1.5E13/ratio
K_1 = 100
V_1 = ratio*1E-9

R_2 = 1.5E13/(1-ratio)
K_2 = 100
V_2 = (1-ratio)*1E-9
#### CONSTANTES


#### FONCTIONS
def F (x_0, C_0) :
    """Flux de C d'une source vers deux organes puits"""
    delta = R_0*(R_1+R_2)+R_1*R_2

    eq_0 = C_0*((R_2*(C_0-x_0[0]) + R_0*(x_0[1]-x_0[0])) / delta)
    eq_1 = C_0*((R_1*(C_0-x_0[1]) + R_0*(x_0[0]-x_0[1])) / delta)
 
    eq = np.array([eq_0,eq_1])

    return eq


def U (x_0) : 
    """Utilisation du C"""
    eq_0 = (V_1 * x_0[0]) / (K_1 + x_0[0])
    eq_1 = (V_2 * x_0[1]) / (K_2 + x_0[1])
 
    eq = np.array([eq_0,eq_1])

    return eq  


def equa (x_0, C_0) :
    """Equation d'intérêt"""
    eq = F(x_0, C_0) - U(x_0)

    return eq
    

def rapportF1F2 (x_0, C_0) :
    """Rapport F1/F2"""
    return F(x_0, C_0)[0]/F(x_0, C_0)[1]
#### FONCTIONS

## CONDITIONS INITIALES

# Initial guess avec K_1 << C_1 et K_2 << C_2
#C1_ini = (-R_0*V_1 - R_0*V_2 - R_1*V_1 + C_0**2)/(C_0)
#C2_ini = (-R_0*V_1 - R_0*V_2 - R_2*V_2 + C_0**2)/(C_0)

# # Initial guess avec K_1 >> C_1 et K_2 >> C_2
#C1_ini = (R_2*C_0**2*K_1*V_2 + C_0**3*K_1*K_2)/(R_0*R_1*V_1*V_2 + R_0*R_2*V_1*V_2 + R_0*C_0*V_1*K_2 + R_0*C_0*K_1*V_2 + R_1*R_2*V_1*V_2 + R_1*C_0*V_1*K_2 + R_2*C_0*K_1*V_2 + C_0**2*K_1*K_2)
#C2_ini = (R_1*C_0**2*V_1*K_2 + C_0**3*K_1*K_2)/(R_0*R_1*V_1*V_2 + R_0*R_2*V_1*V_2 + R_0*C_0*V_1*K_2 + R_0*C_0*K_1*V_2 + R_1*R_2*V_1*V_2 + R_1*C_0*V_1*K_2 + R_2*C_0*K_1*V_2 + C_0**2*K_1*K_2)

# Initial guess avec C_1 = K_1 et C_2 = K_2
#C1_ini = (-R_0*V_1*K_2 - R_0*K_1*V_2 - R_1*V_1*K_2 + C_0**2*K_1*K_2)/(C_0*K_1*K_2)
#C2_ini = (-R_0*V_1*K_2 - R_0*K_1*V_2 - R_2*K_1*V_2 + C_0**2*K_1*K_2)/(C_0*K_1*K_2)
#x_0 = np.array([C1_ini, C2_ini], dtype=float)

x_0 = np.array([0, 0], dtype=float)
## CONDITIONS INITIALES


#### SOLVING
## LOOP POUR C_0
range_C_0 = np.arange(0, 1000, 25)
valeurs_F1 = np.empty(0)
valeurs_F2 = np.empty(0)
F1F2 = np.empty(0)


for C_0 in range_C_0 :

    x_result = scipy.fsolve(equa, x0=x_0, args=(C_0), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)

    ## CREATION ARRAY FLUX
    fluxes = F(x_result, C_0)
    valeurs_F1 = np.append(valeurs_F1, fluxes[0])
    valeurs_F2 = np.append(valeurs_F2, fluxes[1])
    F1F2 = np.append(F1F2, rapportF1F2(x_result, C_0))
    ## CREATION ARRAY FLUX

    # ## RANDOM PRINTS
    # print('valeur C_0 :',C_0)    
    # print('x_results :',x_result)
    # print('Flux 1 et 2 :',fluxes)
    # print('Utilisation 1 et 2 :',U(x_result))
    # print('Equation équilibre :',equa(x_result, C_0))

    # print('F1/F2', rapportF1F2(x_result, C_0))
    # print('______________________')
    # ## RANDOM PRINTS



#### ZONE DE TEST
### TEXTE DE CONTEXTE
# On veut créer une simulation où il y aura une source et deux puits :
# Puits 1 = Tous les organes qui sont en croissance et qaui sont pas hyper interessants
# Puits 2 = Les bourgeons axillaires Z2 et Z3  où l'on veut observer l'accumulation de C ainsi que la baisse de la croissance de l'axe primaire sous traitement LH 

# But = Simuler les Cm dans les organes puits avec des valeurs fixées
### TEXTE DE CONTEXTE


### DEFINITION DE FONCTION
# Il faut que la taille des arrays K_ et V_ soit identique à celle de l'array des solutions inconnues aka C_
K_ = np.array([1,2])
V_ = np.array([1,2])
C_t0 = np.array([0,0]) 
def nom_fonction_UTILISATION (parametres) :
    """Utilisation du C pour la croissance des deux organes puits"""
    eq = (V_ * C) / (K_ + C)

    return eq



V_t0 = 0
V_ = np.array([1,2])
RER = 0 
area = 0
epsilon = 0
def nom_fonction_VOLUME_ORGANE (area) :
    """Variation de l'utilisation en fonction du volume"""
    Vi = epsilon * area
    dVi_dt = Vi + RER 

    return dVi_dt

EVOLUTION_UTILISATION = scipy.fsolve(nom_fonction_VOLUME_ORGANE, x0=Vi_t0, args=(), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)



# On par sur des constantes si on a qu'un organe capable de photosynthèse
alpha = 0
Ai_max = 0
def nom_fonction_PHOTOSYNTHESE (PPFD) :
    """explication"""

    photosynthese = (alpha * PPFD * Ai_max)/(alpha * PPFD + Ai_max)
    photosynthese = VALEUR_A_AJOUTER

    return photosynthese 



# constante_resistance à 20°C
# constante des gaz parfaits R = 8.314J/mol/K
# longueurs = l'ensemble des longueurs entre les deux organes
# rayons = les rayons internes des conduits vasculaires
# Minchin proposent une valeur de viscosité dépendant de la concentration en solutés
viscosity = VALEUR_A_AJOUTER 
temp_20 = 293
gaz_p = 8.314


longueur_ = VALEUR_A_AJOUTER
rayon_ = VALEUR_A_AJOUTER

def nom_fonction_RESISTANCES (longueur_, rayon_) :
    """Calcul du jeu de resistances du modèle - Minchin et al. 1993"""
    # viscosité fixe ? Si non besoin de rajouter une fonction de fitting pour viscosity
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 
    resistances = constante_resistance * (longueur_)/(rayon_**4)

    # return un array de resistances 
    return resistances


C_inconnues_t0 = np.array([0,0]) 
C_0 = 0
def nom_fonction_FLUX (C_inconnues, C_0) :
    """Flux de C d'une source vers deux organes puits"""
    delta = R_0*(R_1+R_2)+R_1*R_2

    eq_0 = C_0*((R_2*(C_0-C_inconnues[0]) + R_0*(C_inconnues[1]-C_inconnues[0])) / delta)
    eq_1 = C_0*((R_1*(C_0-C_inconnues[1]) + R_0*(C_inconnues[0]-C_inconnues[1])) / delta)
 
    eq = np.array([eq_0,eq_1])

    return eq


coeff_stockage = 0
coeff_mobilisation = 0
Cs_ = np.array()

def nom_fonction_CARBONE_STOCKE (Cs_) : 
    """Variation du carbone de réserve"""


    return  eq







### DEFINITION DE FONCTIONS
#### ZONE DE TEST










## LOOP POUR C_0
#### SOLVING

#### PLOTTING
# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.plot(range_C_0, valeurs_F1, color="blue", marker="+")
ax.plot(range_C_0, valeurs_F2, color="green", marker="+")
# set x-axis label
ax.set_xlabel("Concentration C_0", fontsize=14)
# set y-axis label
ax.set_ylabel("Flux F1 et F2", color="black", fontsize=14)

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(range_C_0, F1F2, color="black", linestyle= "-.")
ax2.set_ylabel("F1/F2", color="black", fontsize=14)
ax2.set_ylim([0, 1])
plt.show()
#### PLOTTING

