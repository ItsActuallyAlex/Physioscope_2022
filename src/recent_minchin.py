#### MODULES
import scipy.optimize as scipy 

import matplotlib.pyplot as plt
import numpy as np
import math

#### MODULES



#### ZONE DE TEST
### TEXTE DE CONTEXTE
# On veut créer une simulation où il y aura une source et deux puits :
# Puits 1 = Tous les organes qui sont en croissance et qaui sont pas hyper interessants
# Puits 2 = Les bourgeons axillaires Z2 et Z3  où l'on veut observer l'accumulation de C ainsi que la baisse de la croissance de l'axe primaire sous traitement LH 

# But = Simuler les Cm dans les organes puits avec des valeurs fixées
### TEXTE DE CONTEXTE


### DEFINITION DE FONCTION
## PHOTOSYNTHESE _________________________________________________________________________________________________________________________
# On par sur des constantes si on a qu'un organe capable de photosynthèse
alpha = 0
Ai_max = 0
def nom_fonction_PHOTOSYNTHESE (PPFD) :
    """explication"""

    photosynthese = (alpha * PPFD * Ai_max)/(alpha * PPFD + Ai_max)
    photosynthese = VALEUR_A_AJOUTER

    return photosynthese 
## PHOTOSYNTHESE _________________________________________________________________________________________________________________________















## PHOTOSYNTHESE _________________________________________________________________________________________________________________________
# Il faut que la taille des arrays K_ et V_ soit identique à celle de l'array des solutions inconnues aka C_
K_ = np.array([1,2])
V_ = np.array([1,2])
C_t0 = np.array([0,0]) 
def nom_fonction_UTILISATION (parametres) :
    """Utilisation du C pour la croissance des deux organes puits"""
    eq = (V_ * C) / (K_ + C)

    return eq
## PHOTOSYNTHESE _________________________________________________________________________________________________________________________

## PHOTOSYNTHESE _________________________________________________________________________________________________________________________
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
## PHOTOSYNTHESE _________________________________________________________________________________________________________________________





## RESISTANCE AU TRANSPORT _________________________________________________________________________________________________________________________
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
## RESISTANCE AU TRANSPORT _________________________________________________________________________________________________________________________


## FLUX _________________________________________________________________________________________________________________________
C_inconnues_t0 = np.array([0,0]) 
C_0 = 0
def nom_fonction_FLUX (C_inconnues, C_0) :
    """Flux de C d'une source vers deux organes puits"""
    delta = R_0*(R_1+R_2)+R_1*R_2

    eq_0 = C_0*((R_2*(C_0-C_inconnues[0]) + R_0*(C_inconnues[1]-C_inconnues[0])) / delta)
    eq_1 = C_0*((R_1*(C_0-C_inconnues[1]) + R_0*(C_inconnues[0]-C_inconnues[1])) / delta)
 
    eq = np.array([eq_0,eq_1])

    return eq
## FLUX _________________________________________________________________________________________________________________________


## EQUILIBRE STOCK MOBILE _________________________________________________________________________________________________________________________
coeff_stockage = 0
coeff_mobilisation = 0
Cs_ = np.array()

def nom_fonction_CARBONE_STOCKE (Cs_) : 
    """Variation du carbone de réserve"""


    return  eq
## EQUILIBRE STOCK MOBILE _________________________________________________________________________________________________________________________






### DEFINITION DE FONCTIONS
#### ZONE DE TEST