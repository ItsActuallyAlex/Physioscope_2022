#### MODULES

#### MODULES



#### CONDITIONS INITIALES

#### CONDITIONS INITIALES



#### CONSTANTES
# Constantes pour la résistance
temp_20 = 293
gaz_p = 8.314

# Constantes pour le RER
v_ = VALEUR_
k_ = VALEUR_

# Pour l'utilisation
delta = VALEUR_
"Est-ce que delta sera bien une moyenne (donc une seule valeur par traitement) ?"

# Pour la variation de la masse en C
surface_ = VALEUR_ 
"J'ai oublié l'approche exacte pour la mesure"

# On garde un volume fixe moyen en regardant les organes qui ne sont plus en extension ?
volume_fixe = VALEUR_



#### CONSTANTES



#### FONCTIONS
def PHOTOSYNTHESE_complete (alpha, PPFD, Ai_max, P) :
    """Equation de la photosynthèse, complète"""

    eq = (alpha * PPFD * Ai_max)/(alpha * PPFD + Ai_max) * P 

    return eq


def PHOTOSYNTHESE_simple (VALEUR_PHOTOSYNTHESE) :
    """Equation de la photosynthèse, simplifiée"""

    eq = VALEUR_PHOTOSYNTHESE 

    return eq


def RESISTANCE_complete (viscosity, longueur_, rayon_) :
    """Equation de la résistance au transport (Minchin et al. 1993), complète"""

    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    eq = constante_resistance * (longueur_)/(rayon_**4)

    return eq


# Ci_m obtenu par la résolution du système à l'équilibre au pas de temps d'avant
def INCREMENT_VOLUME_simple (Ci_m) :
    """Equation pour incrément dvi_dt, simplifiée"""

    # 
    eq = (v_ * Ci_m)/(k_ + Ci_m)

    return eq


def UTILISATION_simple (param) :
    """Equation de l'utilisation, simplifiée"""

    eq = delta * INCREMENT_VOLUME_simple (param)

    return eq 


def FLUX_complet (param) :
    """Equations solutionnées comme précédemment"""

    R_ = RESISTANCE_complete (viscosity, longueur_, rayon_)

    denominateur = R_[1]*(R_[2]+R_[3])+R_[2]*R_[3]

    "Les x vont se jumeler avec C_0 pour n'avoir à la fin qu'un seul vecteur "
    eq_0 = C_0*((R_[3]*(C_0-x_0[0]) + R_[1]*(x_0[1]-x_0[0])) / denominateur)
    eq_1 = C_0*((R_[2]*(C_0-x_0[1]) + R_[1]*(x_0[0]-x_0[1])) / denominateur)
 
    eq = np.array([eq_0,eq_1])



def VOLUME_simple (param) :
    """Variation du volume des organes, simplifiée"""

    eq = Vi * INCREMENT_VOLUME_simple (param)
    
    return eq

def MASSE_MOBILE_simple (param) :
    """Variation de la masse mobile en C dans un organe, simplifiée"""

    eq = FLUX + PHOTOSYNTHESE_simple (VALEUR_PHOTOSYNTHESE) * surface_ - UTILISATION_simple (param)

    return eq


def CONCENTRATION_MOBILE_simple (param) : 
    """Concentration en C mobile dans un organe sans la variation en volume, simplifiée"""

    eq = (MASSE_MOBILE_simple (param)) / (volume_fixe) 

    return eq

def RESOLUTION_MINCHIN (param) :
    """equation à résoudre """

    eq = FLUX + PHOTOSYNTHESE_simple (VALEUR_PHOTOSYNTHESE) - UTILISATION_simple (param)

    return eq
#### FONCTIONS



#### RESOLUTION
Cm_i = scipy.fsolve(RESOLUTION_MINCHIN, x0=x_0, args=(param), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)
#### RESOLUTION








































































































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
# P = Durée d'éclairement journalière
P = 0
def nom_fonction_PHOTOSYNTHESE (PPFD, area) :
    """Taux de production de photosynthétats"""
    # Formule complète
    photosynthese = (alpha * PPFD * Ai_max)/(alpha * PPFD + Ai_max) * P * area 

    # Etape 1 = Constante 
    photosynthese = VALEUR_A_AJOUTER

    return photosynthese 
## PHOTOSYNTHESE _________________________________________________________________________________________________________________________



## RESISTANCE AU TRANSPORT ______________________________________________________________________________________________________
# viscosité fixe ? 
viscosity = VALEUR_A_AJOUTER 
#  temp = température en K pour 20°C
temp_20 = 293
# R = constante des gaz parfaits
gaz_p = 8.314
# longueur_ = l'ensemble des longueurs entre les deux organes
longueur_ = np.array([VALEUR_A_AJOUTER])
# rayon_ = les rayons internes des conduits vasculaires
rayon_ = np.array([VALEUR_A_AJOUTER])

def nom_fonction_RESISTANCES (longueur_, rayon_) :
    """Calcul du jeu de resistances au transport dans le modèle"""
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 
    resistances = constante_resistance * (longueur_)/(rayon_**4)

    # return un array de resistances 
    return resistances
## RESISTANCE AU TRANSPORT ______________________________________________________________________________________________________


 
## FLUX _________________________________________________________________________________________________________________________
C_inconnues_t0 = np.array([VALEUR_A_AJOUTER]) 
C_0 = 0
# R_[1] = R_0
# R_[2] = R_1
# R_[3] = R_2
R_ = nom_fonction_RESISTANCES (longueur_, rayon_)

def nom_fonction_FLUX (C_inconnues, C_0) :
    """Flux de C d'une source vers deux organes puits"""
    delta = R_[1]*(R_[2]+R_[3])+R_[2]*R_[3]

    # Sans doute pas la bonne forme car le modèle sera pas en râteau 
    eq_0 = C_0*((R_[3]*(C_0-C_inconnues[0]) + R_[1]*(C_inconnues[1]-C_inconnues[0])) / delta)
    eq_1 = C_0*((R_[2]*(C_0-C_inconnues[1]) + R_[1]*(C_inconnues[0]-C_inconnues[1])) / delta)
 
    eq = np.array([eq_0,eq_1])

    return eq
## FLUX _________________________________________________________________________________________________________________________



## DILUTION Cm_ _________________________________________________________________________________________________________________
def nom_fonction_DILUTION_Cm_ (parametres_maybe) :
    """explications"""

    eq = fonction_variation_masse / dVi_dt

    return eq
## DILUTION Cm_ _________________________________________________________________________________________________________________



## UTILISATION POUR LA CROISSANCE _______________________________________________________________________________________________
# # Il faut que la taille des arrays K_ et V_ soit identique à celle de l'array des solutions inconnues aka C_
# K_ et V_ pas utilisées pour l'instant
K_ = np.array([1,2])
V_ = np.array([1,2])
Cm_ = np.array([VALEUR_A_AJOUTER])
delta = VALEUR_A_AJOUTER
def nom_fonction_UTILISATION (parametres_maybe) :
    """Utilisation du C pour la croissance des deux organes puits"""

    # Etape 1 = Fitting sur données 
    eq = np.array([VALEUR_A_AJOUTER])

    # dVi_dt sera un array qu'on aura obtenu depuis la croissance en volume
    eq = delta * dVi_dt

    return eq
## UTILISATION POUR LA CROISSANCE _______________________________________________________________________________________________



## CROISSANCE EN VOLUME _________________________________________________________________________________________________________
# V_t0 = Volume de l'organe à t0
Vi_t0 = 0
# RER sera un array obtenu par la data
# Michaelienne pour le RER pour le moment, estimer "v et K" 
RER = np.array([VALEUR_A_AJOUTER])
# area sera un array obtenu par la data 
area = np.array([VALEUR_A_AJOUTER])
# Epsilon facteur de conversion obtenu par des mesures experimentales d'une feuille 
epsilon = VALEUR_A_AJOUTER
def nom_fonction_VOLUME_ORGANE (Volume, signaux, C) :
    """Variation de l'utilisation en fonction du volume"""
    Vi = epsilon * area
    eq = Vi * RER 

    return eq

dVi_dt = scipy.fsolve(nom_fonction_VOLUME_ORGANE, x0=Vi_t0, args=(), col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)
## CROISSANCE EN VOLUME _________________________________________________________________________________________________________



## VARIATION CARBONE STOCK ______________________________________________________________________________________________________
coeff_stockage = 0
coeff_mobilisation = 0
Cs_ = np.array([VALEUR_A_AJOUTER])
Cm_ = np.array([VALEUR_A_AJOUTER])
def nom_fonction_CARBONE_STOCKE (Cs_, Cm_) : 
    """Variation du carbone de réserve"""
    eq = coeff_stockage * Cm_ - coeff_mobilisation * Cs_

    return  eq
## VARIATION CARBONE STOCK ______________________________________________________________________________________________________



## VARIATION MASSE EN C MOBILE __________________________________________________________________________________________________
# gamma = respiration de maintenance
gamma = 0
# On le met à 0 pour cette étape
EQUILIBRE_STOCK_MOBILE = 0
def nom_fonction_MASSE_Cm_ (parametres_maybe) :
    """Variation de la masse en C mobile au sein d'un organe"""
    eq = FLUX + PHOTOSYNTHESE * SURFACE + VOLUME * (EQUILIBRE_STOCK_MOBILE) - UTILISATION - gamma

    eq = FLUX + PHOTOSYNTHESE * SURFACE - UTILISATION

    return eq 
## VARIATION MASSE EN C MOBILE __________________________________________________________________________________________________



### DEFINITION DE FONCTIONS
#### ZONE DE TEST