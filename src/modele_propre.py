#### MODULES
#### MODULES

#### CONSTANTES QUI VARIENT PAS DU TOUT
viscosity =
temp_20 =
gaz_p = 
#### CONSTANTES QUI VARIENT PAS DU TOUT

#### CONSTANTES QUI VARIENT ENTRE TRAITEMENTS
photosynthese_journaliere =
longueur_entrenoeuds =
rayon_entrenoeuds = 
#### CONSTANTES QUI VARIENT ENTRE TRAITEMENTS

#### FONCTIONS DU MODELE
def PHOTOSYNTHESE_SOURCE () :
    return photosynthese_journaliere


def RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds) :
    constante_resistance =  (8*viscosity)/(math.pi * temp_20 * gaz_p) 

    resistances_modele = constante_resistance * (longueur_entrenoeuds)/(rayon_entrenoeuds**4)
    
    return resistances_modele 


def FLUX_CARBONE (longueur_entrenoeuds, rayon_entrenoeuds) :
    resistances_modele = RESISTANCES_TRANSPORT (longueur_entrenoeuds, rayon_entrenoeuds)

    denominateur_delta = 

    F01 = 
    F02 = 

    return [F01, F02]
    

def VOLUME_PUITS () :
    return volumes_puits


def RER_PUITS (v, k, Ci_m) :

    RER_puits = (v*Ci_m)/(k+Ci_m)

    return RER_puits
"PAS UTILISE POUR L'INSTANT"

def UTILISATION_CARBONE (delta) :

    utilisation_puits = delta * VOLUME_PUITS()

    return utilisation_puits

def A_RESOUDRE () :
    return equilibre_concentrations 
    

for condition in range(3) :
    constante_variable [condition]
    
    RESOLUTION via fsolve
#### FONCTIONS DU MODELE