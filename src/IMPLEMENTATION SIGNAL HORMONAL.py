# IMPLEMENTATION SIGNAL HORMONAL


plus on a d'utilisation 
plus la resistance baisse
l'utilisation est liée au RER
donc RER first
ensuite UTILISATION
ensuite baisse RESISTANCE


def RER (Ci_m, v_i, k_i, signal) :
    
    eq = ((v_i * Ci_m)/(k_i + Ci_m))*signal

    return eq


if nom_condition == "HH" :
    signal = np.array([1, 1])
else :
    if nom_condition == "LH" :
        signal = np.array([1, 0.5])

    else :
        signal = np.array([1, 0.25])

