#### MODULES
from scipy.optimize import fsolve
#  sympy pour test equa diff
from sympy import *
from scipy.integrate import *
import scipy as scipy 

from autograd import jacobian
import autograd.numpy as np
from numpy import *

from matplotlib import *
import matplotlib.pyplot as plt
#### MODULES




#### TEST FSOLVE
ratio = 1/3.
## CONSTANTES
R_0 = 0.5E13
R_1 = 1.5E13/ratio
R_2 = 1.5E13/(1-ratio)

V_1 = ratio*1E-9
V_2 = (1-ratio)*1E-9

K_1 = 100
K_2 = 100

C_0 = 100
## CONSTANTES


## EQUIVALENCES
a = R_0
b = R_1
c = R_2
d = C_0
g = V_1
h = K_1
i = V_2
j = K_2
## EQUIVALENCES


## X0
# Initial guess avec K_1 << C_1 et K_2 << C_2
# C1_ini = (-a*g - a*i - b*g + d**2)/(d)
# C2_ini = (-a*g - a*i - c*i + d**2)/(d)

# # Initial guess avec K_1 >> C_1 et K_2 >> C_2
C1_ini = (c*d**2*h*i + d**3*h*j)/(a*b*g*i + a*c*g*i + a*d*g*j + a*d*h*i + b*c*g*i + b*d*g*j + c*d*h*i + d**2*h*j)
C2_ini = (b*d**2*g*j + d**3*h*j)/(a*b*g*i + a*c*g*i + a*d*g*j + a*d*h*i + b*c*g*i + b*d*g*j + c*d*h*i + d**2*h*j)

# Initial guess avec C_1 = K_1 et C_2 = K_2
# C1_ini = (-a*g*j - a*h*i - b*g*j + d**2*h*j)/(d*h*j)
# C2_ini = (-a*g*j - a*h*i - c*h*i + d**2*h*j)/(d*h*j)

x_0 = np.array([C1_ini, C2_ini], dtype=float)
x_0 = np.array([0,0])
## X0

delta = R_0*(R_1+R_2)+R_1*R_2

def F (x_0):
    eq_0 = C_0*((R_2*(C_0-x_0[0]) + R_0*(x_0[1]-x_0[0])) / delta)
    eq_1 = C_0*((R_1*(C_0-x_0[1]) + R_0*(x_0[0]-x_0[1])) / delta)
 
    eq = np.array([eq_0,eq_1])

    return eq


def U (x_0) : 
    eq_0 = (V_1 * x_0[0]) / (K_1 + x_0[0])
    eq_1 = (V_2 * x_0[1]) / (K_2 + x_0[1])
 
    eq = np.array([eq_0,eq_1])

    return eq    

# JACOBIENNE / FPRIME
def equa (x_0):
    eq = F(x_0) - U(x_0)
 
    print('EQ', x_0, eq)
    return eq




#### SOLVING
x_result = scipy.optimize.fsolve(equa, x0=x_0, col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)
print('x_results :',x_result)

## Random prints
print('Flux :',F(x_result))
print('Utilisation :',U(x_result))
print('Equation Eéquilibre :',equa(x_result))
## Random prints
#### SOLVING




















# #### TEST SCIPY.OPTIMIZE
# # Dérivées partielles 
# x[0], x[1] = symbols('x y')

# # Use sympy.diff() method
# df1_dx = diff(eq_inconnue_1, x[0])
# df1_dy = diff(eq_inconnue_1, x[1])
# df2_dx = diff(eq_inconnue_2, x[0])
# df2_dy = diff(eq_inconnue_2, x[1])

# # Poser la matrice jacobienne
# jac = np.empty(shape=(0, 0))
# jac = jac.reshape(2, 2)

# # Dérivées partielles 

# scipy.optimize.newton(func=[eq_inconnue_1,eq_inconnue_2], x0=x_0, fprime=None, args=(), tol=1.48e-08, maxiter=50, fprime2=None, x1=None, rtol=0.0, full_output=False, disp=True)
# #### TEST SCIPY.OPTIMIZE















# #### ALGO METHODE DE NEWTON EQUATION LINEAIRES ET NON LINEAIRES
# # Définition des équations du système
# func_1 = lambda x: ((R_2*(C_0-x[0]) + R_0*(x[1]-x[0])) / (R_0*(R_1+R_2)+R_1*R_2))*C_0 - (V_1 * x[0]) / (K_1 + x[0])
# func_2 = lambda x: ((R_1*(C_0-x[1]) + R_0*(x[0]-x[1])) / (R_0*(R_1+R_2)+R_1*R_2))*C_0 - (V_2 * x[1]) / (K_2 + x[1])



# jac_func1 = jacobian(func_1)
# jac_func2 = jacobian(func_2)

# # variable itérative
# i = 0
# # erreur
# error = 100
# # tolérance pour la solution
# tol = 0.0001
# # itérations max
# maxiter = 100
# # paramètres de dimensions des matrices
# # M le nombre d'équations en jeu
# M = 2
# # N le nombre d'inconnues
# N = 2

# ## INITIAL GUESS
# # Array avec notre initial guess ici x[0] et x[1] = 1

# # Initial guess avec K_1 << C_1 et K_2 << C_2
# C1_ini = (-a*i - a*g - b*g + d**2)/(d)
# C2_ini = (-a*g - a*i - c*i + d**2)/(d)

# # Initial guess avec K_1 >> C_1 et K_2 >> C_2
# # C1_ini = (c*d**2*h*i + d**3*h*j)/(a*b*g*i + a*c*g*i + a*d*g*j + a*d*h*i + b*c*g*i + b*d*g*j + c*d*h*i + d**2*h*j)
# # C2_ini = (b*d**2*g*j + d**3*h*j)/(a*b*g*i + a*c*g*i + a*d*g*j + a*d*h*i + b*c*g*i + b*d*g*j + c*d*h*i + d**2*h*j)

# # Initial guess avec C_1 = K_1 et C_2 = K_2
# # C1_ini = (-a*g*j - a*h*i - b*g*j + d**2*h*j)/(d*h*j)
# # C2_ini = (-a*g*j - a*h*i - c*h*i + d**2*h*j)/(d*h*j)


# # Création de l'array avec le couple d'initial guess
# x_0 = np.array([C1_ini, C2_ini], dtype=float).reshape(N,1)
# ## INITIAL GUESS

# # np.any teste une relation et retourne TRUE ou FALSE
# # ici n pose erreur > tolérance
# while np.any(abs(error) > tol) and i < maxiter :
    
#     # On évalue nos fonctions au guess intial et on stocke les valeurs dans fun_evaluate
#     fun__evaluate = np.array([func_1(x_0), func_2(x_0)]).reshape(M,1)
#     # fun_evaluate en fait c'est F(X)

    # Maintenant on veut créer la matric Jacobienne F'(X)
    # on fait une liste plate des valeurs initiales
    # flat_x_0 = x_0.flatten()

    # jac_func1 = jacobian(func_1)
    # jac_func2 = jacobian(func_2)

    # jac = np.array([jac_func1(flat_x_0), jac_func2(flat_x_0)])
    # jac = jac.reshape(N, M)
    # jac c'est notre matrice jacobienne

#     # Donc on a tout il reste plus qu'a appliquer la formule
#     x_new = x_0 - np.linalg.inv(jac)@fun__evaluate
#     # np.linalg.inv() permet d'avoir l'inverse de la totalité d'une matrice
#     # le @ est en fait l'opérateur permettant de multiplier des matrices entre elles

#     # On update la valeur de l'erreur
#     error = x_new - x_0
#     # C'est facile à comprendre car c'est la longueur de chaque pas et plus elle devient petite plus on se rapproche de la solution
#     x_0 = x_new

#     # Facultatif, utile pour checker les itérations
#     # print(i)
#     # print(error)
#     # print("-----------")

#     # On update l'itérateur 
#     i = i + 1

# # printing de la solution
# print("the solution is")
# print(x_new)

# # vérifications 
# print("Eq-1", np.around(func_1(x_new),3))
# print("Eq-2", np.around(func_2(x_new),3))
# # np.around fait un arrondi à 3 décimales près
# #### ALGO METHODE DE NEWTON EQUATION LINEAIRES ET NON LINEAIRES















# #### Dérivées partielles PRINCIPE
# x, y = symbols('x y')

# func = x + y

# print("Before Differentiation : {}".format(func))
  
# # Use sympy.diff() method
# func = diff(func, x)

# print("After Differentiation : {}".format(func))
# #### Dérivées partielles PRINCIPE