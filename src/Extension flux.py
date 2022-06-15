from sympy import *

q1, q2, q3, r0, r1, r2, r3 = symbols('q1 q2 q3 r0 r1 r2 r3')

#### EQUATIONS
# s01 = (q1+q2+q3)*r0 + q1*r1
s02 = (q1+q2+q3)*r0 + q2*r2
s03 = (q1+q2+q3)*r0 + q3*r3


s01 = (q1+q2)*r0 + q1*r1
# s02 = (q1+q2)*r0 + q2*r2
#### EQUATIONS


print(solve(s01, (q1, q2)))

