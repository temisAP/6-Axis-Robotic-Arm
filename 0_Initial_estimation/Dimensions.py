import numpy as np
from math import pi
import matplotlib.pyplot as plt

""" Define components """

class element():
    def __init__(self,description):
        self.description = description

class components():
    def __init__(self):

        self.NEMA   = element('NEMA 17 HS4401')
        self.NEMA.weight = 280e-3 * 9.81  # kg -> N
        self.NEMA.torque = 0.4 * 0.75     # N·m

        self.servo  = element('Micro servo')
        self.servo.weight = 16e-3 * 9.81  # kg -> N

        self.bar    = element('A bar of steal')
        self.bar.diameter = 8e-3          # m
        self.bar.rho      = 7850 * 9.81   # kg/m³ > N/m³
        self.bar.lamb     = self.bar.rho * (self.bar.diameter/2)**2 * pi

        self.polea  = element('Correa y polea')
        self.polea.weight = np.array([0.042, 0.042]) * 9.81  # kg -> N [small side, big side]

        self.weight = element('The weight carried by the hand')
        self.weight.weight = 1*9.81  # kg -> N

        self.delta = 0.6 # Adimensional distance where 4 motor is placed
        self.f2 = 3      # Performance factor of motor 2
        self.f3 = 1      # Performance factor of motor 3

c = components()

""" One motor dimensioning """
from scipy.optimize import fsolve

def equations(p):

    Dse = p
    Dwe = Dse

    # Weights for each position
    Ww = c.weight.weight + c.servo.weight    + c.polea.weight[0]
    We = c.polea.weight[0] + 2*c.polea.weight[1]
    W4 = c.NEMA.weight   + c.polea.weight[0]

    # Equations
    eq1 = Ww * (Dwe + Dse) + We * Dse- c.NEMA.torque * c.f2

    return eq1

Dse = fsolve(equations, (1))[0]
print('')
print(f'For Dwe = Dse the solution is')
print(f'    Dwe = {(Dse*1e3):.2f} mm')
print(f'    Dse = {Dse*1e3:.2f} mm')
print('')


""" Define equation system """


def equations(p):

    Dwe, Dse = p

    # Weights for each position
    Ww = c.weight.weight + c.servo.weight    + c.polea.weight[0]
    We = c.NEMA.weight   + c.polea.weight[0] + c.polea.weight[1]
    W4 = c.NEMA.weight   + c.polea.weight[0]

    # Equations
    eq1 = Ww * (Dwe + Dse) + We * Dse           + c.bar.lamb * (Dwe**2/2 + (1+c.delta) * Dse**2/2) - c.NEMA.torque * c.f2
    eq2 = Ww * Dwe         - W4 * Dse * c.delta + c.bar.lamb * (1-c.delta)* 2*Dse**2               - c.NEMA.torque * c.f3 *0.1

    return (eq1, eq2)


""" Solve and represent in function of NEMA reduction system """

deltas = np.linspace(0,1,11)
f2s = np.linspace(1,51,4)
f3s = np.linspace(1,51,10)

fig = plt.figure()
for marker,f2 in enumerate(f2s):
    for color,f3 in enumerate(f3s):
        c.f2 = f2
        c.f3 = f3
        Dwe, Dse = fsolve(equations, (0.1, 0.1))
        print('')
        print(f'For f2 = {f2:.2f} and f3 = {f3:.2f} the solution is')
        print(f'    Dwe = {Dwe:.2f} m')
        print(f'    Dse = {Dse:.2f} m')
        print('')
        plt.scatter(Dwe,Dse,marker=f'{marker+1}',label=f'f2 = {f2:.2f} | f3 = {f3:.2f}', linewidths=1,c=plt.cm.jet(np.linspace(0,1,len(f3s))[color]))

plt.xlabel('Dwe')
plt.ylabel('Dse')
plt.legend()
plt.grid()
plt.show()



c.f2 = 51
c.f3 = 51

Dwe, Dse =  fsolve(equations, (0.5, 0.5))
print('***** From graphics *****')
print(f'For f2 = {c.f2:.2f} and f3 = {c.f3:.2f} the solution is')
print(f'Dwe = {Dwe:.2f} m')
print(f'Dse = {Dse:.2f} m')
print('**********')


""" Maximize arm size """

from scipy.optimize import minimize

def arm_size(p):
    f2, f3 = p
    c.f2 = f2
    c.f3 = f3
    Dwe, Dse = fsolve(equations, (0.5, 0.5))
    if Dse < 0 or Dwe < 0:
        return 0
    else:
        return -(Dse+Dwe)


x0 = [10,51]
bnds = ((51, 51), (1, 1))
f2, f3 = minimize(arm_size, x0, method='nelder-mead',bounds=bnds, options={'xatol': 1e-8, 'disp': True}).x

c.f2 = f2
c.f3 = f3

Dwe, Dse =  fsolve(equations, (0.5, 0.5))
print('***** From optimization *****')
print(f'For f2 = {c.f2:.2f} and f3 = {c.f3:.2f} the solution is')
print(f'Dwe = {Dwe:.2f} m')
print(f'Dse = {Dse:.2f} m')
print('**********')
