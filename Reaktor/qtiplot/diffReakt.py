# -*- coding: utf-8 -*-
"""
Dieses Programm berechnet die differentielle Reaktivität am AKR-2 in DD
zu gegebenen Messwerten für die Verdopplungszeit.
"""

import numpy as np

T_2 = np.array([127,97,72,96,121])              # Verdopplungszeit
T_s = T_2/np.log(2)                             # stabile Reaktorperiode
# Zerfallskonstante für beteiligte Mutterkerne
zerfKonst = np.array([0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01])
# relativer Anteil der verzögerten Neutronen der einzelnen Mutterkerne
# bzgl. des Gesamtanteils der verzögerten Neutronen
a = np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])

# parameter, der Verhältnis aus mittlerer Lebensdauer der prompten Neutronen
# und dem Anteil der verzögerten Neutronen ist (keine phys. Bedeutung)
c = 0.0051 # in Sekunden
dT_2 = 1.                                   # Fehler auf Verdopplungszeit in s
dT_s = dT_2/np.log(2)                       # Fehler auf Reaktorperiode in s

rho = np.zeros(5)                           # differentielle Reaktivität in $
drho = np.zeros(5)                          # Fehler auf diff Reaktivität in $
for i in range(5): 
    rho[i] = c/T_s[i] + np.sum(a/(1+zerfKonst*T_s[i]))
    drho[i] = (c/T_s[i]**2 + np.sum(a*zerfKonst/(1+zerfKonst*T_s[i])**2))*dT_s
print rho, drho

dIntRho = np.zeros(5)
for i in range(5):
    dIntRho[i] = np.sqrt(np.sum(drho[:i+1])**2)
    print drho[:i+1]
print dIntRho


