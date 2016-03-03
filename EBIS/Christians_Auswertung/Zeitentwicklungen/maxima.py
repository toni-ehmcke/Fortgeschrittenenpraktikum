#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys

out = open("maxima.txt", "w")

for mdat in range(1, len(sys.argv)):
    l = len(sys.argv[mdat])
    spacenumber = 12 - l
    spaces = " "                         #erzeuge einheitliche, übersichliche Tabellenstruktur
    for character in range(0, spacenumber+1):   #erzeuge dynamisch viele Leerzeichen zwischen den Tabelleneinträgen
        spaces += " "
    
    data = np.loadtxt(sys.argv[mdat]) #lade die Messdaten, die als Argumente an das Skript übergeben wurden
    

    sum = 0
    for i in range(0, len(data)):  #summiere alle Einträge der Ladungszahl
        sum += data[i,2]
        
    Nmax = max(data[:,2])   #suche Maxium
    w = Nmax/sum            #normiere die Maxima
    
    index = 0               #finde die Ionisationszeit zum Maximum    
    while data[index,2] <> Nmax:
        index += 1
    t_ion = data[index,1]
    spacesnumber = 5 - len(str(t_ion))
    spacess = " "
    for character in range(0, spacesnumber+1):
        spacess += " "
    
    out.write(sys.argv[mdat] + spaces + str(t_ion) + spacess + str(w)+"\n") #gibt Daten in Datei aus
        
    print(sys.argv[mdat],"W = %f , Zeit = %f, N_max = ", w, t_ion, Nmax)
out.close()