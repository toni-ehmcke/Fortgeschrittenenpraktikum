#!/usr/bin/python
# -*- coding: utf-8 -*-



import numpy as np
import sys
import matplotlib.pyplot as plt

out = open(sys.argv[3], "w") #Benennung der Ausgabedatei, Öffnen im Schreibmodus

def absorption(trans, reflex): 
	""" Diese Funktion gibt d """
	return 100-reflex-trans
	
def spaces(string):		
	""" Diese Funktion berechnet für eine wohloptimierte den Abstand zwischen zwei Einträgen
	    und gibt diesen zurück """
	numspace = 10 - len(str(string))
	space = " "
	for i in range(0, numspace):
		space += " "
	return space
	
Reflex = np.loadtxt(sys.argv[1])
Trans  = np.loadtxt(sys.argv[2])

for index in range(0, len(Reflex)):				#Anmerkung: Dieser Abschnitt ist sehr speziell und könnte gut verallgemeinert werden.
	index2 = 2 * index
	blau   = absorption(Reflex[index,1], Trans[index2,1])
	#gelb   = absorption(Reflex[index,], Trans[index2,1])
	#orange = absorption(Reflex[index,3], Trans[index2,3])
	#rot    = absorption(Reflex[index,4], Trans[index2,4])
	
	out.write(str(Reflex[index,0]) + "     " + str(blau)+"\n") #+ spaces(blau) + str(gelb)+"\n") #+ spaces(gelb) + str(orange) + spaces(orange) + str(rot)+"\n")
out.close
