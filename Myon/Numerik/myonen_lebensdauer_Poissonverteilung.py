# -*- coding: utf-8 -*-
"""Dieses Programm wendet die Likelihood-Methode auf die Datei mdat.txt an, bei der die relevanten Daten in Zeile 85 beginnen.
Die Lebensdauer wird anhand der Poissonverteilung bestimmt. Siehe hierzu in der Anleitung zum Versuch Myonen mehr.
Das Minimum der Kurve gibt Zerfallszeit an.
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

data = np.loadtxt("mdat.txt", skiprows=85)					# Daten fangen erst in Zeile 85 an

N_ges = data[:,1].sum()									    # Anzahl aller Ereignisse
t1 = data[0,0]												# floatIndex erster Kanal
tk = data[-1,0]												# floatIndex letzter Kanal
dt = 1./24													# Zeitschritt zwischen Kanal in micro Sekunden		
tau = np.linspace(0.1,3.4,1000)								# Kandidaten für die Lebensdauer
															# Zeit in Mikrosekunden
															# L hat gleiche Länge wie tau
print "dt =%s"%(tau[1]-tau[0])
def fetchChannelTime(channel, dt):
	"""function to change easily the computation of channel-number -> middle channel time"""
#	print channel * dt + 0.75
	return channel * dt + 0. 

def getN0(N_ges, t1, tk, dt, tau):
	""" bestimmt Normierungskonstante der Exponentialfunktion nach Glg. 14 in der Anleitung Myonen"""
	return N_ges / ( np.e**(-fetchChannelTime(t1, dt)/tau) - np.e**(-(fetchChannelTime(tk, dt)+dt)/tau))

def L_eff(getN0, data, t1, tk, dt, tau):
	"""
	getN0(N_ges, t1, tk, dt, tau) = Funktion die aus Gesamtanzahl der Ereignisse N_ges, 
									Anfangs- bzw Endkanal t1 und tk und Binbreite dt 
									Die Normierungskonstante N0 bestimmt (glg. 14 Skript) 
	data = 2dim array mit Kanälen und zugehörigen Zählereignissen Nk
		tau = momentaner Schätzwert für die Lebensdauer"""
	N0 = getN0(N_ges, t1, tk, dt, tau)					# Daniel arbeitet mit Fkt get_N0_2
	Nk = data[:, 1]
	k =	data[:, 0]										# bei Daniel geht k von 0 bis 255 hier von 1 bis 256							
	L = 0.
	for i in range(Nk.size):
		N_i = Nk[i]
		t_i = fetchChannelTime(k[i], dt)				# Kanäle in Zeiten umrechnen 
		f_i = N0/tau * np.e**(-(t_i+ dt/2.)/tau) * dt     # bei Daniel im exponent t_i+ 0.75 + dt/2 
		L += -2 * N_i * np.log(f_i) 							  # Glg. 22 im Skript	
	return L
def fehler(getN0, data, t1, tk, dt, tau):
	"""Bestimmt die Messunsicherheit auf Leff nach Gausscher Fehlerfortpflanzung auf Glg. 22"""
	N0 = getN0(N_ges, t1, tk, dt, tau)					# Daniel arbeitet mit Fkt get_N0_2
	Nk = data[:, 1]
	k =	data[:, 0]										# bei Daniel geht k von 0 bis 255 hier von 1 bis 256							
	dL_tilde = 0.
	for i in range(Nk.size):
		N_i = Nk[i]
		t_i = fetchChannelTime(k[i], dt)				# Kanäle in Zeiten umrechnen 
		f_i = N0/tau * np.e**(-(t_i+ dt/2.)/tau) * dt     # bei Daniel im exponent t_i+ 0.75 + dt/2 
		dL_tilde += (2* np.log(f_i))**2 * N_i 			# Glg. 22 im Skript fehler auf N_i ist np.sqrt(N_i)	
	return np.sqrt(dL_tilde)


def main():
	L = L_eff(getN0, data, t1, tk, dt, tau)
	dl = fehler(getN0, data, t1, tk, dt, tau)			# Messunsicherheit nach Fehlerfortpflanzung
	# Berechne Lebensdauer und Unsicherheut
	tau_mean = tau[L.argmin()]								  # Mittelwert der Lebensdauer
	tau_max = tau[(L+dl).argmin()]
	tau_min = tau[(L-dl).argmin()]
	error = max(np.abs(tau_mean-tau_min), np.abs(tau_mean-tau_max))# Maximale Abweichung ist Fehler
	print("Lebensdauer: tau[L.argmin()]=%fe-6 s"%tau[L.argmin()])
	print "Messunsicherheit dtau = %se-6 s"%error
	fig = plt.figure("exp")
	ax = plt.subplot(111)
	ax.set_title(r"Poisson $\tau = %s \pm %s \mu s $"%(round(tau_mean, 3), round(error, 3)))
	ax.plot(tau, L)
	ax.plot(tau, L+dl, ls ="--", color="blue")			# Abweichung nach oben	
	ax.plot(tau, L-dl, ls ="--", color="blue")
	# Beschriftung
	ax.text(2.75, -400000., r'$2\cdot\sum(-N_i\ \mathrm{ln}f_i)$')
	ax.set_xlabel(r"Lebensdauer in $\mu s$ ")
	ax.set_ylabel(r"$\mathrm{L}_\mathrm{eff}$")	
	# mehr Achsenstriche
	majorLocator = MultipleLocator(1)								   # große Striche in 1 Abständen
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(0.25)
	ax.xaxis.set_major_locator(majorLocator)
	ax.xaxis.set_major_formatter(majorFormatter)
	# for the minor ticks, use no labels; default NullFormatter
	ax.xaxis.set_minor_locator(minorLocator)

	plt.draw()
	plt.show()

main()
