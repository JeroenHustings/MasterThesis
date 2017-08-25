import os
import numpy as np
import matplotlib.pyplot as plt 
from sys import argv
import satlas as s 


# massnumber of the isptope
mass =47
#scan_number of the isotope
scan_no_1 = 135

# Fitting doesn't work well with big numbers, so use an offset
# the first number is the NIST value of the transition in MHz 
#https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=K+I&low_wl=&upp_wn=&upp_wl=&low_wn=&unit=1&submit=Retrieve+Data&de=0&java_window=3&java_mult=&format=0&line_out=0&en_unit=0&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&order_out=0&max_low_enrg=&show_av=2&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on
OFFSET = 389286074.578447

# speed of ligth
c = 299792458

# function for the doppler shift (look it up the thesis)
def doppler(V,m):
	beta = np.sqrt(1-((m* 931.494095*10**6)**2/(V+m* 931.494095*10**6)**2))
	dop = np.sqrt((1-beta)/(1+beta))
	return dop


# directory of the file you want to fit
directory= 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_online_binned\\' 
#directory where you want to save the results
directory_save = 'D:\\Agi\\Documents\\PhD\\CRIS\\K\\K data\\K_online_results\\iscool_calibration\\'
# these two lines creat the directory_save if the folder doesn't exist
if not os.path.exists(directory_save):
    os.makedirs(directory_save)
	

m = 47	   # add the correct mass, as many decimal palces as possible (for example from Audi database: https://www-nds.iaea.org/amdc/ame2016/mass16.txt)
A_L = 12   # these values you can look up in Data_exl.csv
B_L = 12   # these values you can look up in Data_exl.csv
A_U = 12   # these values you can look up in Data_exl.csv
B_U = 12   # these values you can look up in Data_exl.csv

# loops over files in the forlder. You can play a bit with this to understand how it works
# for example you can print the filename 
for filename in os.listdir(directory):
	if filename.endswith(".txt"): # checks if the file is .txt

		# using numpy to load in data from a .txt file
		data = np.loadtxt(directory + filename, skiprows=1, delimiter = ',')
		# I leave it for you to figure out what this one does
		scan_no = filename.strip('.txt')

		# checks if the file corresponds r=to the scan number you are looking for
		if    int(scan_no)==scan_no_1: 

			# making lists with data you want to use in your fit
			wavelength=data[:,2]
			wavelength_err=data[:, 3]
			counts=data[:,13]
			n = data[:,6]
			rate= counts/n
			v_iscool = data[:,7]
			yerr = np.sqrt(counts)
			yerr[yerr == 0] = 1
			yerr = yerr / n*2.3
			freq= c*wavelength*10**(-4)
			frequency  = freq * doppler(v_iscool, m) # doppler shift the frequency


	
			I=0 # add here the correct spin
			J=[0.5,0.5] # J number in the case of K 0.5, 0.5
			ABC=[A_L, A_U, B_L, B_U,0, 0] # A and B parameters, C is 0. That would be the octupole deformation

			centroid =3500 # you can change this nu,ber to get the correct centroid. Once you have a plot, 
							# you can eastimate where the centroid is
			frequency -= OFFSET # making the fitt easier using small numbers

			scale=np.max(rate)
			bkg=[rate[0]]
			scale -= bkg[0]

			fwhm=[20,20] 
			# satlas part: http://woutergins.github.io/satlas/tutorial.html

			spectrum = s.HFSModel(I=I,J=J,ABC=ABC,centroid=centroid,fwhm = fwhm, 
			                      background_params =[0.01], scale = scale, 
			                      use_racah=False
			                      , shape = 'voigt')

			spectrum.set_variation({'Al':False, 'Au':False,'Cl':False,'Cu':False, 'FWHMG':True, 'FWHML':True})
			spectrum.plot(x=frequency,y=rate,yerr=yerr)
			#do the fit!
			s.chisquare_fit(spectrum,x=frequency,y=rate,yerr=yerr)

			# show the fit
			spectrum.display_chisquare_fit()
			print(spectrum.get_result_frame(vary=True))


			# plotting data
			plt.errorbar(x=frequency,y=rate,yerr=yerr, fmt='bo')
			plt.plot(frequency,spectrum(frequency), 'g-', linewidth=2.5)
			plt.xlabel('MHz')
			plt.ylabel('Rate')
			plt.title('{}'.format(filename))
			plt.show()


# saving the results in a txt file:

			# scan_no = filename.strip('.txt')
			# with open(directory_save+'{}.txt'.format(scan_no), 'w') as f:
			# 	for p in spectrum.params.values():
			# 		f.write('{}\t{}\t{}\n'.format(p.name,p.value,p.stderr))
			# 	f.write('{}\t{}\n'.format('dof',spectrum.ndof_chi))
			# 	f.write('{}\t{}\n'.format('chisqr',spectrum.chisqr_chi))
			# 	f.write('{}\t{}\n'.format('redchisqr',spectrum.redchi_chi))


