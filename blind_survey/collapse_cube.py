import numpy
import astropy.io.fits as pyfits

def cube_collapse(CUBENAME, velocity_width_collapse, sign):
	'''
	This file contains the function which will collapse the cube, with the user input defining if it is positive or negative, and inputting
	the name of the cube
	:param CUBENAME: the reduced ALMA fits cube    velocity_width_collapse: This was set to 74 km s-1 for CR7 sign: MUST = +- 1
	:return: the script produces image files with velocity width defined by the user.
	'''

	if sign==1:
		print("collapsing the positive mapping")
	elif sign == -1:
		print("collapsing the negative mapping")
	else:
		print("you have not put in a valid value for the sign! (+-1)")
		break
	beam_area_pixels=300.19 #pixelscale 0.05
	pixscale=0.05
	Lya_redshift=6.606
	############################################################################################
	# This first section of the script reads in the value of delta_v to be read in as velocity_range
	# which is the array read into the for loop seen below. This means that the velocity width
	# boundaries are read in as an input, which enables the code to self determine the start
	# point for where to begin 'slicing' from
	############################################################################################
	data=pyfits.getdata(CUBENAME)
	hd=pyfits.getheader(CUBENAME)



	freq=numpy.arange(float(len(data[0,:,:,:])))

	freq=freq*hd['CDELT3']+hd['CRVAL3']
	c=299792.458 #speed of light in km/s
	redshift=((1.90054E12)/freq)-1.
	delta_v=c*((1+redshift)/(1+Lya_redshift)-1)
	print 'Frequencies in the cube:', freq
	print 'Correspond do velocity offsets:', delta_v
	#################
	velocity_range = delta_v
	# This for loop will repeat the cut of velocity_width_collapse from each value of delta_v that exists in the cube.
	for maximum_velocity in velocity_range: #km/s

		#maximum_velocity=1289. #+100 km/s w.r.t. Lya redshift

		#FIRST FROM COLLAPSED FIG FROM FREDERIQUE
		data=pyfits.getdata(CUBENAME)
		hd=pyfits.getheader(CUBENAME)



		freq=numpy.arange(float(len(data[0,:,:,:])))

		freq=freq*hd['CDELT3']+hd['CRVAL3']
		c=299792.458 #speed of light in km/s
		redshift=((1.90054E12)/freq)-1.
		delta_v=c*((1+redshift)/(1+Lya_redshift)-1)
		print 'Frequencies in the cube:', freq
		print 'Correspond do velocity offsets:', delta_v

		############
		#just an alternative way to compute velocity, doesnt do anything at this moment
		z =Lya_redshift
		freqcii =  1.90054e+03/(1.+z)
		frest=freqcii*1e9
		velo = -299792.458 * ((freq-frest)/freq) #another way of writing it
		############
		dv=abs(delta_v[0]-delta_v[1])
		print 'Velocity resolution', dv

		number_of_slices_to_collapse=int(round(velocity_width_collapse/dv,0))
		print 'We will collapse # slices:',number_of_slices_to_collapse

		print delta_v
		print numpy.argmin(abs(maximum_velocity-delta_v))
		array_index_maximum_velocity=numpy.argmin(abs(maximum_velocity-delta_v))
		print delta_v[array_index_maximum_velocity]

		array_index_minimum_velocity=array_index_maximum_velocity+number_of_slices_to_collapse

		print 'Now collapsing between velocities',delta_v[array_index_maximum_velocity],delta_v[array_index_minimum_velocity]


		data=sign*data[:,array_index_maximum_velocity:array_index_minimum_velocity,:,:] #Frederic used: 25:43
		#print numpy.shape(data)


		collapse=dv*numpy.nansum(data,axis=1)
		collapse=collapse[0,:,:]
		print numpy.shape(collapse)

		pyfits.writeto('CR7_collapse_dv_%skms_%skms.fits'%(int(delta_v[array_index_maximum_velocity]),int(delta_v[array_index_minimum_velocity])),collapse,hd,clobber=True)
