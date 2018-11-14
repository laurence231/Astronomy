import numpy
import math
import astropy.io.fits as pyfits
from matplotlib import pyplot
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
import random
from ascii_read import *

import sys

def cii_flux(COLLAPSEFILE, noise):
	sources=ascii.read("sources.dat")
	print >>open('randomrealisations.cat','w'),'#x y flux lum'

	data=pyfits.getdata(COLLAPSEFILE)
	hd=pyfits.getheader(COLLAPSEFILE)
	beam_area_pixels=300.19 #pixelscale 0.05

	X=len(data[0,:])
	Y=len(data[:,0])
	pixscale=0.05
	radius=0.35/pixscale

	Xrandom=numpy.arange(radius,X-radius,1.)
	Yrandom=numpy.arange(radius,Y-radius,1.)

	FLUXES=[]

	for i in range(noise):
		xcenter=random.choice(Xrandom)
		ycenter=random.choice(Yrandom)

		y,x=numpy.ogrid[-ycenter:Y-ycenter,-xcenter:X-xcenter]
		circle_mask=(x*x+y*y<=radius**2)
		totflux=numpy.nansum(data[circle_mask])/beam_area_pixels
		FLUXES.append(totflux)
		print >>open('randomrealisations.cat','a'),xcenter,ycenter,totflux,totflux*1143195867.5586417/1E8


	####################################################################
	#Below is the cii flux script which              #
	#we want to use to calculate the cii flux given an input from a    #
	#text file of x and y coordinates, and refer to a slice which will #
	#refer to the collapse file variable we are considering.           #
	####################################################################

	xcenter=sources['col1'] 		#these centers open the x and y columns from the sources table
	ycenter=sources['col2']

	FLUX=[]

	radius=0.35/pixscale 				#0.35 arcsec is one beam radius
	for i in range(len(xcenter)):
		y,x=numpy.ogrid[-ycenter[i]:Y-ycenter[i],-xcenter[i]:X-xcenter[i]]
		circle_mask=(x*x+y*y<=radius**2)
		#data[-circle_mask]=numpy.nan

		totflux=numpy.nansum(data[circle_mask])/beam_area_pixels
		FLUX.append(totflux)
		B=totflux*1143195867.5586417/1E8
		C=totflux/numpy.std(FLUXES)
		D=10**(((numpy.log10(B*1E8)-5.720095873)/1.18438764))
		print >> open('slice_table.cat','a'), xcenter[i], ycenter[i], totflux, B, C, D
