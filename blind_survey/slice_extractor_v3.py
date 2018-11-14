import os 
from astropy.io import fits as pyfits
from ascii_read import *
from __main__ import *
from collapse_cube import *
from noise_and_ciiV3 import *


def test_for_empty_slice(filename):
	"""
	This function is an ugly way of checking that the slice contains meaningful data, by testing
	to see if two random pixels near the center are zero valued or not. Some slices come out blank if at extremities
	so this was the easiest way to check in my eyes. This could be improved 
	"""
	imagefile = pyfits.open('%s'%(filename))    #the purpose of this is to test 
	pixel = imagefile[0].data
	region = pixel[500,550]
	region2 = pixel[412,528]
	if region != 0 and region2 !=0:
		return True
	else:
		return False

def sextractor(sign, noise, redshift, luminosity_distance):
	'''
	This function will apply sextractor to each of the slices that are produced and within the current directory. It will output
	full tables of the resulting analysis.
	It requires this script to be in a directory with
	-the cube for analysis
	-the sextractor config and param files, as well as any filters if used
	-the python scripts which define the functions imported and used within this script.
	:param sign:
	:return: A table with all of the excellent data for analysis!
	'''
	print >> open('%s_cube_table.cat'%(sign), 'w'), '# X Y Flux Luminosity S/N SFR ra dec filename '
	for filename in os.listdir('.'):
		if filename.endswith("kms.fits") == True:
			
			is_valid_slice = test_for_empty_slice(filename)
    		
    		if is_valid_slice == True and filename.endswith("kms.fits") == True:
				print >> open('slice_table.cat', 'w'), '# X Y Flux Luminosity S/N SFR '
				print(filename)
				os.system('sex -c config.sex %s ' % (filename))
				catalog = 'test.cat'

				ra = read_col(1, catalog)  # These all correspond to the sextractor detections
				dec = read_col(2, catalog)
				x = read_col(7, catalog)
				y = read_col(8, catalog)

				print >> open('sources.dat', 'w')
				for i in range(len(x)):
					print >> open('sources.dat', 'a'), x[i], y[i]

				cii_flux(filename, noise, redshift, luminosity_distance)

				analysis = 'slice_table.cat'

				x_a = read_col(1,analysis)  # Use a subscript to denote its from the internship analysis script (the good data)
				y_a = read_col(2, analysis)
				flux_a = read_col(3, analysis)
				lum_a = read_col(4, analysis)
				thresh_a = read_col(5, analysis)
				sfr_a = read_col(6, analysis)

				for i in range(len(x_a)):
					print >> open('%s_cube_table.cat'%(sign), 'a' ), x_a[i], y_a[i], flux_a[i], lum_a[i], thresh_a[i], sfr_a[
						i], ra[i], dec[i], filename


os.system('sudo mkdir positive negative')
cubename = str(input("please enter the name of the ALMA cube (enclosed in ' '):", ))
output_name = str(input('what output name do you want for the catalogue (e,g, the spw number enclosed in ' '): ', ))
collapse_width = input("please enter the collapse width (CR7 we used 74km/s): ", )
noise = input('how many iterations for noise? (we use 5000 for study, could use lower for testing): ', )
redshift = input('Enter the redshift value for the center of the cube',  )
luminosity_distance = input('Enter the luminosity distance (use ned wrights cosmo calculator):' , )

cube_collapse(cubename, collapse_width,1)
sextractor('positive', noise, redshift, luminosity_distance)                        #collapse cube, SExtract
os.system('sudo mv positive_cube_table.cat %s_positives.cat'%(output_name))

for filename in os.listdir('.'):          
	if filename.endswith("kms.fits"):     #important file management - stops positives and negatives getting mixed up
		os.system('sudo mv %s positive/'%(filename))

cube_collapse(cubename, collapse_width,-1)
sextractor('negative', noise, redshift, luminosity_distance)
os.system('sudo mv negative_cube_table.cat %s_negatives.cat'%(output_name))

for filename in os.listdir('.'):
	if filename.endswith("kms.fits"):
		os.system('sudo mv %s negative/'%(filename))