import numpy
import numpy as np
from astropy.io import fits as pyfits
from scipy.optimize import curve_fit
import math
from scipy.integrate import quadpack
import matplotlib.pyplot as plt    #imports the matplotlib.pyplot package for the actual graph. Renames it plt for easier use later
from ascii_read import *
import matplotlib

####################################
# Gaussian Function
def gaussian(x, x0, sigma):
    return np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
####################################

def double_normal(phi, err_down, err_up, size):
    '''
Double gaussian, asymmetric
phi_dummy generates the space of values to pick from
phi_probs generates the PDF for a 1sigma error around each phi
    '''
    phi_dummy = np.linspace(phi - 5. * err_down, phi + 5. * err_up, 100000)
    phi_probs = np.append(gaussian(phi_dummy[phi_dummy < phi], phi, err_down),
                          gaussian(phi_dummy[phi_dummy >= phi], phi, err_up))
    return np.random.choice(phi_dummy, size, p=phi_probs / np.sum(phi_probs))

def difference_in_bins(positives, negatives, bin_width, bins):
    luminosity_pos = read_col(4, positives)
    luminosity_neg = read_col(4, negatives)
    log_l_pos = np.log10(luminosity_pos) + 8.
    log_l_neg = np.log10(luminosity_neg) + 8.
    number = bins - bins
    for i in range(len(bins)):
        for j in range(len(log_l_pos)):
            if bins[i] - bin_width / 2. < log_l_pos[j] < bins[i] + bin_width / 2.:
                number[i] = number[i] + 1
    print(number)
    number_neg = bins - bins
    for i in range(len(bins)):
        for j in range(len(log_l_neg)):
            if bins[i] - bin_width / 2. < log_l_neg[j] < bins[i] + bin_width / 2.:
                number_neg[i] = number_neg[i] + 1
    print(number_neg)
    Difference = number - number_neg
    print(bins)
    for i in range(len(number)):
        if number[i] - number_neg[i] > 0:
            Difference[i] = number[i] - number_neg[i]
        else:
            Difference[i] = 1
    return Difference, number, number_neg

def calculate_phi(Difference, volume, bin_width):
    Error_difference = np.sqrt(Difference)
    log_difference = np.log10(Difference / (volume * bin_width))
    log_diff_error_up = np.log10((Difference + Error_difference) / (volume * bin_width)) - np.log10(
        Difference / (volume * bin_width))
    log_diff_error_down = np.log10(Difference/(volume * bin_width)) - np.log10(
        (Difference - Error_difference) / (volume * bin_width))
    return log_difference, log_diff_error_down, log_diff_error_up

def compute_overdensity(bins, phi, err_up, err_down, ref_phi, ref_err_down, ref_err_up ):
    iterations = 1000
    med_overdensity = numpy.zeros(len(bins))  # This is the median overdensity in the bins instead of overall
    med_overdensity_erup = numpy.zeros(len(bins))
    med_overdensity_erdown = numpy.zeros(len(bins))
    print >> open('file_check.cat', 'w')
    NEW_CR7 = numpy.zeros(len(phi))
    NEW_REF = numpy.zeros(len(ref_phi))
    for j in range(len(bins)):
        if np.isnan(phi[j]) == False and np.isinf(phi[j]) == False and np.isnan(err_down[j]) == False and np.isinf(err_down[j]) == False and np.isnan(err_up[j]) == False and np.isinf(err_up[j]) == False:
            Overdensity = []
            for i in range(iterations):
                NEW_CR7[j] = double_normal(phi[j], err_down[j], err_up[j], 1)[0]
                NEW_REF[j] = double_normal(ref_phi[j], ref_err_down[j], ref_err_up[j], 1)[0]
                Overdensity.append(pow(10., NEW_CR7[j]) / pow(10., NEW_REF[j]))
                print >> open('file_check_%s.cat' % (j), 'a'), NEW_CR7[j], NEW_REF[j], \
                (pow(10., NEW_CR7) / pow(10., NEW_REF))[0]
            med_overdensity[j] = numpy.median(Overdensity)  # This is the median overdensity in the bin [j]
            med_overdensity_erup[j] = numpy.percentile(Overdensity, 84) - numpy.median(Overdensity)
            med_overdensity_erdown[j] = numpy.median(Overdensity) - numpy.percentile(Overdensity, 16)
        else:
            print("there was NAN's or INFS, take a look at the data please!")
            STOP
    print med_overdensity, 'Determined overdensity'
    print med_overdensity_erup, 'error up'
    print med_overdensity_erdown, 'error down'
    print pow(10., phi) / pow(10., ref_phi), 'Expected Overdensity spw3'
    return med_overdensity, med_overdensity_erup, med_overdensity_erdown

#####################
#This is where the bins are configured
#####################
bin_width = 0.3
bins=np.array([7.15])

positives = 'positive_cube_table_spw3.cat'
negatives = 'negative_cube_table_spw3.cat'
volume=93.43       #cube volume - leave in units of Mpc3

spw3_bins = difference_in_bins(positives, negatives, bin_width, bins)
spw3_data = calculate_phi(spw3_bins[0], volume, bin_width)

spw2_vol = 95 #ish, this needs to be updated
positives = 'positive_cube_table_spw2.cat'
negatives = 'negative_cube_table_spw2.cat'
spw2_bins = difference_in_bins(positives, negatives, bin_width, bins)
print('this is spw2 bins', spw2_bins[0])

spw1_vol = 100
positives = 'positive_cube_table_spw1.cat'
negatives = 'negative_cube_table_spw1.cat'
spw1_bins = difference_in_bins(positives, negatives, bin_width, bins)
print('this is spw1 bins', spw1_bins[0])

spw0_vol = 101
positives = 'positive_cube_table_spw0.cat'
negatives = 'negative_cube_table_spw0.cat'
spw0_bins = difference_in_bins(positives, negatives, bin_width, bins)
print('this is spw0 bins', spw0_bins[0])

reference_volume_difference = spw2_bins[0] + spw1_bins[0] + spw0_bins[0]
volume_of_reference_cubes = spw2_vol + spw1_vol + spw0_vol
print('this is CR7 bins', spw3_bins[0])
print('this is the ref vol difference: ',reference_volume_difference)
reference_volume_data = calculate_phi(reference_volume_difference, volume_of_reference_cubes, bin_width)

overdensity = compute_overdensity(bins, spw3_data[0], spw3_data[1], spw3_data[2], reference_volume_data[0], reference_volume_data[1], reference_volume_data[2])

print(overdensity)