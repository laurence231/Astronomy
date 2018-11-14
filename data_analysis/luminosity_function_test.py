import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from ascii_read import *


#matplotlib.rc('text', usetex = True)
matplotlib.rcParams['lines.linewidth'] = 2.0
matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['xtick.major.size'] = 9
matplotlib.rcParams['xtick.minor.size'] = 4
matplotlib.rcParams['xtick.major.width'] = 1.9
matplotlib.rcParams['xtick.minor.width'] = 1.3
matplotlib.rcParams['ytick.major.size'] = 9
matplotlib.rcParams['ytick.minor.size'] = 4
matplotlib.rcParams['ytick.major.width'] = 1.9
matplotlib.rcParams['ytick.minor.width'] = 1.3


def difference_in_bins(positives, negatives, bin_width, bins):
    luminosity_pos = read_col(3, positives)
    luminosity_neg = read_col(3, negatives)
    log_l_pos = np.log10(luminosity_pos)
    log_l_neg = np.log10(luminosity_neg)
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

def define_bins(positives, bin_width):
    luminosity_pos = read_col(3, positives)
    log_l_pos = np.log10(luminosity_pos)
    bins = np.arange(min(log_l_pos), max(log_l_pos), bin_width)
    return bins

def calculate_phi(Difference, volume, bin_width):
    Error_difference = np.sqrt(Difference)
    log_difference = np.log10(Difference / (volume * bin_width))
    log_diff_error_up = np.log10((Difference + Error_difference) / (volume * bin_width)) - np.log10(
        Difference / (volume * bin_width))
    log_diff_error_down = np.log10((Difference) / (volume * bin_width)) - np.log10(
        (Difference - Error_difference) / (volume * bin_width))
    return log_difference, log_diff_error_down, log_diff_error_up



#################
#define the bin width for plots
#################
bin_width = 0.2

###################
# THE DATA FOR SPW3
###################
volume=93.43
positives = 'positive_cube_table_spw3.cat'
negatives = 'negative_cube_table_spw3.cat'
bins = define_bins(positives, bin_width)        #THESE bins must be used invariantly throughout the code.
spw3_bins = difference_in_bins(positives, negatives, bin_width, bins)
spw3_data = calculate_phi(spw3_bins[0], volume, bin_width)

####################
# THE DATA FOR SPW0,1,2 combined
####################
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




#######################
# Plots below here
#######################


CR7_spw3_luminosity_function=plt.errorbar(bins ,spw3_data[0], yerr=[spw3_data[1],spw3_data[2]],xerr=bin_width/2.,marker='s',linestyle='None',color='k', label='CR7')

reference_volume_plot =plt.errorbar(bins,reference_volume_data[0], yerr=[reference_volume_data[1],reference_volume_data[2]],xerr=bin_width/2.,marker='s',linestyle='None',color='g', label='ref vol')


# plt.xlabel('$\log_{10}$(Flux) Jy km s$^{-1}$', fontsize=16)
plt.xlabel('Log 10 (flux) Jy km /s')
# plt.ylabel('$\log_{10}(\Phi/Mpc$^{-3}$)$\,dLogL$^{-1}$',fontsize=16)


#plt.axis([-2.04,-1.2,-4,1.5])
plt.ylabel( 'log10 Phi Mpc-3')
# plt.legend((CR7_LF,CR7_NC_LF,CR7_spw_0_1_2_LF,LAG_LF,POP_LF, pop_3_2, pop_4_3, pop_5_4, pop_6_5, total_line),('This work (CR7 cube) z=6.6','This work (CR7 cube - CR7 clumps)','This work: Reference volumes z \sim 6.6-7.2','Lagache et al. [CII] z=6.7 (prediction)','Popping et al. [CII] z=6.0 (prediction)','Popping CO(3-2) z \sim 0.5', 'Popping CO(4-3) z \sim 1','Popping CO(5-4) z \sim 1.5', 'Popping CO(6-5) z \sim 2', 'Sum of CO and [CII] predictions'))
plt.savefig('LF_full.png', format='png')
plt.show()