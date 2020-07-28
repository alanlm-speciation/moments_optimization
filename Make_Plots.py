# Library importation
import sys
import os
import getopt
import numpy
from numpy import array
import moments
import pylab
import matplotlib
from datetime import datetime
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import Model_2pop_folded
import Model_2pop_unfolded


#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================




#Help function
def usage():
    """ Function for help """
    print("# This script allow you to plot the results of the best model infer from previous model comparison using moments on python 3 \n"+
    "# It uses the plotting routine scripted by Portick et al, 2017 and adapted to moments by Le Moan et al., 2020 \n\n"+
    "# This is an exemple of the most complete command line :\n"+
    "# python Script_moments_model_optimisaton.py -f SNP_file -y Pop1 -x Pop2 -p 30,30 -m IM_b -s True -r 1 -z \n\n"+
    "# that will performed the optimisaton for the IM model with a folded SFS (-s True) with singleton masked (-z) obtained from the file named SNP_file.data \n\n"+
    "# -h --help : Display the help you are looking at.\n"+
    "# -f SNPs input data file with .data as extention \n"+
    "# -y --population1 : Take the name of the first population in the sfs (y-axis)\n"+
    "# -x --population2 : Take the name of the second population in the sfs (x-axis)\n"+
    "# -p --projection : tak two number for the size of the projecition"+
    "# -m --model : model to test, which are found in Model_2pop and comprise 4 basic (SI, IM, AM, SC) and their extention:\n"+
    "# SC_b for SC with bottleneck SC_ae for ancestral expension SC_2N for Hill-Robertson effect SC_2M for heterogeneous gene flow and combinaton of thereof  \n"+
    "# For more information on models see paper by \n"+
    "# -s --folded : True or False: True if SFS is folded, or False if SFS unflolded using an outgroup \n"+
    "# -z : mask the singletons (no argument to provide) \n"+
    "# -b : coma-separated list of parameter obtains from the best run of previous model optimisation #" +
    "########################## Enjoy ###########################")
    return()
	

def takearg(argv):
   fs_fil_name = ''
   popname1=''
   popname2=''
   projection=''
   model=''
   folded=''
   masked=False
   best_params=''
   try:
      opts, args = getopt.getopt(argv,"hf:x:y:p:m:s:z,b:",["fs_fil_name=","popname1=","popname2=","projection=","model=","folded=","masked","best_params="])
   except getopt.GetoptError:
    print("ERROR : INCORRECT INPUT PROVIDED.\nDescription of the script usage bellow: \n\n")
    print(usage())
    sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-f", "--fs_fil_name"):
         fs_fil_name = arg
      elif opt in ("-x", "--popname1"):
         popname1 = arg
      elif opt in ("-y", "--popname2"):
         popname2 = arg
      elif opt in ("-p", "--projection"):
         projection = arg.split(',')
      elif opt in ("-m", "--model"):
         model = arg
      elif opt in ("-s", "--folded"):
         folded = arg
      elif opt in ("-z", "--masked"):
         masked = True
      elif opt in ("-b", "--best_params"):
         best_params = arg.split(',')
   print ('Input file is ', fs_fil_name)
   print ('popname1 is ', popname1)
   print ('popname2 is ', popname2)
   print ('projection of', projection,' is folded? ',folded, 'and is masked? ', masked)
   print ('Model ploted is  ', model)
   return(fs_fil_name, popname1, popname2, projection,model, folded,masked,best_params)
if __name__ == "__takearg__":
   takearg(sys.argv[1:])

###########################
###      actual script         #####
###########################

# Load parameters
fs_fil_name, popname1, popname2, projection,model, folded, masked,best_params= takearg(sys.argv[1:])


# Load the data
#Create python dictionary from snps file
input_file=fs_fil_name+'.data'
dd = moments.Misc.make_data_dict(input_file)


#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=[popname2, popname1]
print (pop_ids)


#projection sizes, in ALLELES not individuals
proj=list(map(int,projection))

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
if folded == 'True':
    data = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
    func=getattr(Model_2pop_folded, model)

if folded == 'False':
    data = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)
    func=getattr(Model_2pop_unfolded, model)

if masked:
	data.mask[1,0] = True
	data.mask[0,1] = True


#Summary statistic of the SFS
print ("\nSFS summary \n")
print ("sample sizes", data.sample_sizes)
sfs_sum = numpy.around(data.S(), 2)
print ("Sum of SFS = ", sfs_sum)
Fst = numpy.around(data.Fst(),4)
print ("FST =", Fst, '\n')



#**************
#Provide best optimized parameter set for empirical data.
#These will come from previous analyses you have already completed.
print(best_params)
emp_params =list(map(float,best_params))

#**************
#Fit the model using these parameters and return the model SFS.
#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model, but
#everything else can stay as it is. See above for argument explanations.
#================================================================================
# Fit the empirical data based on prior optimization results, obtain model SFS
#================================================================================

mod = func(emp_params, proj)
# Likelihood of the data given the model AFS.
ll_model = moments.Inference.ll_multinom(mod, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
theta = moments.Inference.optimal_sfs_scaling(mod, data)
print('Optimal value of theta: {0}'.format(theta))

# Plot a comparison of the resulting fs with the data and save it as svg
import pylab
svg_name= popname2+'_'+popname1+'_'+model +'.svg'
fig = pylab.figure(1)
moments.Plotting.plot_2d_comp_multinom(mod, data, vmin=0.1, resid_range=3,
                                    pop_ids =pop_ids)
fig.savefig(svg_name)
                                    

