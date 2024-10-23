#adapted form here: https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb
#incorporated multiprocessing and calculations from a trajectory
import numpy as np
from numpy.linalg import norm
import shutil
import wget
import subprocess
import wget
import os
import pandas as pd
import scipy.stats as scs
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from kneed import KneeLocator
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import scipy

#mol mod
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.lib.log import ProgressBar


#plot
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib_inline.backend_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'svg')

#custom libs
import BME as BME
from main import BlockAnalysis


######################DOWNLOAD files. Only remove global comment if needd to download##########
'''
#rm -r sample_data

#wget https://cssb.biology.gatech.edu/skolnick/files/PULCHRA/pulchra304.tgz &> /dev/null
#tar -zxf pulchra304.tgz &> /dev/null
#rm pulchra304.tgz
#mv ./pulchra304/bin/linux/pulchra .
#chmod +x pulchra
#rm -r pulchra304/
wget https://github.com/fpesceKU/EnsembleLab/raw/main/utils/pulchra &> /dev/null
chmod +x pulchra'''

'''wget https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/Linux/3.0/Pepsi-SAXS-Linux.zip
unzip Pepsi-SAXS-Linux.zip &> /dev/null
rm Pepsi-SAXS-Linux.zip

wget.download('https://raw.githubusercontent.com/fpesceKU/BLOCKING/main/block_tools.py') 
wget.download('https://raw.githubusercontent.com/fpesceKU/BLOCKING/main/main.py')

wget.download('https://raw.githubusercontent.com/KULL-Centre/BME/main/BME_tools.py')
wget.download('https://raw.githubusercontent.com/KULL-Centre/BME/main/BME.py')

#wget https://github.com/ehb54/GenApp-BayesApp/raw/main/bin/source/old_versions/bift_v33.f
# &> /dev/null
#gfortran bift_v33.f -march=native -O2 -o bift'''
################################# End block ################

def cm_dist(u,sel1,sel2):
    cm1 = u.select_atoms(sel1).center_of_mass()
    cm2 = u.select_atoms(sel2).center_of_mass()
    return norm(cm1 - cm2)


def dmax (u,sel):
    CA_selection = u.select_atoms(sel)
    CA_coordinates = CA_selection.positions #expose numpy array of coords
    distance_matrix_pool = scipy.spatial.distance.cdist(CA_coordinates, CA_coordinates)
    maximum_distance_pool = distance_matrix_pool.max()
    return norm(maximum_distance_pool)

def eed(u):
    nterm = u.select_atoms('protein and name N')[0]
    cterm = u.select_atoms('protein and name C')[-1]
    return norm(cterm.position - nterm.position)


def calc_traj(t):
    #rgarray = md.compute_rg(t)
    results = []
    
    #last_frame=u_eom.trajectory.n_frames-1
    
    for ts in ProgressBar(t.trajectory[:]):
        results.append({'frame':ts.frame+1, 'rg':t.atoms.radius_of_gyration(), 'dmax':dmax(t,"protein"), 'cm_dist':cm_dist(t,"resid 1-10","resid 20-30"), 'eed':eed(t)})

    df=pd.DataFrame.from_dict(results)
    return df


def kde(a, w=None, min_=None, max_=None):
    if type(w) == 'NoneType':
        w = np.full(len(a), 1)
    if min_ == None:
        min_ = np.min(a)
    if max_ == None:
        max_ = np.max(a)
    x = np.linspace( min_, max_, num = 50)
    d = scs.gaussian_kde( a, bw_method = "silverman", weights = w ).evaluate(x)
    u = np.average(a, weights = w)
    return x,d/np.sum(d),u
    

def plot_dist(ax,x,p,av):
    ax.plot(x,p,c='k')
    ax.vlines(av,0,100, color='k')
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0,np.max(p)+0.8*np.max(p))

def plot_rew_dist(ax,x,p,av):
    ax.plot(x,p,c='tab:red',ls='dashed')
    ax.vlines(av,0,100, color='tab:red', linestyle='dashed')
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0,np.max(p)+0.8*np.max(p))


def autoblock(cv, multi=1, plot=False):
    block = BlockAnalysis(cv, multi=multi)
    block.SEM()

    if plot == True:
        plt.errorbar(block.stat[...,0], block.stat[...,1], block.stat[...,2], fmt='', color='k', ecolor='0.5')
        plt.scatter(block.bs, block.sem,zorder=10,c='tab:red')
        plt.xlabel('Block size')
        plt.ylabel('SEM')
        plt.show()

    return block.av, block.sem, block.bs

def traj2saxs(frame_index):
    frame_file = f'tmp_frame_{frame_index}.pdb'
    output_file = f'tmp_saxs_{frame_index}.dat'

    traj.trajectory[frame_index]  # zero indexed. It updates the traj and outputs the current state no neeed to assign variable

    # Save the frame to a new file
    with mda.Writer(frame_file, traj.atoms.n_atoms) as W:
        W.write(traj.atoms)
    

    #flag --dist calculates p(r)
    pepsi_comm = f'/Users/tiago/software_science/pepsisaxs/Pepsi-SAXS {frame_file} {NAME}_bift.dat -o {output_file} -cst -x'
    #pepsi_comm = f'/Users/tiago/software_science/pepsisaxs/Pepsi-SAXS {frame_file} {NAME}_bift.dat -o {output_file} -cst -x --cstFactor 0 --I0 1.0 --dro 1.0 --r0_min_factor 1.025 --r0_max_factor 1.025 --r0_N 1'
    subprocess.run(pepsi_comm.split(), stdout=subprocess.DEVNULL)
    
    calc_saxs = np.loadtxt(output_file)[..., 3]
    
    os.remove(frame_file)  # Clean up temporary files
    os.remove(output_file)
    
    return calc_saxs

'''#using mdanalysis
def traj2saxs(frame_index):
    traj.trajectory[frame_index]  # This sets the universe to the desired frame
    frame_file = f'tmp_frame_{frame_index}.pdb'
    output_file = f'tmp_saxs_{frame_index}.dat'
    
    # Save the current frame to a PDB file
    traj.atoms.write(frame_file)
    
    # Prepare and run the command for Pepsi-SAXS
    pepsi_comm = f'/Users/tiago/software_science/pepsisaxs/Pepsi-SAXS {frame_file} {NAME}_bift.dat -o {output_file} -cst -x'
    subprocess.run(pepsi_comm.split(), stdout=subprocess.DEVNULL)
    
    # Load the calculated SAXS data
    calc_saxs = np.loadtxt(output_file)[..., 3]
    
    # Clean up temporary files
    os.remove(frame_file)
    os.remove(output_file)
    
    return calc_saxs'''

def iBME(theta, exp_file, calc_file):
    print('Reweighting with theta={}'.format(theta))
    rew = BME.Reweight('ibme_t{}'.format(theta))
    rew.load(exp_file, calc_file)
    rew.ibme(theta=theta, iterations=25, ftol=0.001)
    weights = rew.get_ibme_weights()[-1]
    stats = rew.get_ibme_stats()[-1]
    #print('chi2={:.2f}, phi_eff={:.2f}'.format(stats[-1][1],stats[-1][2]))
    return theta, stats, weights


# Function to find the optimal theta using the KneeLocator
def find_optimal_theta(thetas, stats):
    knee = KneeLocator(stats[:, 2], stats[:, 1], S=1, curve="convex", direction="increasing")
    if knee.knee is None:
        return None  # or choose a default value or handle this case as needed
    optimal_indices = np.where(stats[:, 2] == knee.knee)[0]
    if len(optimal_indices) > 0:
        optimal_index = optimal_indices[0]
        choice = thetas[optimal_index]
        return choice
    else:
        return None  # or choose a default value or handle this case as needed

############## LOAD traj and experimental SAXS #################

NAME = "Ubq4" 

#mdtraj
#traj = md.load_xtc('./traj/{:s}_trajectory.xtc'.format(NAME), top='./traj/top_{:s}.pdb'.format(NAME))

#mdanalysis
traj = mda.Universe('./traj/top_{:s}.pdb'.format(NAME), './traj/{:s}_trajectory.xtc'.format(NAME))

saxs_file = '/Users/tiago/Desktop/Students/Nuno_Fernandes/bme_tests/exp_saxs/bift_{:s}.dat'.format(NAME)

frame_indices = range(len(traj.trajectory))


num_workers = 12  # Adjust based on your system's resources


if __name__ == '__main__':
############################### Prepare experimental SAXS profile ###############

    EXPERIMENT = "SAXS" 

    if EXPERIMENT == "SAXS":
        #print('SAXS data must be in a file containing 3 columns, which are q, I and sigma. Commented lines (#) are allowed.')

        #check data
        try:
            np.loadtxt(saxs_file)
        except:
            print("Unable to read file. Make sure the file only contains 3 columns (q,I,sigma) and #commented lines")
        assert np.shape(np.loadtxt(saxs_file))[1] == 3, "Expected file with 3 columns (q,I,sigma)"

        exp_saxs = np.loadtxt(saxs_file)
        if exp_saxs[...,0][-1] <  1:
            print('q is in Å units. Converting to nm.')
            exp_saxs[...,0] = exp_saxs[...,0]*10
            np.savetxt(saxs_file, exp_saxs)

        #To cut the SAXS profile change value here (e.g. noise or aggregation)
        if (exp_saxs[...,0] >= 5).sum() > 0:
            print('Found {} q-values above 5 nm^(-1). SAXS calculations are not reliable in that region of the spectrum. Those datapoints will be removed'.format((exp_saxs[...,0] >= 5).sum()))
            exp_saxs = exp_saxs[(exp_saxs[...,0] < 5)]
            np.savetxt(saxs_file, exp_saxs)

        shutil.copy(saxs_file, 'saxs_input.dat')


    #correct experimental errors with BIFT

    #f = open('inputfile.dat','w')
    #f.write("saxs_input.dat\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
    #f.close()
    #subprocess.run('./bift < inputfile.dat'.split())
    np.savetxt('{:s}_bift.dat'.format(NAME), np.loadtxt('saxs_input.dat'), header=' DATA=SAXS') #quick_fix to bypass bift. Experimental SAXS profile
    #print('Experimental errors on SAXS intensities have been corrected with BIFT')
    #print('Factor used for rescaling errors is: {}'.format(np.loadtxt('scale_factor.dat')[0,1]))
    #print('SAXS data with corrected errors is in {:s}_bift.dat\n'.format(NAME))


    ################################## Back-calculate from trajectory #############################
    print('Calculating metrics from trajectory...')

    #Rg
    #metrics.rg = Rg(traj)
    #print(metrics.rg)
    #rg_av, rg_err, rg_bs = autoblock(metrics.rg, plot=False) #block calc if ensemble from MD sim

    metrics = calc_traj(traj)

    rg_av = np.average(metrics.rg)
    x_rg, p_rg, _ = kde(metrics.rg)


    #SAXS
    print('Calculating SAXS from trajectory...')

    #multiprocessor
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        all_calc_saxs = list(tqdm(executor.map(traj2saxs, frame_indices), total=len(frame_indices)))


    # Save all calculated SAXS data
    all_calc_saxs = np.vstack(all_calc_saxs)
    #all_calc_saxs.shape()
    col0 = np.arange(1, len(all_calc_saxs) + 1).reshape(len(all_calc_saxs), 1)
    all_calc_saxs = np.hstack((col0, all_calc_saxs))
    #all_calc_saxs.shape()

    np.savetxt('calc_saxs.dat', all_calc_saxs)

    ################################ END of back-calculate ###############################


    ################### Maxent optimization block ####################
    
    thetas = np.array([1, 10, 20, 50, 75, 100, 200, 400, 750, 1000, 5000, 10000])
    exp_file = '{:s}_bift.dat'.format(NAME)
    calc_file = 'calc_saxs.dat'

    # Parallelize the computation across the theta values
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(iBME, theta, exp_file, calc_file) for theta in thetas]
        results = [f.result() for f in futures]

    # Extract the results
    results.sort()  # Ensure the results are sorted by theta
    thetas, stats, weights = zip(*results)

    # Convert to numpy arrays for easier indexing
    stats = np.array(stats)

    # Find the optimal theta
    choice = find_optimal_theta(thetas, stats)
    print(choice)

    # Plotting the results
    mpl.rcParams.update({'font.size': 10})
    fig = plt.figure(layout='constrained', dpi=300, figsize=(8, 3))
    plt.scatter(stats[..., 2], stats[..., 1], c='k')
    ndx = np.where(thetas == choice)[0][0]
    plt.scatter(stats[..., 2][ndx], stats[..., 1][ndx], c='tab:red', label=r'Chosen $\theta$')
    plt.xlabel(r'$\phi_{eff}$', fontsize=10)
    plt.ylabel(r'$\chi^2_r$', fontsize=10)
    plt.legend()
    plt.title(NAME)
    plt.show()


    ###############Plot after and before reweight. Red is reweighted. Blue is unreweighted ##############
    #weights_z = np.sum(weights[ndx]) // np.min(weights[ndx])



    print(f'Reweighting using theta = {choice}')

    #Linear fitting of scale and offset from SAXS profiles
    q, I_exp, err = np.loadtxt('{}_bift.dat'.format(NAME), unpack=True)
    I_prior = np.average(np.loadtxt('calc_saxs.dat')[...,1:],axis=0)
    I_post = np.average(np.loadtxt(list(filter(lambda x: x.startswith('ibme_t{}_'.format(choice)) and x.endswith('.calc.dat'), os.listdir('.')))[0])[...,1:], axis=0, weights=weights[ndx])
    wlr = 1/(err**2)
    model = LinearRegression()
    model.fit(I_prior.reshape(-1,1),I_exp,wlr)
    a = model.coef_[0]
    b = model.intercept_
    I_prior = a*I_prior+b


    mpl.rcParams.update({'font.size': 14})
    fig, axs = plt.subplots(2, 2, figsize=(12,8), facecolor='w', dpi=300, layout='constrained')
    axs = axs.flatten()


    #Exp SAXS vs SAXS reweighted (I_post) and unreweighted. No log scale
    axs[0].errorbar(q,I_exp,err, lw=1,c='0.5',alpha=0.5)
    axs[0].plot(q,I_prior, lw=1, zorder=500)
    axs[0].plot(q,I_post, lw=1, color='tab:red', ls='dashed', zorder=1000)
    axs[0].set_xlabel(r'q [nm$^{-1}$]')
    axs[0].set_ylabel('Intensity')


    #Exp SAXS vs SAXS reweighted and unreweighted. yy axis log scale
    axs[1].errorbar(q,I_exp,err, lw=1, c='0.5',alpha=0.5)
    axs[1].plot(q, I_exp, lw=3, c='grey')  # plot exp in higher lw
    axs[1].plot(q,I_prior, lw=2, zorder=500)
    axs[1].plot(q,I_post, lw=2, color='tab:red', ls='dashed', zorder=1000)
    #axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xlabel(r'q [nm$^{-1}$]')
    axs[1].set_ylabel('Intensity')


    #Kratky plots
    kratky_exp = (q**2)*I_exp
    kratky_err = (q**2)*err
    axs[2].errorbar(q,kratky_exp,kratky_err, lw=1,c='0.5',alpha=0.5)
    axs[2].plot(q,kratky_exp, lw=3,c='grey',alpha=0.5)
    axs[2].plot(q,(q**2)*I_prior, lw=2, zorder=500)
    axs[2].plot(q,(q**2)*I_post, lw=2, color='tab:red', ls='dashed', zorder=1000)
    axs[2].set_xlabel(r'q [nm$^{-1}$]')
    axs[2].set_ylabel(r'q$^2$I')
    axs[2].set_ylim(0,np.max((q**2)*I_prior)+ 0.1*np.max((q**2)*I_prior))

    #Residuals
    axs[3].plot(q, (I_exp-I_prior)/err, color='tab:blue', lw=1 )
    axs[3].plot(q, (I_exp-I_post)/err, color='tab:red', lw=1, alpha=0.6 )
    axs[3].set_xlabel(r'q [nm$^{-1}$]')
    axs[3].set_ylabel(r'(I$^{EXP}$-I$^{CALC}$)/$\sigma$')
    

    ##########################################################################################
    #            Plot and calculate distribuitions before and after reweight                 #                                                          #
    ########################################################################################## 

    mpl.rcParams.update({'font.size': 10})
    fig, axs = plt.subplots(2, 2, figsize=(10,6), facecolor='w', dpi=300, layout='constrained')
    axs = axs.flatten()


    #Rg
    x_rg_rew, p_rg_rew, rg_av_rew = kde(metrics.rg, w=weights[ndx], min_=np.min(metrics.rg), max_=np.max(metrics.rg))
    rg_av_rew = np.average(metrics.rg, weights=weights[ndx])


    plot_dist(axs[0], x_rg, p_rg, rg_av)
    plot_rew_dist(axs[0], x_rg_rew, p_rg_rew, rg_av_rew)
    

    axs[0].set_xlabel(r'$R_g$ [nm]')
    axs[0].set_ylabel(r'p($R_g$)')
    axs[0].text(0.75, 0.9, r'$\langle R_g \rangle$ reweight ={:.2f} nm'.format(rg_av_rew), color='tab:red', horizontalalignment='center',verticalalignment='center', transform=axs[0].transAxes, fontsize=10)
    axs[0].text(0.75, 0.8, r'$\langle R_g \rangle$ uniform ={:.2f} nm'.format(rg_av), color='black', horizontalalignment='center',verticalalignment='center', transform=axs[0].transAxes, fontsize=10)

    #eed
    eed_av = np.average(metrics.eed)
    x_eed, p_eed, _ = kde(metrics.eed)

    x_eed_rew, p_eed_rew, eed_av_rew = kde(metrics.eed, w=weights[ndx], min_=np.min(metrics.eed), max_=np.max(metrics.eed))
    eed_av_rew = np.average(metrics.eed, weights=weights[ndx])

    plot_dist(axs[1], x_eed, p_eed, eed_av)
    plot_rew_dist(axs[1], x_eed_rew, p_eed_rew, eed_av_rew)        

    axs[1].set_xlabel(r'$EED$ [nm]')
    axs[1].set_ylabel(r'p($EED$)')
    axs[1].text(0.75, 0.9, r'$\langle EED \rangle$ reweight ={:.2f} nm'.format(eed_av_rew), color='tab:red', horizontalalignment='center',verticalalignment='center', transform=axs[1].transAxes, fontsize=10)
    axs[1].text(0.75, 0.8, r'$\langle EED \rangle$ uniform ={:.2f} nm'.format(eed_av), color='black', horizontalalignment='center',verticalalignment='center', transform=axs[1].transAxes, fontsize=10)

    #dmax
    dmax_av = np.average(metrics.dmax)
    x_dmax, p_dmax, _ = kde(metrics.dmax)

    x_dmax_rew, p_dmax_rew, dmax_av_rew = kde(metrics.dmax, w=weights[ndx], min_=np.min(metrics.dmax), max_=np.max(metrics.dmax))
    dmax_av_rew = np.average(metrics.dmax, weights=weights[ndx])

    plot_dist(axs[2], x_dmax, p_dmax, dmax_av)
    plot_rew_dist(axs[2], x_dmax_rew, p_dmax_rew, dmax_av_rew)        

    axs[2].set_xlabel(r'$dmax$ [nm]')
    axs[2].set_ylabel(r'p($dmax$)')
    axs[2].text(0.75, 0.9, r'$\langle dmax \rangle$ reweight ={:.2f} nm'.format(dmax_av_rew), color='tab:red', horizontalalignment='center',verticalalignment='center', transform=axs[2].transAxes, fontsize=10)
    axs[2].text(0.75, 0.8, r'$\langle dmax \rangle$ uniform ={:.2f} nm'.format(dmax_av), color='black', horizontalalignment='center',verticalalignment='center', transform=axs[2].transAxes, fontsize=10)


    #cm distance
    cm_dist_av = np.average(metrics.cm_dist)
    x_cm_dist, p_cm_dist, _ = kde(metrics.cm_dist)

    x_cm_dist_rew, p_cm_dist_rew, cm_dist_av_rew = kde(metrics.cm_dist, w=weights[ndx], min_=np.min(metrics.cm_dist), max_=np.max(metrics.cm_dist))
    cm_dist_av_rew = np.average(metrics.cm_dist, weights=weights[ndx])

    plot_dist(axs[3], x_cm_dist, p_cm_dist, cm_dist_av)
    plot_rew_dist(axs[3], x_cm_dist_rew, p_cm_dist_rew, cm_dist_av_rew)        

    axs[3].set_xlabel(r'$cm_dist$ [nm]')
    axs[3].set_ylabel(r'p($cm_dist$)')
    axs[3].text(0.75, 0.9, r'$\langle cm_dist \rangle$ reweight ={:.2f} nm'.format(cm_dist_av_rew), color='tab:red', horizontalalignment='center',verticalalignment='center', transform=axs[3].transAxes, fontsize=10)
    axs[3].text(0.75, 0.8, r'$\langle cm_dist \rangle$ uniform ={:.2f} nm'.format(cm_dist_av), color='black', horizontalalignment='center',verticalalignment='center', transform=axs[3].transAxes, fontsize=10)


    '''sns.kdeplot(metrics.rg, weights=weights[ndx], ax=axs[2])
    sns.kdeplot(metrics.rg, ax=axs[2])'''