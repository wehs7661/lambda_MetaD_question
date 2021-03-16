import os 
import time 
import glob
import copy
import yaml
import argparse 
import numpy as np 
import matplotlib.pyplot as plt 
from tqdm.auto import tqdm
from collections import OrderedDict
from matplotlib import rc 


def represent_dictionary_order(self, dict_data):
    """
    Works with the setup_yaml() are to preserve the order of the dictionary 
    to be printed to the trajectory yaml file.
    """
    return self.represent_mapping('tag:yaml.org,2002:map', dict_data.items())

def setup_yaml():
    """
    Setup yaml files
    """
    yaml.add_representer(OrderedDict, represent_dictionary_order)


def initialize():
    """
    This function is an argument parser to read the user input arguments. 
    """
    parser = argparse.ArgumentParser(
        description='This code calculates the free energy profile and assess \
                    the influence of the block size on its uncertainty.')
    parser.add_argument('-i',
                        '--input',
                        default='plumed.dat',
                        help='The PLUMED input file for running metadynamics. Default: plumed.dat')
    parser.add_argument('-ir',
                        '--input_re',
                        default='plumed_reweight.dat',
                        help='The PLUMED input file for reweighting. Default: plumed_reweight.dat')
    parser.add_argument('-ih',
                        '--hills',
                        default=glob.glob('HILLS*')[0],
                        help='The filename of the HILLS file. Default: HILLS*')
    parser.add_argument('-b',
                        '--b_sizes',
                        required=True,
                        type=int,
                        nargs='+',
                        help='The min, max and the spacing of the block sizes, in the units of \
                              the STRIDE used to print COLVAR in metadynamics. For example, if \
                              the this argument is set to 5000 and STRIDE=10 was used in the PLUMED\
                              input file for running the metadynamics, then the real block size is \
                              50000 simulation steps. If dt = 0.002 ps, this means 100 ps. For this\
                              argument, if only one value is given, only one free energy profile based \
                              on the given value will be calcualted. (Only one fes_*.dat file.)')
    parser.add_argument('-p',
                        '--periods',
                        type=float,
                        required=True,
                        nargs='+',
                        help='The periods of the CV values. For example, if a dihedral angle and \
                              a coordination number are biased in a metadynamics, period = [2 * np.pi, 0]. Default: 0.')
    parser.add_argument('-dt',
                        '--dt',
                        type=float,
                        default=0.002,
                        help='The simulation time step (ps). Default: 0.002 ps.')
    parser.add_argument('-t',
                        '--truncation',
                        type=float,
                        default=0,
                        help='The truncation fraction. Default: 0.')
    args_parse = parser.parse_args()

    return args_parse 

def logger(file_name='free_energy_results.txt', *args, **kwargs):
    print(*args, **kwargs)
    with open(file_name, "a") as f:
        print(file=f, *args, **kwargs)

def clear_directory(file_name):
    if type(file_name) == str:
        if os.path.isfile(file_name):
            os.remove(file_name)
    elif type(file_name) == list:
        for i in file_name:
            clear_directory(i)

def prepare_plumed_input(template, b_size):
    f = open(template, 'r')
    lines = f.readlines()
    f.close()
    
    line_n = -1  # line number
    for line in lines:
        line_n += 1
        if 'CLEAR' in line:
            lines[line_n] = f'CLEAR={b_size}\n'
        if 'DUMPGRID' in line:
            line_split = line.split('STRIDE=')
            line_split[-1] = f'STRIDE={b_size}\n'
            lines[line_n] = ''.join(line_split)

    output = open(template, 'w')
    output.write(''.join(lines))       
    output.close()

def read_histogram(hist_file):
    data = np.loadtxt(hist_file)
    with open(hist_file, "r") as f:
        for line in f:
            if line.startswith('#! SET normalisation'):
                norm = float(line.split()[-1])
    hist = data[:, -1]  # normalized probability: the last column
    
    return norm, hist


def calculate_free_energy(hist_dir, hist_files):
    # Step 1: Calculate the average of the weighted histogram for each gridded CV value
    w_sum, avg = 0, 0
    for f in hist_files:
        norm, hist = read_histogram(f'{hist_dir}/{f}')
        w_sum += norm
        avg += norm * hist
    avg = avg / w_sum
    
    # Step 2: Calculate the uncertainty of each gridded CV value
    error = 0
    for f in hist_files:
        norm, hist = read = read_histogram(f'{hist_dir}/{f}')
        error += norm * norm * (hist - avg) ** 2
    error = np.sqrt(error / (w_sum ** 2))
        
    # Step 3: Conver to the uncertainty in free energy 
    fes = -np.log(avg)    # units: kT
    # f_err = error / avg   # units: kT

    return fes

class SimulationParmeters:
    """
    A class for getting the simulation parameters of well-tempered metadynamics.
    """
    def __init__(self, plumed_input, dt):
        """
        This function gets the simulation parameters of well-tempered metadynamics.
        
        Parameters
        ----------
        plumed_input : str 
            The filename of the plumed input file for running metadynamics.
        dt : float
            The simulation time step in ps. 

        Attributes
        ----------
        dt : float
            The simulation time step (units: ps).
        dt_c : float
            The factor converting the block size to ps. / The time step in the COLVAR file in ps.
        dt_g : float 
            The Gaussian deposition stride in ps. / The time step in the HILLS file in ps. 
        stride_c : float 
            The stride for writing data to COLVAR. (units: steps)
        stride_g : float 
            The stride for depositing Gaussians (units: steps).
        T : float
            The simulation temperature (units: K).
        gamma : float 
            The bias factor used for well-tempered metadynamics.
        delta_T : float 
            The \Delta T parameter in well-tempered metadynamics.
        k_delta_T : float
            k * \Delta T in kJ/mol.
        """
        f = open(plumed_input, 'r')
        lines = f.readlines()
        f.close()

        NA = 6.02214076E23   # Avogadro constant
        kB = 1.38064852E-23  # Boltzmann constant
            
        for l in lines:
            if 'PACE' in l:
                self.stride_g = float(l.split('#')[0].split('=')[1])  # stride for depositing gaussian (units: steps)
            if 'TEMP' in l:
                self.T = float(l.split('#')[0].split('=')[1])
            if 'BIASFACTOR' in l:
                self.gamma = float(l.split('#')[0].split('=')[1])
            if 'COLVAR' in l:
                for i in l.split():
                    if 'STRIDE=' in i:  # stride for writing data to COLVAR
                        self.stride_c = float(i.split('=')[1])
        self.dt = dt 
        self.dt_c = self.stride_c * self.dt   
        self.dt_g = self.stride_g * self.dt
        self.delta_T = (self.gamma - 1) * self.T 
        self.k_delta_T = kB * self.delta_T * NA / 1000  # k\Delta T in kJ/mol

class UncertaintyEstimator(SimulationParmeters):
    """
    A class for estimating the uncertainty of the free energy surface.
    """
    def __init__(self, plumed_input, hills, dt):
        """
        This function gets the data of the HILLS file for later calculations.
        
        Parameters
        ----------
        plumed_input : str 
            The filename of the plumed input file for running metadynamics.
        hills : str
            The filenname of the HILLS file output from metadynamics
        dt : float
            The simulation time step in ps. 

        Attributes
        ----------
        n_CVs : int 
            The number of CVs. 
        n_gaussians : int
            The number of Gaussians added in the simulation
        length : float
            The simulation length in ps.
        data: numpy.array
            The content of the HILLS file.
        mu: numpy.array
            The time series of the CV values visited when depositing Gaussians in the form of array([[CV_1], [CV_2], ...]). / Centers of the Gaussians. 
        sigma : numpy.array 
            The time series of the Gaussian width in the form of array([[sigma_1], [sigma_2], ...]).
        h : numpy.array
            The time series of the Gaussian height (units: kJ/mol).
        """
        SimulationParmeters.__init__(self, plumed_input, dt)
        self.data = np.transpose(np.loadtxt(hills))
        self.n_CVs = int((len(self.data) - 3) / 2)
        self.mu = self.data[1 : self.n_CVs + 1]
        self.sigma = self.data[self.n_CVs + 1 : 2 * self.n_CVs + 1]
        self.h = self.data[-2] * ((self.gamma - 1) / self.gamma)
        self.n_gaussians = len(np.transpose(self.data))
        self.length = self.n_gaussians * self.dt_g  # simulation length (ps)

    def get_uni_gaussian(self, x, mu, sigma, period=0, h=1):
        """
        This function calculates the value of a univarite Gaussian at a given CV value.

        Parameters
        ----------
        x : float
            The value of a certain CV (one-dimensional).
        mu : float
            The center of the Gaussian distribution.
        sigma : float
            The width of the Gaussian distribution.
        period : float
            The period of the CV, if any. For example, if the CV is a dihedral angel, 
            then period=2 * np.pi. Default: 0.
        h : float
            The height of the Gaussian distribution. Default: 1.
        
        Returns
        -------
        prob : float
            The probability at the given CV value. 
        """
        delta_x = np.min([(period - np.abs(x - mu)), np.abs(x - mu)])
        prob = h * np.exp(- delta_x ** 2 / (2 * sigma ** 2))

        return prob

    def get_multi_gaussian(self, x, mu, sigma, period, h):
        """
        This function calculates the value of a multi-variate Gaussian at 
        a given set of CV values. The validity of this function can be tested by 
        the last column of a COLVAR file.

        Parameters
        ----------
        x : array-like
            The set of CV values.
        mu : array-like
            The centers of the Gaussians in different dimensions. For example, 
            [0.5, 02] can be used to specify the width a 2D Gaussian. Note the 
            the multivare Gaussians in metadynamics are composed of independent 
            univariate Guassians.
        period : array-like
            The periods of the CV values. For example, if a dihedral angle and 
            a coordination number are biased in a metadynamics. period = [2 * np.pi, 0].
        height : float
            The height of the multi-variate Gaussian.

        Returns
        -------
        prob : float
            The probability at the given set of CV values.
        """
        prob = 1  
        for i in range(self.n_CVs):
            p = self.get_uni_gaussian(x[i], mu[i], sigma[i], period[i])
            prob *= p
        prob *= h

        return prob

    def get_accumulated_bias(self, x, period, t):
        """
        This function calculates the bias accumulated to time t (in ps), but not including
        the bias added at time t, if any. The validity of this function can be tested by 
        the last column of a COLVAR file.
        
        Parameters
        ----------
        x : array-like
            The set of CV values.
        period : array-like
            The periods of the CV values. For example, if a dihedral angle and 
            a coordination number are biased in a metadynamics. period = [2 * np.pi, 0].
        t : float
            The value of time (in ps).

        Returns
        -------
        bias : float
            The accumulated bias.
        """
        bias = 0

        # First get the index of t of interest. Note that we need to consider the deposition stride.
        t_index = int(np.floor(t / self.dt_g) - 1)   # np.floor(t / self.dt_g) is the number of Gaussian deposited before t

        for i in tqdm(range(t_index + 1)):  # tqdm: progress bar
            mu = np.transpose(self.mu)[i]         # multi-dimensional mu at time t
            sigma = np.transpose(self.sigma)[i]   # multi-dimensional sigma at time t
            bias += self.get_multi_gaussian(x, mu, sigma, period, self.h[i])
        
        return bias

    def estimate_error(self, x, period, t_fill, L_b, N_b):
        """
        This function estimates the uncertainty of the free energy at a given 
        set of CV values x. This is the implementation of Equation 10 in the paper
        "Using metadynamics to explore complex free-energy landscapes" by Dr. Bussi.

        Parameters
        ----------
        x : array-like
            The set of CV values of interest.
        period : array-like
            The periods of the CV values. For example, if a dihedral angle and 
            a coordination number are biased in a metadynamics. period = [2 * np.pi, 0].
        t_fill : float
            The simulation length to be truncated (units: ps).
        L_b : float
            Block size (units: ps). 
        N_b : int
            The number of blocks.
        
        Returns
        -------
        f_err : float
            The uncertainty of the free energy at a given set of CV values.
        """
        log_term = []
        for i in tqdm(range(N_b)):  # tqdm: progress bar
            b1 = self.get_accumulated_bias(x, period, t_fill + (i + 1) * L_b)  # kJ/mol
            b2 = self.get_accumulated_bias(x, period, t_fill + i * L_b)        # kJ/mol
            if b1 - b2 != 0:
                # it is possible that b1 == b2 when the system does not sample during one block
                exp1 = np.exp(b1 / self.k_delta_T)   
                exp2 = np.exp(b2 / self.k_delta_T)
                log_term.append(np.log(exp1 - exp2))

        f_err = np.var(log_term) * (1 / np.sqrt(N_b))

        return f_err


if __name__ == '__main__':
    t1 = time.time()
    # Setting up parameters and initializing relevant classes
    args = initialize()
    setup_yaml()
    
    if len(args.b_sizes) == 1:
        multi_b_size = False
        blocks = [args.b_sizes[0]]
    else:  
        multi_b_size = True
        try:
            os.system('rm -r FES')
        except:
            pass
        os.makedirs('FES')
        blocks = np.arange(args.b_sizes[0], args.b_sizes[1] + args.b_sizes[2], args.b_sizes[2])    

    SP = SimulationParmeters(args.input, args.dt)
    UE = UncertaintyEstimator(args.input, args.hills, args.dt)
    t_fill = args.truncation * UE.length  # ps

    # Document the parameters for the data analysis
    with open('parameters.yaml', 'w') as params:
        params.write('Section 1: User inputs\n')
        params_args = copy.copy(vars(args))
        params_args.pop('b_sizes')
        params.write(f'b_sizes: {args.b_sizes}\n')
        yaml.dump(params_args, params, default_flow_style=False)    
        
        params.write('\nSection 2: Simulation parameters\n')
        yaml.dump(vars(SP), params, default_flow_style=False)
        params_ue = copy.copy(vars(UE))
        for i in ['data', 'mu', 'sigma', 'h']:
            params_ue.pop(i)
        yaml.dump(params_ue, params, default_flow_style=False)

    delta_f_err = []
    for  b_size in blocks:
        # Step 1: Clear directory
        rm_list = [
            'COLVAR_REWEIGHT', 
            'weights.dat', 
            glob.glob('bck.*'), 
            glob.glob('fes*'), 
            glob.glob('histograms/*')
        ]
        for i in rm_list:
            clear_directory(i)
        
        # Step 2: Prepare the PLUMED input file
        prepare_plumed_input(args.input_re, b_size)

        # Step 3: Run plumed driver to get weighted histograms
        os.system(f'plumed driver --plumed {args.input_re} --noatoms')
        
        # Step 4: Free energy calculations 
        """
        files = []
        for f in os.listdir('histograms'):
            if f.endswith('hist.dat'):
                files.append(f)
        files = natsort.natsorted(files, reverse=False)
        N_b = len(files)    # the number of blocks
        fes = calculate_free_energy('histograms', files)
        """

        files = glob.glob('histograms/*hist.dat')
        N_b = len(files)    # the number of blocks
        fes = calculate_free_energy('.', files)

        CV_points = []
        
        for i in range(UE.n_CVs):
            CV_points.append(np.transpose(np.loadtxt('histograms/hist.dat'))[i])  # time series of the i-th CV
        
        if multi_b_size is False:   # only one fes_*.dat
            fes_file = open(f'fes_bsize_{b_size}.dat', 'w')
        else:
            fes_file = open(f'FES/fes_bsize_{b_size}.dat', 'w')

        f_err = []
        for i in tqdm(range(len(CV_points[0]))):  # tqdm: progress bar
            # calculate the uncertainty for each CV point
            f_err.append(UE.estimate_error(np.transpose(CV_points)[i], args.periods, t_fill, b_size * UE.dt_c, N_b))

            CV_str = ''
            for j in range(UE.n_CVs):
                CV_str += f'    {CV_points[j][i]: .3f}'  # the value of the j-th dimension of 
            fes_file.write(f'{CV_str}    {fes[i]: .6f}    {f_err[i]: .6f}\n')
            
        delta_f_err.append(np.sqrt(f_err[-1] ** 2 + f_err[0] ** 2))
        if multi_b_size:
            logger('FES/free_energy_results.txt', f'The free energy difference is {fes[-1] - fes[0]: .6f} +/- {delta_f_err[-1]: .6f} kT. (Block size: {b_size} steps, or {b_size * UE.dt_c} ps.)')
        else:
            print(f'The free energy difference is {fes[-1]-fes[0]: .6f} +/- {delta_f_err[-1]: .6f} kT. (Block size: {b_size} steps, or {b_size * UE.dt_c} ps.)')

    # Plotting the results if needed
    rc('font', **{
       'family': 'sans-serif',
       'sans-serif': ['DejaVu Sans'],
       'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    if multi_b_size:
        plt.figure()
        plt.plot(blocks * args.factor, delta_f_err)
        plt.xlabel('Block size (ps)')
        plt.ylabel('Uncertainty in the free energy difference ($k_{B}T$)')
        plt.title('The uncertainty as a function of block size')
        plt.grid()
        plt.savefig('FES/delta_f_blocks.png', dpi=600)

    t2 = time.time()
    print(f'Time elapsed: {t2-t1: .2f} seconds.')
    



