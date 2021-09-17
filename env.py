"""Environment class simulates a volume of air in which temperature, pressure
and humidity are constant and homogeneous. Its elements take these parameters
as inputs, to determine a set of properties for the gas, in particular the
speed of sound.
"""

# from pathlib import Path
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.misc as misc
from matplotlib import pyplot as plt
import seaborn as sns

TEMP_0 = 273.15 #K
ATM_P = 101325 #Pa
CP1_STP = 1.0057E3 #CP0 = 1.006E3 #j Kg^-1 K^-1
R = 8.31446261815324 #m^3 Pa K^-1 mol^-1

#fitting functions

#Does not converge for Dry Air, to fix
def lennard (temp: float or np.ndarray, par_m: int, par_e: float, par_s: float):
    """Ab initio definition of second virial potential B, as the integral
    for r€[0,inf[ of (1-exp(U(r)/KbT)); where U is the interatomic
    potential. Valid description for Dry Air without CO2.
    Taken from Sengers, Klein and Gallagher, (1971)
    'Pressure-volume-temperature relationships of
    gases-virial coefficients.'

    The interatomic potential is treated as a (m-6) potential,
    which for m=12 becomes the Lennard-Jones Potential

    Parameters
    ----------
    temp : Float / Numpy.ndarray of floats
        Temperature of the gas, expressed in K.
    par_m : Int > 0
        m parameter of the interatomic potential, dimensionless quantity;
        determines the entity of repulsive potential and the amplitude.
    par_e : Float
        epsilon/kb parameter of the interatomic potential, expressed in K;
        determines the amplitude.
    par_s : Float
        sigma parameter of the interatomic potential, expressed in m;
        determines the scalelenght of the interaction.

    Returns
    -------
    Float / Numpy.ndarray of floats
        second virial coefficient, expressed in cm^3/mol.

    """

    def integrand(radius: float, temp: float or np.ndarray,
                  par_m: int, par_e: float, par_s: float):
        amplitude = (par_m*par_e) / (par_m-6) * (par_m/6) ** (6 / (par_m-6))
        attractive = (par_s/radius) ** 6
        repulsive = (par_s/radius) ** par_m
        potential = amplitude * (repulsive - attractive)
        return 1 - np.exp(-potential/temp)
    int_0_inf = (lambda t: integrate.quad(
            (lambda x: integrand(x, t, par_m, par_e, par_s)), 0.1, np.inf)[0])
    if isinstance(temp, float):
        return 0.5 * int_0_inf(temp)
    return np.array([0.5 * int_0_inf(t) for t in temp])

def parabole(temp: float, par_0: float, par_1: float, par_2: float):
    """Simple parabolic function, good enough to fit fast b_aa and b_cc,
    expressed in cm^3/mol."""
    return par_0 + par_1*temp + par_2*temp**2

def exponential(temp: float, par_0: float, par_1: float, par_2: float):
    """exponential function, useful to evaluate b_aa in cm^3/mol,
    taken from Cramer DOI: 10.1121/1.405827.

    Parameters
    ----------
    temp : Float or numpy.ndarray
        Temperature data
    par_0 : Float
        0 Temperature value of b_aa
    par_1 : Float
        Linear b_aa(T) dependance
    par_2 : Float
        Temperature lenght-scale of b_aa

    Returns
    -------
    Float
        exponential description of b_aa(T)

    """
    return par_0 - par_1*np.exp(-par_2/temp)

def hyland(temp: float, par_0: float, par_1: float, par_2: float):
    """Approximated function to describe b_ww in cm^3/mol, taken from hyland
    DOI: 10.6028/jres.079A.017.

    Parameters
    ----------
    temp : Float or numpy.ndarray
        Temperature data
    par_0 : Float
        0K Temperature value of b_ww
    par_1 : Float
        Linear decrease b_aa(T) dependance
    par_2 : Float
        Temperature lenght-scale of b_ww

    Returns
    -------
    Float
        exponential description of b_ww(T)

    """
    return par_0 - par_1/temp * 10**(par_2/temp**2)

def b_aw_fit(temp: float, par_0: float, par_1: float, par_2: float, # pylint: disable=too-many-arguments
             par_3: float, par_4: float, par_5: float):
    """Approximated function to evaluate b_aw(T) in cm^3/mol, taken from
    hyland"""
    return (par_0 + par_1*temp + par_2*temp**2 + par_3*temp**3
            + par_4*temp**4 + par_5*temp**5)

def b_fit(b_data, draw: bool = False):
    """Function which is executed as an instance of the class is created:
    it reads data from databases and fits them with the appropriate
    functions, to store the optimized parameters.

    Parameters
    ----------
    draw : bool, optional
        If True, the function will draw plots of the fits, default = False.

    Returns
    -------
    optimized_parameters: 3x3 np.ndarray of floats
        Optimal values for the parameters so that the sum of the squared
        residuals of f(xdata, *popt) - ydata is minimized; of respectively
        b_aa, b_ww, b_cc, b_aw.

    optimized_covariances: 3x3x3 np.ndarray of floats
        The estimated covariance of optimized_parameters.
        The diagonals provide the variance of the parameter estimate.

    """

    function_list = [exponential, hyland, parabole, b_aw_fit]
    p0_list = [[1, 1, 1], [33.97, 55306, 72000], [21, 147, 100],
               [36.98928, -0.331705, 0.139035E-2,
                -0.574154E-5, 0.326513E-7,-0.142805E-9]]
    title = ['Air', 'Water', 'CO2', 'Moist']
    optimized_parameters = []
    optimized_deviations = []
    for i, data in enumerate(b_data):
        popt, pcov = curve_fit(function_list[i], data['T,K'], # pylint: disable=unbalanced-tuple-unpacking
                               data['B,cm^3/mol'], p0_list[i],
                               data['uB,cm^3/mol'])
        optimized_parameters.append(popt)
        optimized_deviations.append(np.sqrt(np.diagonal(pcov)))
        if draw:
            new_temp = np.arange(273, 313, 0.5)
            new_b = function_list[i](new_temp, *popt)
            fig = plt.figure(figsize=(6,6))
            axis = fig.add_subplot(111)
            axis.set(xlabel = 'Temperature (K)', ylabel = 'B(T) (cm^3/mol)',
                     title = title[i]+str(popt))
            axis.tick_params(direction = 'in')
            plt.errorbar(x = data['T,K'], y = data['B,cm^3/mol'],
                         yerr=data['uB,cm^3/mol'], fmt='o',
                            color = 'blue', label = 'data')
            sns.lineplot(x = new_temp, y = new_b,
                         color = 'orange', label = 'fit')
    return optimized_parameters, optimized_deviations

def giacomo_func(temp, par_0: float, par_1: float, par_2: float, par_3: float):
    """
    Fitting function for saturated vapor pressure taken from P.Giacomo,
    Equation for the Determination of the Density of Moist Air (1981).
    """
    return np.exp(par_0 +par_1*temp + par_2*temp**2 + par_3*temp**(-1))

def read_svp(path: str):
    """
    Reads data of saturated vapor pressure for T in range 0-26°C and RH in
    range 0-90%, taken from P.Giacomo, Equation for the Determination of
    the Density of Moist Air (1981).

    Returns
    -------
    svp_tidy : pd.DataFrame
        Table with columns Temperature (K), Humidity (%), Pressure_SV (Pa)

    """
    svp_data = pd.read_csv(path, '	')
    svp_tidy = pd.DataFrame()
    t_increments = svp_data.columns.unique()[1:]
    t_increments_wide = np.array(np.tile(t_increments,27),dtype=float)
    t_integers = np.repeat(TEMP_0+svp_data['Temperature(C)'],10)
    svp_tidy['Temperature (K)'] = t_integers + t_increments_wide
    svp_data = svp_data.drop(columns='Temperature(C)')
    svp_tidy['Pressure_SV (Pa)'] = np.concatenate([svp_data.T[col]
                                            for col in svp_data.T.columns])
    return svp_tidy

@dataclass
class Environment ():#pylint: disable=R0902
    """Instances of this class represents a thermodynamic system
    with homogeneous temperature, humidity and pressure.
    It will later implement compatibility with uncertainties for
    errors evaluation and with pyroomacoustics for geometrical considerations.
    Initialization data are stored inside folder 'Data/'.

    Parameters
    ----------
    t_input: Float, optional
        Temperature inside the Environment, expressed in K
    h_input: Float, optional
        Humidity inside the Environment, expressed in %
    p_input: Float, optional
        Pressure inside the Environment, expressed in Pa
    """

    t_input: float = field(default = TEMP_0 + 25, metadata={'unit' : 'K'})
    h_input: float = field(default = 0.0, metadata={'unit' : '%'})
    p_input: float = field(default = ATM_P, metadata={'unit' : 'Pa'})
    molar_fraction = pd.DataFrame()
    svp_results = list
    b_results = list
    molar_mass = float
    xww = float

    def __post_init__(self):
        """
        When the first instance of class Environment is created:
            -it will store data from the paths provided into databases
            -it will fit those data to obtain the right parameters of Bi(T)

        Raises
        ------
        FileNotFoundError
            This error raises if it has not been possible to find the specified
            data in the locations provided.

        """
        if Environment.molar_fraction.empty:
            molar_fraction_path = 'Data/Molar_Fraction.txt'
            air_data_path = 'Data/Dry_Air.txt'
            water_data_path = 'Data/H2O.txt'
            co2_data_path = 'Data/CO2.txt'
            moist_data_path = 'Data/Moist_Air.txt'
            saturated_vapor_path = 'Data/SVP.txt'
            try:
                Environment.molar_fraction = pd.read_csv(molar_fraction_path,
                                                  sep=' ', header=0)
                svp_data = read_svp(saturated_vapor_path)
                b_data = [pd.read_csv(air_data_path, sep=' ', header=0),
                          pd.read_csv(water_data_path, sep=' ', header=1),
                          pd.read_csv(co2_data_path, sep=' ', header=0),
                          pd.read_csv(moist_data_path, sep=' ', header=3)]
                #Added https://pubs.acs.org/doi/pdf/10.1021/je60069a019
                Environment.b_results = b_fit(b_data)
                Environment.svp_results = curve_fit(giacomo_func,
                                                    svp_data['Temperature (K)'],
                                                    svp_data["Pressure_SV (Pa)"],
                                                    [34, -2E-2, 1E-5, 6E3])
                self.xww = self.rh_to_xw()
                Environment.molar_mass = self.set_molar_mass()
            except FileNotFoundError as fnne:
                raise FileNotFoundError("Check that input data are in " +
                                        "Data/ folder.") from fnne
            except (ValueError, RuntimeError):
                print("Fitting not succeded, something wrong in the data.")
                raise
        else:
            self.xww = self.rh_to_xw()

    def set_molar_mass(self):
        """Evaluates molar mass of the air components from Molar Fraction
        Database, taking into account the most relevant components.
        Called as an instance of the class is created.

        Returns
        -------
        mass : float
            Molar mass of dry air CO2 free, water vapor and CO2 in kg/mol
        """
        molecules = ['N2', 'O2', 'Ar', 'Ne','CO']
        emf = Environment.molar_fraction
        molar_masses_xi = np.array([float(emf[emf['Constituent'] == m]['xiMi'])
                                    for m in molecules])
        x_ad = [float(emf[emf['Constituent'] == m]['xi']) for m in molecules]
        air_molar_mass = np.sum(molar_masses_xi) / np.sum(x_ad)
        water_molar_mass = float(emf[emf['Constituent'] == 'H2O']['Mi'])
        co2_molar_mass = float(emf[emf['Constituent'] == 'CO2']['Mi'])
        xcc = float(emf[emf['Constituent'] == 'CO2']['xi'])
        xaa = 1-self.xww-xcc
        mass = air_molar_mass*xaa + water_molar_mass*self.xww + co2_molar_mass*xcc
        return 1E-3*mass

    def rh_to_xw(self):
        """Conversion from relative humidity in % to water molar fraction,
        taken from 'A simple expression for the saturation vapour pressure of
        water in the range −50 to 140°C' J M Richards
        """
        svp_pars, svp_cov = self.svp_results
        svp_err = np.sqrt(np.diag(svp_cov))# pylint: disable=unused-variable
        saturated_vapor_p = giacomo_func(self.t_input, *svp_pars)
        enhance_f = (1.00062
                    + 3.14E-8*self.p_input
                    + 5.6E-7*(self.t_input-TEMP_0)**2
                    ) #Enhance factor, error = 1E-4
        return enhance_f * self.h_input/100 * saturated_vapor_p/self.p_input

    def b_mix_function(self, temp: float or np.ndarray): #Decreases too much with RH
        """Evaluates the composed B at the Environment temperature.
        Something about the theory.

        Parameters
        ----------
        temp : Float or numpy.ndarray
            Temperature data in K

        Returns
        -------
        b_mix : float or np.ndarray
            Second virial coefficient of the mixed gas, in m^3/mol
        """
        emf = Environment.molar_fraction
        eb_vals, eb_dev  = Environment.b_results# pylint: disable=unused-variable
        xcc = float(emf[emf['Constituent'] == 'CO2']['xi'])
        xaa = (1-xcc-self.xww)
        b_aa = exponential(temp, *eb_vals[0])
        b_ww = hyland(temp, *eb_vals[1])
        b_cc = parabole(temp, *eb_vals[2])
        b_aw = b_aw_fit(temp, *eb_vals[3])
        b_mix = b_aa*xaa**2 + b_cc*xcc**2 + 2*b_ww*self.xww**2 + 5E2*b_aw*xaa*self.xww
        return 1E-6*b_mix         #conversion from cm^3/mol to m^3/mol

    def gamma_adiabatic(self):
        """This function determines the heat capacity ratio γ for the
        Environment, starting from its temperature, pressure, molar mass,
        B coefficient and its derivatives up to second order.

        Returns
        -------
        gamma Float
            Heat capacity ratio, adimensional.

        """
        b_prime = misc.derivative(self.b_mix_function,
                                  self.t_input, dx=0.1, n=1) #m^3 mol^-1 K^-1
        b_second = misc.derivative(self.b_mix_function,
                                   self.t_input, dx=0.1, n=2) #m^3 mol^-1 K^-2
        b_second_0 = misc.derivative(self.b_mix_function,TEMP_0, dx=0.1, n=2)
        mass = Environment.molar_mass #kg mol^-1
        cp1 = CP1_STP - self.p_input/mass * (self.t_input*b_second
                                             - TEMP_0*b_second_0)
        cv1 = cp1 - (R + 2*self.p_input*b_prime) / mass
        gamma = cp1/cv1
        return gamma

    def sound_speed_0(self):
        """Evaluates 0 frequency sound speed c approximated to the second
        virial term, following Cramer's tractation.
        Returns
        -------
        c_0 : float or np.ndarray
            c at 0 frequency, in m/s
        """
        b_temp = self.b_mix_function(self.t_input)
        gamma = self.gamma_adiabatic()
        mass = Environment.molar_mass
        c_0 = np.sqrt(gamma/mass * (self.t_input*R + 2*self.p_input*b_temp))
        return c_0

    def freq_relax_nitro(self):
        """
        Evaluates the relaxation frequency of nitrogen at the environmental
        conditions.

        Returns
        -------
        f_nitro : float
            Frequency of the wave in Hz

        """
        t_red = self.t_input/(TEMP_0+20)
        p_red = self.p_input/ATM_P
        h_enh = self.xww*1E3
        f_nitro = p_red*t_red**(-0.5)*(9 + 28*h_enh*np.exp(-4.170*(t_red**(-1/3)-1)))
        return f_nitro

    def freq_relax_oxy(self):
        """
        Evaluates the relaxation frequency of nitrogen at the environmental
        conditions.

        Returns
        -------
        f_oxy : float
            Frequency of the wave in Hz

        """
        p_red = self.p_input/ATM_P
        h_enh = self.xww*1E3
        f_oxy = p_red*(24+4.04E3*h_enh*(0.2+h_enh)/(3.91+h_enh))
        return f_oxy

    def attenuation_corrections(self, frequency: float or np.ndarray):
        """
        Evaluates the corrections to speed of sound due to attenuation of
        oxygen and nitrogen, for a given frequency. These corrections are
        expressed as:
            α_i/2πf_i with i=N,O
        where α_i is the attenuation coefficient and f_i the attenuation
        frequency of the chemical species.

        Parameters
        ----------
        frequency : floar or np.ndarray
            Frequency of the wave in Hz

        Returns
        -------
        corrections : list
            List of the corrections for N and O, in s/m

        """
        relax_freqs = np.array([self.freq_relax_nitro(),self.freq_relax_oxy()])
        t_red = np.array([3352, 2239.1]) / self.t_input
        coeffs = np.array([0.781, 0.209]) / 35
        if isinstance(frequency, np.ndarray):
            frequency = frequency.reshape(len(frequency),1)
            relax_freqs = relax_freqs.reshape(1,2)
            t_red = t_red.reshape(1,2)
            coeffs = coeffs.reshape(1,2)
        f_red = frequency/relax_freqs
        corrections = (coeffs*f_red/ self.sound_speed_0() * f_red/(1+f_red**2)
                         * t_red**2 * np.exp(-t_red))
        return corrections
    
    def attenuation_factor(self, frequency:np.ndarray):
        """
        Atmospheric atenuation factor for the environment, evaluated for a
        range of frequencies. 

        Parameters
        ----------
        frequency : np.ndarray
            Frequency of the wave in Hz
            
        Returns
        -------
        attenuation : np.ndarray
            Attenuation of the wave for each frequency in input, in m^-1.

        """
        relax_freqs = np.array([self.freq_relax_nitro(),
                                self.freq_relax_oxy()]).reshape(1,2)
        absorption_coeff = self.attenuation_corrections(frequency)*relax_freqs
        absorption = 2 * np.pi* absorption_coeff.sum(axis=1)
        attenuation = 0.001 * absorption
        return attenuation

    def sound_speed_f(self, frequency: float = 0 or np.ndarray):
        """
        Evaluates frequency dependant speed of sound

        Parameters
        ----------
        frequency : float or np.ndarray
            Frequency of the wave in Hz, if 0 (default) it calls directly
            function Environment.sound_speed_0()

        Returns
        -------
        sound_speed : float or np.ndarray
            c at the input frequency, in m/s

        """
        if isinstance(frequency,float):
            if frequency==0:
                return self.sound_speed_0()
            c_nitro, c_oxy = self.attenuation_corrections(frequency)
        else:
            c_nitro, c_oxy = self.attenuation_corrections(frequency).T
        sound_speed = 1 / (1/self.sound_speed_0() - c_nitro - c_oxy)
        return sound_speed
    