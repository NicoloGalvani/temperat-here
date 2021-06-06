# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 10:05:06 2021

@author: Nicolò Galvani
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


"""Environment class simulates a volume of air in which temperature, pressure
and humidity are constant and homogeneous. Its elements take these parameters
as inputs, to determine a set of properties for the gas, in particular the
speed of sound.
"""

TEMP_0 = 273.15 #K
ATM_P = 101325 #Pa
CP0 = 1.006E3 #j Kg^-1 K^-1
R = 8.31446261815324 #m^3 Pa K^-1 mol^-1
KL_URL = ('https://web.archive.org/web/20190508003406/http://www.kayelaby.npl'
          +'.co.uk/general_physics/2_4/2_4_1.html#speed_of_sound_in_air')

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
            (lambda x: integrand(x, t, par_m, par_e, par_s)), 0.01, np.inf)[0])
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

def b_aw_fit(temp: float or np.ndarray):
    """Approximated function to evaluate b_aw(T) in cm^3/mol, taken from
    hyland"""
    return (
            36.98928
            - 0.331705*temp
            + 0.139035E-2*temp**2
            - 0.574154E-5*temp**3
            + 0.326513E-7*temp**4
            - 0.142805E-9*temp**5
            )

def b_fit(draw: bool = True):
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
        b_aa, b_ww, b_cc.

    optimized_covariances: 3x3x3 np.ndarray of floats
        The estimated covariance of optimized_parameters.
        The diagonals provide the variance of the parameter estimate.

    """

    data_list = Environment.b_data
    function_list = [exponential, hyland, parabole]
    p0_list = [[1, 1, 1], [33.97, 55306, 72000], [21, 147, 100]]
    title = ['Air', 'Water', 'CO2']
    optimized_parameters = []
    optimized_covariances = []
    for i, data in enumerate(data_list):
        popt, pcov = curve_fit(function_list[i], data['T,K'], data['B,cm^3/mol'],
                               p0_list[i])
        optimized_parameters.append(popt)
        optimized_covariances.append(pcov)
        if draw:
            new_temp = np.arange(273, 313, 0.5)
            new_b = function_list[i](new_temp, *popt)
            plt.figure(figsize=(6,6))
            axis = plt.axis()
            axis.set(xlabel = 'Temperature (K)', ylabel = 'B(T) (cm^3/mol)',
                     title = title[i]+str(popt))
            axis.tick_params(direction = 'in')
            sns.scatterplot(x = data['T,K'], y = data['B,cm^3/mol'],
                            color = 'blue', label = 'data')
            sns.lineplot(x = new_temp, y = new_b,
                         color = 'orange', label = 'fit')
            # ax.text(1, 1, str(popt), fontsize=10)
    return optimized_parameters, optimized_covariances

@dataclass
class Environment ():
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
    b_data = pd.DataFrame()
    b_values = pd.DataFrame()
    b_covariances = pd.DataFrame()
    molar_mass = float

    ##Possibility to insert also the link to pyroomsound?

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
            molar_fraction_path = 'Data/molar_fraction.txt'
            air_data_path = 'Data/Dry_Air.txt'
            water_data_path = 'Data/H2O.txt'
            co2_data_path = 'Data/CO2.txt'
            try:
                Environment.molar_fraction = pd.read_csv(
                                            molar_fraction_path, sep=' ')
                Environment.b_data = [pd.read_csv(air_data_path, sep=' '),
                                      pd.read_csv(water_data_path, sep=' '),
                                      pd.read_csv(co2_data_path, sep=' ')]
                #Added https://pubs.acs.org/doi/pdf/10.1021/je60069a019
                Environment.b_values, Environment.b_covariances = b_fit()
                Environment.molar_mass = self.set_molar_mass()
            except FileNotFoundError as fnne:
                raise FileNotFoundError("Unable to find the reference data"
                                        ) from fnne

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
        air_molar_mass = np.sum(molar_masses_xi)/np.sum(x_ad)
        water_molar_mass = float(emf[emf['Constituent'] == 'H2O']['Mi'])
        co2_molar_mass = float(emf[emf['Constituent'] == 'CO2']['Mi'])
        xcc = float(emf[emf['Constituent'] == 'CO2']['xi'])
        xww = self.rh_to_xw()
        xaa = 1-xww-xcc
        mass = air_molar_mass*xaa + water_molar_mass*xww + co2_molar_mass*xcc
        return 1E-3*mass

    def rh_to_xw(self):
        """Conversion from relative humidity in % to water molar fraction,
        taken from 'A simple expression for the saturation vapour pressure of
        water in the range −50 to 140°C' J M Richards
        """
        temp_reduced = 1 - (100 + TEMP_0) / self.t_input
        saturated_vapor_pressure = ATM_P * np.exp(
                                                  13.3185*temp_reduced
                                                  - 1.9760*temp_reduced**2
                                                  - 0.6445*temp_reduced**3
                                                  - 0.1299*temp_reduced**4
                                                  )
        en_f = 1.00062 + 3.14E-8*self.p_input + 5.6E-7*self.t_input**2 #Enhance factor
        return self.h_input/100 * en_f * saturated_vapor_pressure/self.p_input

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
        ebv = Environment.b_values
        xcc = float(emf[emf['Constituent'] == 'CO2']['xi'])
        xww = self.rh_to_xw()
        xaa = (1-xcc-xww)
        b_aa = exponential(temp, *ebv[0])
        b_ww = hyland(temp, 33.97, 55306, 72000)#ebv[1][0],ebv[1][2],ebv[1][2])
        b_cc = parabole(temp, *ebv[2])
        b_aw = b_aw_fit(temp)
        b_mix = b_aa*xaa**2 + b_cc*xcc**2 + b_ww*xww**2 + 2*b_aw*xaa*xww
        return 1E-6*b_mix         #conversion from cm^3/mol to m^3/mol

    def gamma(self):
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
        mass = Environment.molar_mass #kg mol^-1
        cp1 = CP0 - self.p_input * self.t_input * b_second / mass
        cv1 = cp1 - (R + 2*self.p_input*b_prime) / mass
        gamma = cp1/cv1
        return gamma #Decreases too much with RH

    def sound_speed(self):
        """Mockup function, used now for testing, later implemented.
        """
        b_temp = self.b_mix_function(self.t_input)
        gamma = self.gamma()
        mass = Environment.molar_mass
        temp = self.t_input
        c_0 = np.sqrt(gamma/mass * (temp*R + 2*self.p_input*b_temp))
        return c_0
    