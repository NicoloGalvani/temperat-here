"""
Module with the definition of Environment class and the functions needed to
fit experimental data.
Notable dependencies:
    module uncertainties for handling of error propagation

"""


from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.misc as misc
from matplotlib import pyplot as plt
import seaborn as sns
from uncertainties import ufloat
from uncertainties import umath
import uncertainties




TEMP_0 = 273.15 #K
ATM_P = 101325 #Pa
CP0 = 1.006E3 #j Kg^-1 K^-1
R_GAS = 8.31446261815324 #m^3 Pa K^-1 mol^-1



#################Functions needed for fitting##############
#Does not converge for Dry Air, to fix
def lennard (temperature: float or np.ndarray, m_par: int,
             e_par: float, s_par: float):
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
    temperature : Float / Numpy.ndarray of floats
        Temperature of the gas, expressed in K.
    m_par : Int > 0
        m parameter of the interatomic potential, dimensionless quantity;
        determines the entity of repulsive potential and the amplitude.
    e_par : Float
        epsilon/kb parameter of the interatomic potential, expressed in K;
        determines the amplitude.
    s_par : Float
        sigma parameter of the interatomic potential, expressed in m;
        determines the scalelenght of the interaction.

    Returns
    -------
    Float / Numpy.ndarray of floats
        second virial coefficient, expressed in cm^3/mol.

    """

    def integrand(distance: float, temperature: float,# or np.ndarray,
                  m_par: int, e_par: float, s_par: float):
        amplitude = (m_par*e_par) / (m_par-6) * (m_par/6) ** (6 / (m_par-6))
        attractive = (s_par/distance) ** 6
        repulsive = (s_par/distance) ** m_par
        potential = amplitude * (repulsive - attractive)
        return 1 - np.exp(-potential/temperature)

    def integration(temperature: float, m_par: float, e_par:float,s_par:float):
        return 0.5 * integrate.quad(lambda x:
                                    integrand(x, temperature, m_par,
                                              e_par, s_par),
                                    0.01, np.inf)[0]

    if isinstance(temperature, float):
        return integration(temperature, m_par, e_par, s_par)
    return np.array([integration(t_0, m_par, e_par, s_par)
                       for t_0 in temperature])

def parabole(temperature: float, c_0: float, c_1: float, c_2: float):
    """Simple parabolic function, good enough to fit fast Baa and Bcc,
    expressed in cm^3/mol."""
    return c_0 + c_1*temperature + c_2*temperature**2

def exponential(temperature: float, c_0: float, c_1: float, c_2: float):
    """Exponential function, useful to evaluate Baa in cm^3/mol,
    taken from Cramer DOI: 10.1121/1.405827.

    Parameters
    ----------
    temperature : Float or numpy.ndarray
        Temperature data
    c_0: Float
        0 Temperature value of Baa
    c_1 : Float
        Linear Baa(T) dependance
    c_2 : Float
        Temperature lenght-scale of Baa

    Returns
    -------
    Float
        Exponential description of Baa(T)

    """
    return c_0 - c_1*np.exp(-c_2/temperature)

def hyland_b_ww(temperature: float, c_0: float, c_1: float, c_2: float):
    """Approximated function to describe Bww in cm^3/mol, taken from Hyland
    DOI: 10.6028/jres.079A.017.

    Parameters
    ----------
    temperature : Float or numpy.ndarray
        Temperature data
    c_0: Float
        0K Temperature value of Bww
    c_1 : Float
        Linear decrease Baa(T) dependance
    c_2 : Float
        Temperature lenght-scale of Bww

    Returns
    -------
    Float
        Exponential description of Bww(T)

    """
    return c_0 - c_1/temperature * 10**(c_2/temperature**2)

def hyland_b_aw(temperature: float or np.ndarray):
    """Approximated function to evaluate Baw(T) in cm^3/mol, taken from
    Hyland DOI: 10.6028/jres.079A.017.
    Deprecated, substituted by data from Hellmann
    DOI: 10.1021/acs.jced.0c00465"""
    return (
            36.98928
            - 0.331705*temperature
            + 0.139035E-2*temperature**2
            - 0.574154E-5*temperature**3
            + 0.326513E-7*temperature**4
            - 0.142805E-9*temperature**5
            )

def b_graph(x_data, y_data, popt, function, title):
    """Functions which plots virial coefficients fit for a determined dataset


    Parameters
    ----------
    x_data :
        Data for the temperatures.
    y_data : TYPE
        Data for the virial coefficients.
    popt : list of float
        Optimized parameters obtained fitting the function on the dataset.
    function : function
        DESCRIPTION.
    title : string
        Title to apply to the graph.

    Returns
    -------
    None.

    """
    new_x = np.arange(x_data[0], x_data[len(x_data)-1], 0.5)
    new_y = function(new_x, popt[0], popt[1], popt[2])
    fig, axis = plt.subplots()
    axis.tick_params(direction = 'in')
    fig.set_figheight(6)
    fig.set_figwidth(6)
    sns.scatterplot(x = x_data, y = y_data,
                    color = 'blue', label = 'data')
    sns.lineplot(x = new_x, y = new_y,
                 color = 'orange', label = 'fit')
    plt.title(title+str(popt))
    axis.set_ylabel('B(T) (cm^3/mol)')
    axis.set_xlabel('Temperature (K)')

def b_fit():
    """Function which is executed as an instance of the class is created:
    it reads data from databases and fits them with the appropriate
    functions, to store the optimized parameters.

    Returns
    -------
    optimized_parameters: 3D np.ndarray of floats
        Optimal values for the parameters so that the sum of the squared
        residuals of f(xdata, *popt) - ydata is minimized; of respectively
        Baa, Bww, Bcc.

    optimized_deviations: 3D np.ndarray of floats
        The estimated standard deviations of optimized_parameters.

    """

    functions = [exponential, hyland_b_ww, parabole, hyland_b_ww]
    p_0 = [[1, 1, 1], [33.97, 55306, 72000], [21, 147, 100], [21, 147, 100]]
    title = ['Air', 'Water', 'CO2', 'Moist']
    optimized_parameters = []
    optimized_deviations = []
    for i, data in enumerate(Environment.virial_data):
        function = functions[i]
        x_vector = data['T,K']
        y_vector = data['B,cm^3/mol']
        err = data['uB,cm^3/mol']
        popt, pcov = curve_fit(function, x_vector, y_vector, p_0[i], sigma=err) # pylint: disable=unbalanced-tuple-unpacking
        optimized_parameters.append(popt)
        optimized_deviations.append(np.sqrt(np.diag(pcov))/popt)
        if Environment.draw_b_plots_in_b_fit:
            b_graph(x_vector, y_vector, popt, function, title[i])
    return optimized_parameters, optimized_deviations

###############Class Definition###############
@dataclass
class Environment ():
    # pylint: disable=too-many-instance-attributes
    # pylint interprets as instance attributes the parameters handled by
    # dataclass, thus it raises the warning 'too-many-instance-attributes'
    """Instances of this class represents a thermodynamic system
    with homogeneous temperature, humidity and pressure.
    It will later implement compatibility with uncertainties for
    errors evaluation and with pyroomacoustics for geometrical considerations.

    Parameters
    ----------
    t_input: Float, optional
        Temperature inside the environment, expressed in K
    h_input: Float, optional
        Humidity inside the environment, expressed in %
    p_input: Float, optional
        Pressure inside the environment, expressed in Pa
    input_path: Path, optional
        Folder containing experimental data, written as csv-space-separated

    Attributes
    ----------
    molar_fraction : pd.DataFrame()
         Dataframe storing data from molar_fraction_path
    virial_data: list
        List of the four Dataframe containing data of virial coefficients
    draw_b_plots_in_b_fit : bool
        If True, b_fit method will plot graphs comparing experimental data and
        fitting curve for Baa, Bww, Bcc
    b_values : pd.DataFrame()
        Dataframe storing b_fit optimal parameters
    b_deviations : pd.DataFrame()
        Dataframe storing b_fit optimal parameters standard deviations
    molar_mass : float
        Total molar mass of the gas, evaluated through set_molar_mass

    """



    t_input: ufloat = field(default = ufloat(TEMP_0 + 25, .5, 'Temperature'),
                            metadata={'unit' : 'K'})
    h_input: ufloat = field(default = ufloat(0.0, 5, 'Humidity'),
                            metadata={'unit' : '%'})
    p_input: ufloat = field(default = ufloat(ATM_P, 100, 'Pressure'),
                            metadata={'unit' : 'Pa'})
    input_path: Path = field(default = 'Data')
    molar_fraction = pd.DataFrame()
    virial_data = list
    draw_b_plots_in_b_fit = False
    b_values = pd.DataFrame()
    b_deviations = pd.DataFrame()
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
            try:
                Environment.molar_fraction = pd.read_csv(self.input_path +
                                                         '/molar_fraction.txt',
                                                         sep=' ')
                air_data = pd.read_csv(self.input_path + '/Dry_Air.txt',
                                       sep=' ', header=2)
                water_data = pd.read_csv(self.input_path + '/H2O.txt',
                                         sep=' ', header=1)
                co2_data = pd.read_csv(self.input_path + '/CO2.txt',
                                       sep=' ', header=3)
                moist_air_data = pd.read_csv(self.input_path +'/Moist_Air.txt',
                                             sep=' ', header=3)
                Environment.virial_data = (air_data, water_data,
                                           co2_data, moist_air_data)
            except FileNotFoundError as error:
                raise FileNotFoundError("Unable to find the input data"
                                        ) from error
            Environment.b_values, Environment.b_deviations = b_fit()
            Environment.molar_mass = self.set_molar_mass()

    def set_molar_mass(self):
        """Evaluates molar mass of the air components from Molar Fraction
        Database, taking into account the most relevant components.
        Called as an instance of the class is created.

        Returns
        -------
        Mass : float
            Molar mass of dry air CO2 free, water vapor and CO2 in kg/mol
        """
        molecules = ['N2', 'O2', 'Ar', 'Ne','CO']
        emf = Environment.molar_fraction
        molar_masses_xi = np.array([float(emf[emf['Constituent'] == m]['xiMi'])
                                    for m in molecules])
        x_data = [float(emf[emf['Constituent'] == m]['xi']) for m in molecules]
        air_molar_mass = np.sum(molar_masses_xi)/np.sum(x_data)
        water_molar_mass = float(emf[emf['Constituent'] == 'H2O']['Mi'])
        co2_molar_mass = float(emf[emf['Constituent'] == 'CO2']['Mi'])
        x_carbon = float(emf[emf['Constituent'] == 'CO2']['xi'])
        x_water = self.rh_to_xw()
        x_air = 1-x_water-x_carbon
        mass = (air_molar_mass*x_air
                + water_molar_mass*x_water
                + co2_molar_mass*x_carbon)
        return 1E-3*mass

    def rh_to_xw(self, temperature : float = 0.):
        """Conversion from relative humidity in % to water molar fraction,
        taken from 'c_0simple expression for the saturation vapour pressure of
        water in the range −50 to 140°C' J M Richards

        Parameters
        ----------
        temp : Float
            Temperature data in K, to be substituted to self.t_input if
            rh_to_xw is called by a function which want to study different
            temperatures (b_mix_function)

        Returns
        -------
        xw : uFloat
            Water molar fraction in the Environment


        """
        if temperature == 0.:
            temperature = self.t_input
        t_reduced = 1 - (100 + TEMP_0) / temperature
        u_exp = uncertainties.wrap(np.exp)
        saturated_vapor_pressure = ufloat(ATM_P,0) * u_exp(
                                                        13.3815*t_reduced
                                                        - 1.9760*t_reduced**2
                                                        - 0.6445*t_reduced**3
                                                        - 0.1299*t_reduced**4
                                                        )  #Approximation good up to 0.1%
        intrinsic_error = saturated_vapor_pressure.n*0.001
        svp = ufloat(saturated_vapor_pressure.n,
                     intrinsic_error + saturated_vapor_pressure.s)
        enhance_factor = 1.00062 + 3.14E-8*self.p_input + 5.6E-7*(temperature
                                                                  -TEMP_0)**2
        x_water = self.h_input/100 * enhance_factor * svp/self.p_input
        x_water = ufloat(x_water.n, x_water.s, 'Water molar mass')
        return x_water

    def b_mix_function(self, temperature: float = 0.): #Decreases too much with RH
        """Evaluates the composed B at the environment temperature.
        Something about the theory.

        Parameters
        ----------
        temperature : Float
            Temperature data in K

        Returns
        -------
        b_mix : ufloa
            Second virial coefficient of the mixed gas, in m^3/mol
        """
        if temperature == 0.:
            temperature = self.t_input
        ebv = Environment.b_values
        x_carbon = float(Environment.molar_fraction
                         [Environment.molar_fraction[
                             'Constituent'] == 'CO2']['xi'])
        x_water = self.rh_to_xw(temperature)
        x_air = (1-x_carbon-x_water)
        fit_b_aa = uncertainties.wrap(exponential)
        fit_b_ww = uncertainties.wrap(hyland_b_ww)
        fit_b_cc = uncertainties.wrap(parabole)
        fit_b_aw = uncertainties.wrap(hyland_b_ww)
        b_aa = fit_b_aa(temperature,ebv[0][0],ebv[0][1],ebv[0][2])
        b_ww = fit_b_ww(temperature,ebv[1][0],ebv[1][2],ebv[1][2])
        b_cc = fit_b_cc(temperature,ebv[2][0],ebv[2][1],ebv[2][2])
        b_aw = fit_b_aw(temperature,ebv[3][0],ebv[3][1],ebv[3][2])
        b_mix = (b_aa*x_air**2
                + b_cc*x_carbon**2
                + b_ww*x_water**2
                + 2*b_aw*x_air*x_water)
        #conversion from cm^3/mol to m^3/mol
        b_mix = ufloat(1E-6*b_mix.n, 1E-6*b_mix.s, 'Virial coefficient')
        return b_mix

    def gamma_adiabatic(self):
        """This function determines the heat capacity ratio γ for the
        Environment, starting from its temperature, pressure, molar mass,
        B coefficient and its derivatives up to second order.

        Returns
        -------
        gamma : Float or np.ndarray of Float
            Heat capacity ratio, adimensional.

        """
        #gamma_adiabatic too high
        b_prime = misc.derivative(self.b_mix_function,
                                  self.t_input, dx=0.1, n=1) #m^3 mol^-1 K^-1
        b_second = misc.derivative(self.b_mix_function,
                                   self.t_input, dx=0.1, n=2) #m^3 mol^-1 K^-2
        mass = Environment.molar_mass #kg mol^-1
        cp1 = CP0 - self.p_input * self.t_input * b_second / mass
        cv1 = cp1 - (R_GAS + 2*self.p_input * b_prime) / mass
        gamma = cp1/cv1
        gamma = ufloat(gamma.n, gamma.s, 'gamma_adiabatic')
        return gamma

    def sound_speed(self):
        """Evaluates 0 frequency sound speed c approximated to the second
        virial term, following Cramer's tractation.

        Returns
        -------
        c0 : float or np.ndarray
            c at 0 frequency, in m/s
        """
        temperature = self.t_input
        b_tot = self.b_mix_function(temperature)
        gamma = self.gamma_adiabatic()
        mass = Environment.molar_mass
        c_0 = umath.sqrt(gamma/mass * # pylint: disable=no-member
                          (temperature*R_GAS + 2*self.p_input*b_tot))
        return c_0
