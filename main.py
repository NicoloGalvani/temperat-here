import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.misc as misc
from matplotlib import pyplot as plt
import seaborn as sns



"""This is the module docstring
"""

TEMP_0 = 273.15 #K
ATM_P = 101325 #Pa
CP0 = 1.006E3 #j Kg^-1 K^-1
R = 8.31446261815324 #m^3 Pa K^-1 mol^-1
kayelaby = 'https://web.archive.org/web/20190508003406/http://www.kayelaby.npl.co.uk/general_physics/2_4/2_4_1.html#speed_of_sound_in_air'


@dataclass
class Environment ():
    """Instances of this class represents a thermodynamic system
    with homogeneous temperature, humidity and pressure.
    It will later implement compatibility with uncertainties for
    errors evaluation and with pyroomacoustics for geometrical considerations.

    Parameters
    ----------
    Molar_Fraction_path: Path, optional
        Path containing air composition data, written as a csv-space-separated
    Air_Data_path: Path, optional
        Path containing CO2-free-dry-air B parameter at varying T, written as
        a csv-space-separated.
    Water_Data_path: Path, optional
        Path containing water vapor B parameter at varying T, written as
        a csv-space-separated.
    CO2_Data_path: Path, optional
        Path containing CO2 B parameter at varying T, written as a
        csv-space-separated.
    Moist_Air_Data_path: Path, optional
        Path containing moist-air Baw parameter at varying T, written as a
        csv-space-separated.
    T_input: Float, optional
        Temperature inside the environment, expressed in K
    H_input: Float, optional
        Humidity inside the environment, expressed in %
    P_input: Float, optional
        Pressure inside the environment, expressed in Pa

    Attributes
    ----------
    Molar_Fraction : pd.DataFrame()
        Dataframe storing data from Molar_Fraction_path
    Air_Data : pd.DataFrame()
        Dataframe storing data from Air_Data_path
    Water_Data : pd.DataFrame()
        Dataframe storing data from Water_Data_path
    CO2_Data : pd.DataFrame()
        Dataframe storing data from CO2_Data_path
    Moist_Air_Data: pd.DataFrame()
        Dataframe storing data from Moist_Air_Data_path
    draw_B_plots_in_B_fit : bool
        If True, B_fit method will plot graphs comparing experimental data and
        fitting curve for Baa, Bww, Bcc
    B_values : pd.DataFrame()
        Dataframe storing B_fit optimal parameters
    B_covariances : pd.DataFrame()
        Dataframe storing B_fit optimal parameters covariances
    Molar_Mass : float
        Total molar mass of the gas, evaluated through Set_Molar_Mass

    """



    T_input: float = field(default = TEMP_0 + 25, metadata={'unit' : 'K'})
    H_input: float = field(default = 0.0, metadata={'unit' : '%'})
    P_input: float = field(default = ATM_P, metadata={'unit' : 'Pa'})
    Molar_Fraction_path: Path = field(default = 'Data/Molar_Fraction.txt')
    Air_Data_path: Path = field(default = 'Data/Dry_Air.txt')
    Water_Data_path: Path = field(default = 'Data/H2O.txt')
    CO2_Data_path: Path = field(default = 'Data/CO2.txt')
    Moist_Air_Data_path: Path = field(default = 'Data/Moist_Air.txt')
    Molar_Fraction = pd.DataFrame()
    Air_Data = pd.DataFrame()
    Water_Data = pd.DataFrame()
    CO2_Data = pd.DataFrame()
    Moist_Air_Data = pd.DataFrame()
    draw_B_plots_in_B_fit = False
    B_values = pd.DataFrame()
    B_covariances = pd.DataFrame()
    Molar_Mass = float

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
        if Environment.Molar_Fraction.empty:
            try:
                Environment.Molar_Fraction = pd.read_csv(self.Molar_Fraction_path, sep=' ')
                Environment.Air_Data = pd.read_csv(self.Air_Data_path, sep=' ')
                Environment.Water_Data = pd.read_csv(self.Water_Data_path, sep=' ',header=1)
                Environment.CO2_Data = pd.read_csv(self.CO2_Data_path, sep=' ')
                Environment.Moist_Air_Data = pd.read_csv(self.Moist_Air_Data_path, sep=' ',header=3)
                (Environment.B_values,
                 Environment.B_covariances) = Environment.B_fit()
                Environment.Molar_Mass = self.Set_Molar_Mass()
            except:
                raise FileNotFoundError("Unable to find the reference data")



#Does not converge for Dry Air, to fix
    def Lennard (T: float or np.ndarray, m: int, e: float, s: float):
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
        T : Float / Numpy.ndarray of floats
            Temperature of the gas, expressed in K.
        m : Int > 0
            m parameter of the interatomic potential, dimensionless quantity;
            determines the entity of repulsive potential and the amplitude.
        e : Float
            epsilon/kb parameter of the interatomic potential, expressed in K;
            determines the amplitude.
        s : Float
            sigma parameter of the interatomic potential, expressed in m;
            determines the scalelenght of the interaction.

        Returns
        -------
        Float / Numpy.ndarray of floats
            second virial coefficient, expressed in cm^3/mol.

        """

        def integrand(r: float, T: float or np.ndarray,
                      m: int, e: float, s: float):
                amplitude = (m*e) / (m-6) * (m/6) ** (6 / (m-6))
                attractive = (s/r) ** 6
                repulsive = (s/r) ** m
                potential = amplitude * (repulsive - attractive)
                return 1 - np.exp(-potential/T)

        if isinstance(T, float):
            return 0.5 * integrate.quad(lambda x: integrand(x, T, m, e, s),
                                        0.01,
                                        np.inf
                                        )[0]
        else:
            return np.array([0.5 * integrate.quad(
                                        lambda x: integrand(x, t, m, e, s),
                                        0.01,
                                        np.inf
                                        )[0] for t in T])

    def Parabole(T: float, A: float, B: float, C: float):
        """Simple parabolic function, good enough to fit fast Baa and Bcc,
        expressed in cm^3/mol."""
        return A + B*T + C*T**2

    def Exponential(T: float, A: float, B: float, C: float):
        """Exponential function, useful to evaluate Baa in cm^3/mol,
        taken from Cramer DOI: 10.1121/1.405827.

        Parameters
        ----------
        T : Float or numpy.ndarray
            Temperature data
        A : Float
            0 Temperature value of Baa
        B : Float
            Linear Baa(T) dependance
        C : Float
            Temperature lenght-scale of Baa

        Returns
        -------
        Float
            Exponential description of Baa(T)

        """
        return A - B*np.exp(-C/T)

    def Hyland(T: float, A: float, B: float, C: float):
        """Approximated function to describe Bww in cm^3/mol, taken from Hyland
        DOI: 10.6028/jres.079A.017.

        Parameters
        ----------
        T : Float or numpy.ndarray
            Temperature data
        A : Float
            0K Temperature value of Bww
        B : Float
            Linear decrease Baa(T) dependance
        C : Float
            Temperature lenght-scale of Bww

        Returns
        -------
        Float
            Exponential description of Bww(T)

        """
        return A - B/T * 10**(C/T**2)

    def B_aw(T: float or np.ndarray):
        """Approximated function to evaluate Baw(T) in cm^3/mol, taken from
        Hyland DOI: 10.6028/jres.079A.017.
        Deprecated, substituted by data from Hellmann
        DOI: 10.1021/acs.jced.0c00465"""
        return (
                36.98928
                - 0.331705*T
                + 0.139035E-2*T**2
                - 0.574154E-5*T**3
                + 0.326513E-7*T**4
                - 0.142805E-9*T**5
                )

    def B_fit():
        """Function which is executed as an instance of the class is created:
        it reads data from databases and fits them with the appropriate
        functions, to store the optimized parameters.

        Returns
        -------
        Optimized_parameters: 3x3 np.ndarray of floats
            Optimal values for the parameters so that the sum of the squared
            residuals of f(xdata, *popt) - ydata is minimized; of respectively
            Baa, Bww, Bcc.

        Optimized_covariances: 3x3x3 np.ndarray of floats
            The estimated covariance of Optimized_parameters.
            The diagonals provide the variance of the parameter estimate.

        """

        Data = [Environment.Air_Data,
                Environment.Water_Data,
                Environment.CO2_Data,
                Environment.Moist_Air_Data,]
        functions = [Environment.Exponential,
                     Environment.Hyland,
                     Environment.Parabole,
                     Environment.Hyland]
        p0 = [[1, 1, 1], [33.97, 55306, 72000], [21, 147, 100], [21, 147, 100]]
        title = ['Air', 'Water', 'CO2', 'Moist']
        Optimized_parameters = []
        Optimized_covariances = []
        for i,D in enumerate(Data):
            function = functions[i]
            x = D['T,K']
            y = D['B,cm^3/mol']
            popt, pcov = curve_fit(function, x, y, p0[i])
            Optimized_parameters.append(popt)
            Optimized_covariances.append(pcov)
            if Environment.draw_B_plots_in_B_fit:
                new_x = np.arange(x[0], x[len(x)-1], 0.5)
                new_y = function(new_x, popt[0], popt[1], popt[2])
                fig,ax = plt.subplots()
                ax.tick_params(direction = 'in')
                fig.set_figheight(6)
                fig.set_figwidth(6)
                sns.scatterplot(x = x, y = y,
                                color = 'blue', label = 'data')
                sns.lineplot(x = new_x, y = new_y,
                             color = 'orange', label = 'fit')
                plt.title(title[i]+str(popt))
                ax.set_ylabel('B(T) (cm^3/mol)')
                ax.set_xlabel('Temperature (K)')
        return Optimized_parameters, Optimized_covariances



    def Set_Molar_Mass(self):
        """Evaluates molar mass of the air components from Molar Fraction
        Database, taking into account the most relevant components.
        Called as an instance of the class is created.

        Returns
        -------
        Mass : float
            Molar mass of dry air CO2 free, water vapor and CO2 in kg/mol
        """
        Molecules = ['N2', 'O2', 'Ar', 'Ne','CO']
        EMF = Environment.Molar_Fraction
        Molar_Masses_xi = np.array([float(EMF[EMF['Constituent'] == m]['xiMi'])
                                    for m in Molecules])
        xad = [float(EMF[EMF['Constituent'] == m]['xi']) for m in Molecules]
        Air_Molar_Mass = np.sum(Molar_Masses_xi)/np.sum(xad)
        Water_Molar_Mass = float(EMF[EMF['Constituent'] == 'H2O']['Mi'])
        CO2_Molar_Mass = float(EMF[EMF['Constituent'] == 'CO2']['Mi'])
        xcc = float(EMF[EMF['Constituent'] == 'CO2']['xi'])
        xww = self.RH_to_Xw()
        xaa = 1-xww-xcc
        Mass = Air_Molar_Mass*xaa + Water_Molar_Mass*xww + CO2_Molar_Mass*xcc
        return 1E-3*Mass

    def RH_to_Xw(self):
        """Conversion from relative humidity in % to water molar fraction,
        taken from 'A simple expression for the saturation vapour pressure of
        water in the range −50 to 140°C' J M Richards
        """
        t = 1 - (100 + TEMP_0) / self.T_input
        Saturated_Vapor_Pressure = ATM_P * np.exp(
                                                  13.3185*t
                                                  - 1.9760*t**2
                                                  - 0.6445*t**3
                                                  - 0.1299*t**4
                                                  )
        f = 1.00062 + 3.14E-8*self.P_input + 5.6E-7*self.T_input**2 #Enhance factor
        return self.H_input/100 * f * Saturated_Vapor_Pressure/self.P_input

    def B_mix_function(self, T: float or np.ndarray): #Decreases too much with RH
        """Evaluates the composed B at the environment temperature.
        Something about the theory.

        Parameters
        ----------
        T : Float or numpy.ndarray
            Temperature data in K

        Returns
        -------
        B_mix : float or np.ndarray
            Second virial coefficient of the mixed gas, in m^3/mol
        """
        EMF = Environment.Molar_Fraction
        EBV = Environment.B_values
        xcc = float(EMF[EMF['Constituent'] == 'CO2']['xi'])
        xww = self.RH_to_Xw()
        xaa = (1-xcc-xww)
        Baa = Environment.Exponential(T,EBV[0][0],EBV[0][1],EBV[0][2])
        Bww = Environment.Hyland(T,EBV[1][0],EBV[1][2],EBV[1][2])
        Bcc = Environment.Parabole(T,EBV[2][0],EBV[2][1],EBV[2][2])
        Baw = Environment.Hyland(T,EBV[3][0],EBV[3][1],EBV[3][2])#Environment.B_aw(T)#
        B_mix = Baa*xaa**2 + Bcc*xcc**2 + Bww*xww**2 + 2*Baw*xaa*xww
        return 1E-6*B_mix         #conversion from cm^3/mol to m^3/mol

    def Gamma(self):
        """This function determines the heat capacity ratio γ for the
        Environment, starting from its temperature, pressure, molar mass,
        B coefficient and its derivatives up to second order.

        Returns
        -------
        gamma : Float or np.ndarray of Float
            Heat capacity ratio, adimensional.

        """
        #Gamma too high
        B_prime = misc.derivative(self.B_mix_function,
                                  self.T_input, dx=0.1, n=1) #m^3 mol^-1 K^-1
        B_second = misc.derivative(self.B_mix_function,
                                   self.T_input, dx=0.1, n=2) #m^3 mol^-1 K^-2
        M = Environment.Molar_Mass #kg mol^-1
        cp1 = CP0 - self.P_input * self.T_input * B_second / M
        cv1 = cp1 - (R + 2*self.P_input*B_prime) / M
        gamma = cp1/cv1
        return gamma

    def Sound_Speed(self):
        """Evaluates 0 frequency sound speed c approximated to the second virial
        term, following Cramer's tractation.

        Returns
        -------
        c0 : float or np.ndarray
            c at 0 frequency, in m/s
        """
        B = self.B_mix_function(self.T_input)
        gamma = self.Gamma()
        M = Environment.Molar_Mass
        T = self.T_input
        c0 = np.sqrt(gamma/M * (T*R + 2*self.P_input*B))
        return c0



def Read_Kayelaby_Speed():
    """Function which reads speed of sound data from Kaye and Laby website and
    store them into a dataframe for comparisons.

    Returns
    -------
    Speed_DF : pandas.Dataframe
        Table of values of speed of sound at varying temperature (from
        0 °C to 30 °C) and RH (from 10% to 90%), expressed in m/s

    """
    Speed_Table = pd.read_html(kayelaby)[10].dropna()
    Speed_DF = pd.DataFrame({
                             'Speed (m/s)' : [],
                             'Temperature (°C)' : [],
                             'Relative Humidity (%)' : []
                             })
    for x in range(1,10):
        Fixed_RH = pd.DataFrame({
                                'Speed (m/s)' : Speed_Table[x][2:],
                                'Temperature (°C)' : Speed_Table[0][2:],
                                'Relative Humidity (%)' :
                                    [Speed_Table[x][1] for j in range(7)]
                                 })
        Speed_DF = pd.concat([Speed_DF,Fixed_RH]).astype(float)
    Speed_DF = Speed_DF.reset_index(drop = True)
    return Speed_DF

def Read_Kayelaby_Attenuation():
    """Function which reads attenuation frequency dependant data from Kaye and
    Laby website and store them into a dataframe for comparisons.

    Returns
    -------
    Attenuation_DF : pandas.Dataframe
        Table of values of sound attenuation at varying temperature (from
        0 °C to 30 °C) and RH (from 10% to 90%), expressed in dB/km

    """
    Attenuation_Table = pd.read_html(kayelaby)[9].dropna()
    Attenuation_DF = pd.DataFrame({
                                   'Attenuation (dB/km)' : [],
                                   'Frequency (kHz)' : [],
                                   'Relative Humidity (%)' : []
                                   })
    for x in range(1,10):
        Fixed_RH = pd.DataFrame({
                                 'Attenuation (dB/km)' :
                                     Attenuation_Table[x][2:],
                                 'Frequency (kHz)' :
                                     Attenuation_Table[0][2:],
                                 'Relative Humidity (%)' :
                                     [Attenuation_Table[x][1]
                                      for j in range(21)]
                                 })
        Attenuation_DF = pd.concat([Attenuation_DF,Fixed_RH]).astype(float)
        Attenuation_DF = Attenuation_DF.reset_index(drop = True)
        return Attenuation_DF
