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

TEMP_0 = 273.15
ATM_P = 101325
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
    T_input: Float, optional
        Temperature inside the environment, expressed in K
    H_input: Float, optional
        Humidity inside the environment, expressed in %
    P_input: Float, optional
        Pressure inside the environment, expressed in Pa
    """

    T_input: float = field(default = TEMP_0 + 25, metadata={'unit' : 'K'})
    H_input: float = field(default = 0.0, metadata={'unit' : '%'})
    P_input: float = field(default = ATM_P, metadata={'unit' : 'Pa'})
    Molar_Fraction_path: Path = field(default = 'Data/Molar_Fraction.txt')
    Air_Data_path: Path = field(default = 'Data/Dry_Air.txt')
    Water_Data_path: Path = field(default = 'Data/H2O.txt')
    CO2_Data_path: Path = field(default = 'Data/CO2.txt')
    Molar_Fraction = pd.DataFrame()
    Air_Data = pd.DataFrame()
    Water_Data = pd.DataFrame()
    CO2_Data = pd.DataFrame()
    B_values = pd.DataFrame()
    B_covariances = pd.DataFrame()
    Speed_Data = pd.DataFrame()
    Attenuation_Data = pd.DataFrame()
    Molar_Mass = float

    ##Possibility to insert also the link to pyroomsound?

    def __post_init__(self):
        """
        When the first instance of class Environment is created:
            -it will store data from the paths provided into databases
            -it will fit those data to obtain the right parameters of Bi(T)
            -it will read and store data from Kaye and Laby website archive

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
                Environment.Water_Data = pd.read_csv(self.Water_Data_path, sep=' ')
                Environment.CO2_Data = pd.read_csv(self.CO2_Data_path, sep=' ')
                (Environment.B_values,
                 Environment.B_covariances) = Environment.B_fit()
                (Environment.Speed_Data,
                 Environment.Attenuation_Data) = Environment.Read_Kayelaby()
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
        """Simple parabolic function, good enough to fit fast Baa and Bcc."""
        return A + B*T + C*T**2

    def Exponential(T: float, A: float, B: float, C: float):
        """Exponential function, useful to evaluate Baa or its component,
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
        """Approximated function to describe Bww, taken from Hyland
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
        """Approximated function to evaluate Baw(T), taken from Hyland"""
        return (
                36.98928
                - 0.331705*T
                + 0.139035E-2*T**2
                - 0.574154E-5*T**3
                + 0.326513E-7*T**4
                - 0.142805E-9*T**5
                )

    def B_fit(draw: bool = False):
        """Function which is executed as an instance of the class is created:
        it reads data from databases and fits them with the appropriate
        functions, to store the optimized parameters.

        Parameters
        ----------
        draw : bool, optional
            If True, the function will draw plots of the fits, default = False.

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
                Environment.CO2_Data]
        functions = [Environment.Exponential,
                     Environment.Hyland,
                     Environment.Hyland]
        p0 = [[1, 1, 1], [21, 147, 100], [33.97, 55306, 72000]]
        title = ['Air', 'Water', 'CO2']
        Optimized_parameters = []
        Optimized_covariances = []
        for i,D in enumerate(Data):
            function = functions[i]
            x = D['T,K']
            y = D['B,cm^3/mol']
            popt, pcov = curve_fit(function, x, y, p0[i])
            Optimized_parameters.append(popt)
            Optimized_covariances.append(pcov)
            if draw:
                new_x = np.arange(273, 313, 0.5)
                new_y = function(new_x, popt[0], popt[1], popt[2])
                fig,ax = plt.subplots(nrows = 1)
                ax.tick_params(direction = 'in')
                fig.set_figheight(6)
                fig.set_figwidth(6)
                sns.scatterplot(x, y, ax = ax, color = 'blue', label = 'data')
                sns.lineplot(new_x, new_y, ax = ax,
                            color = 'orange', label = 'fit')
                plt.title(title[i])
                ax.set_ylabel('B(T) (cm^3/mol)')
                ax.set_xlabel('Temperature (K)')
        return Optimized_parameters, Optimized_covariances

    def Read_Kayelaby():
        """Function which is called as an instance of the class is created:
        it reads data from Kaye and Laby website and store them into
        dataframes common to all instances.


        Returns
        -------
        Speed_DF : pandas.Dataframe
            Table of values of speed of sound at varying temperature (from
            0 °C to 30 °C) and RH (from 10% to 90%), expressed in m/s
        Attenuation_DF : pandas.Dataframe
            Table of values of sound attenuation at varying temperature (from
            0 °C to 30 °C) and RH (from 10% to 90%), expressed in dB/km

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
        Speed_DF = Speed_DF.reset_index(drop = True)
        Attenuation_DF = Attenuation_DF.reset_index(drop = True)
        return Speed_DF, Attenuation_DF

    def Set_Molar_Mass(self):
        """Evaluates molar mass of the air components from Molar Fraction
        Database, taking into account the most relevant components.
        Called as an instance of the class is created.

        Returns
        -------
        Mass : float
            Combined molar mass of dry air CO2 free, water vapor and CO2
        """
        Molecules = ['N2', 'O2', 'Ar', 'Ne']
        EMF = Environment.Molar_Fraction
        Molar_Masses_xi = np.array([float(EMF[EMF['Constituent'] == m]['xiMi'])
                                    for m in Molecules])
        M = [float(EMF[EMF['Constituent'] == m]['Mi']) for m in Molecules]
        Air_Molar_Mass = np.sum(Molar_Masses_xi)/np.sum(M)
        Water_Molar_Mass = float(EMF[EMF['Constituent'] == 'H2O']['Mi'])
        CO2_Molar_Mass = float(EMF[EMF['Constituent'] == 'CO2']['Mi'])
        xcc = float(EMF[EMF['Constituent'] == 'CO2']['xi'])
        xww = self.RH_to_Xw()
        xaa = (1-xcc-xww)
        Mass = Air_Molar_Mass*xaa + Water_Molar_Mass*xww + CO2_Molar_Mass*xcc
        return Mass

    def RH_to_Xw(self):
        """Conversion from relative humidity to water molar fraction,
        taken from 'A simple expression for the saturation vapour pressure
        of water in the range −50 to 140°C' J M Richards
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

        """
        Molecules = ['N2', 'O2', 'Ar', 'Ne']
        Molar_Masses = np.array([float(Environment.Molar_Fraction
                                        [Environment.Molar_Fraction
                                        ['Constituent'] == m]['xiMi']
                                        ) for m in Molecules
                                ])
        return np.sum(Molar_Masses)

    def Sound_Speed(self):
        """Mockup function, used now for testing, later implemented.
        """
        return 20.5*np.sqrt(self.T_input)
