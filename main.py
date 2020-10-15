import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field

TEMP_0 = 273.15
ATM_P = 101325

@dataclass
class Environment ():
    """Instances of this class represents a thermodynamic system
    with homogeneous temperature, humidity and pressure.
    It will later implement compatibility with uncertainties for
    errors evaluation and with pyroomacoustics for geometrical considerations.

    Parameters
    ----------
    Molar_Fraction_path: Path
        Path from which takes experimental data regarding air composition
    Air_Data_path: Path
        Path from which takes experimental data regarding dry air without CO2's
        B parameter at varying T
    Water_Data_path: Path
        Path from which takes experimental data regarding water vapor's
        B parameter at varying T
    CO2_Data_path: Path
        Path from which takes experimental data regarding CO2's
        B parameter at varying T
    T_input: Float
        Temperature inside the environment, expressed in K
    H_input: Float
        Humidity inside the environment, expressed in %
    P_input: Float
        Pressure inside the environment, expressed in Pa
    """

    T_input: float = field(default = TEMP_0 + 25, metadata={'unit' : 'K'})
    H_input: float = field(default=0.0, metadata={'unit' : '%'})
    P_input: float = field(default=ATM_P, metadata={'unit' : 'Pa'})
    Molar_Fraction_path: Path = field(default='**/Data/Molar_Fraction.txt')
    Air_Data_path: Path = field(default='**/Data/Dry_Air.txt')
    Water_Data_path: Path = field(default='**/Data/H2O.txt')
    CO2_Data_path: Path = field(default='**/Data/CO2.txt')
    Molar_Fraction: pd.DataFrame = field(init=False)
    Air_Data: pd.DataFrame = field(init=False)
    Water_Data: pd.DataFrame = field(init=False)
    CO2_Data: pd.DataFrame = field(init=False)


    ##Possibility to insert also the link to pyroomsound?

    def __post_init__(self):
       try:
           self.Molar_Fraction = pd.read_csv(self.Molar_Fraction_path, sep=' ')
           self.Air_Data = pd.read_csv(self.Air_Data_path, sep=' ')
           self.Water_Data = pd.read_csv(self.Water_Data_path, sep=' ')
           self.CO2_Data = pd.read_csv(self.CO2_Data_path, sep=' ')
       except:
           raise FileNotFoundError("Unable to find the reference data")


    def Lennard (T: float or np.ndarray, m: int, e: float, s: float):
        """Ab initio definition of second virial potential B, as the
        integral for r€[0,inf[ of (1-exp(U(r)/KbT)); where U is the
        interatomic potential. Valid description for Dry Air without CO2.
        Taken from Sengers, Klein and Gallagher, (1971)
        Pressure-volume-temperature relationships of gases-virial coefficients.
        U.S. National Bureau of Standards Report,
        U.S. Air Forcc Technical Publication

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

        def integrand(r: float, T, m: int, e: float, s: float):
                amplitude = (m * e) / (m - 6) * (m / 6)**(6 / (m - 6))
                attractive = (s / r)**6
                repulsive = (s / r)**m
                potential = amplitude * (repulsive - attractive)
                return 1 - np.exp(-potential / T)

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

    def parabole(T, A, B, C):
        """Simple parabolic function, good enough to fit fast Baa, Bcc and Bww
        """
        return A + B*T + C*T**2

    def exponential(T, A, B, C):
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
        return A - B * np.exp(-C / T)

    def Hyland(T, A, B, C):
        """Approximated function to describe Bww, taken from Hyland
        DOI: 10.6028/jres.079A.017.


        Parameters
        ----------
        T : Float or numpy.ndarray
            Temperature data
        A : Float
            0 Temperature value of Bww
        B : Float
            Linear decrease Baa(T) dependance
        C : Float
            Temperature lenght-scale of Bww

        Returns
        -------
        Float
            Exponential description of Bww(T)

        """
        return A - B/T * 10**(C / T**2)
