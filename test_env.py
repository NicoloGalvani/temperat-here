"""Tests for environment class predicitions."""

from hypothesis import given
import hypothesis.strategies as st
import pandas as pd
from env import Environment

KL_URL = ('https://web.archive.org/web/20190508003406/http://www.kayelaby.npl'
          +'.co.uk/general_physics/2_4/2_4_1.html#speed_of_sound_in_air')
PRECISION = 1E-6

def read_kayelaby_speed():
    """Function which reads speed of sound data from Kaye and Laby website and
    store them into a dataframe for comparisons.

    Returns
    -------
    speed_df : pandas.Dataframe
        table of values of speed of sound at varying temperature (from
        0 °C to 30 °C) and rh (from 10% to 90%), expressed in m/s

    """
    path = 'Data/KL_speed.csv'
    try:
        speed_df = pd.read_csv(path)
        return speed_df
    except FileNotFoundError:
        speed_table = pd.read_html(KL_URL)[10].dropna()
        speed_df = pd.DataFrame({
                                 'Speed (m/s)' : [],
                                 'Temperature (°C)' : [],
                                 'Relative Humidity (%)' : []
                                 })
        for col in range(1,10):
            fixed_rh = pd.DataFrame({
                                    'Speed (m/s)' : speed_table[col][2:],
                                    'Temperature (°C)' : speed_table[0][2:],
                                    'Relative Humidity (%)' :
                                        [speed_table[col][1] for j in range(7)]
                                     })
            speed_df = pd.concat([speed_df, fixed_rh]).astype(float)
        speed_df = speed_df.reset_index(drop = True)
        speed_df.to_csv(path)
        return speed_df

def read_kayelaby_attenuation():
    """Function which reads attenuation frequency dependant data from Kaye and
    Laby website and store them into a dataframe for comparisons.

    Returns
    -------
    attenuation_df : pandas.Dataframe
        table of values of sound attenuation at varying temperature (from
        0 °C to 30 °C) and rh (from 10% to 90%), expressed in dB/km

    """
    path = 'Data/KL_attenuation.csv'
    try:
        attenuation_df = pd.read_csv(path)
        return attenuation_df
    except FileNotFoundError:
        attenuation_table = pd.read_html(KL_URL)[9].dropna()
        attenuation_df = pd.DataFrame({
                                       'Attenuation (dB/km)' : [],
                                       'Frequency (kHz)' : [],
                                       'Relative Humidity (%)' : []
                                       })
        for col in range(1,10):
            fixed_rh = pd.DataFrame({
                                     'Attenuation (dB/km)' :
                                         attenuation_table[col][2:],
                                     'Frequency (kHz)' :
                                         attenuation_table[col][2:],
                                     'Relative Humidity (%)' :
                                         [attenuation_table[col][1]
                                          for j in range(21)]
                                     })
            attenuation_df = pd.concat([attenuation_df, fixed_rh]).astype(float)
            attenuation_df = attenuation_df.reset_index(drop = True)
            return attenuation_df

def test_air_molar_mass_compatible():
    """Check: molar mass compatible with exp. data """
    room = Environment()
    assert abs(room.molar_mass()-0.02897) < PRECISION

@given(st.floats(250.,330.), st.floats(0,100))
def test_xw_in_appropriate_range(temp, rel_hum):
    """Check: xw produced in the range of exp. observations """
    xw_min = 0.
    xw_max = 0.19
    room = Environment(t_input = temp, h_input = rel_hum)
    assert room.rh_to_xw() >= xw_min
    assert room.rh_to_xw() <= xw_max

@given(st.floats(250.,330.), st.floats(0,100), st.floats(101325,111325))
def test_sound_speed_positive(temp, rel_hum, pressure):
    """Check: c0 > 0 """
    room = Environment(temp, rel_hum, pressure)
    assert room.sound_speed() > 0

def test_sound_speed_compatible_with_data():
    """Check: c0 produced in the range of exp. observations """
    precision = 3 #will be reduced when a true sound_speed function will be implemented
    data = read_kayelaby_speed()
    for i in range(len(data)):
        temp_k = data['Temperature (°C)'][i] + 273.15
        rel_hum = data['Relative Humidity (%)'][i]
        room = Environment(temp_k,rel_hum)
        assert abs(room.sound_speed() - data['speed (m/s)'][i]) < precision
