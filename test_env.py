"""
Test file for the class Environment defined in env.py

"""

from hypothesis import given
import hypothesis.strategies as st
import pandas as pd
from env import Environment


KAYE = ('https://web.archive.org/web/20190508003406/http://www.kayelaby.' +
        'npl.co.uk/general_physics/2_4/2_4_1.html#speed_of_sound_in_air')
PRECISION = 0.001
ERROR = 4#will be reduced when a true sound_speed function will be implemented


def read_kayelaby_speed():
    """Function which reads speed of sound data from Kaye and Laby website and
    store them into a dataframe for comparisons.

    Returns
    -------
    speed_df : pandas.Dataframe
        Table of values of speed of sound at varying temperature (from
        0 °C to 30 °C) and rel_humidity (from 10% to 90%), expressed in m/s

    """
    speed_table = pd.read_html(KAYE)[10].dropna()
    speed_df = pd.DataFrame({
                             'Speed (m/s)' : [],
                             'Temperature (°C)' : [],
                             'Relative Humidity (%)' : []
                             })
    for column in range(1,10):
        fixed_rel_humidity = pd.DataFrame({
                                'Speed (m/s)' : speed_table[column][2:],
                                'Temperature (°C)' : speed_table[0][2:],
                                'Relative Humidity (%)' :
                                    [speed_table[column][1] for j in range(7)]
                                 })
        speed_df = pd.concat([speed_df,fixed_rel_humidity]).astype(float)
    speed_df = speed_df.reset_index(drop = True)
    return speed_df

def read_kayelaby_attenuation():
    """Function which reads attenuation frequency dependant data from Kaye and
    Laby website and store them into a dataframe for comparisons.

    Returns
    -------
    attenuation_df : pandas.Dataframe
        Table of values of sound attenuation at varying temperature (from
        0 °C to 30 °C) and rel_humidity (from 10% to 90%), expressed in dB/km

    """
    attenuation_table = pd.read_html(KAYE)[9].dropna()
    attenuation_df = pd.DataFrame({
                                   'attenuation (dB/km)' : [],
                                   'Frequency (kHz)' : [],
                                   'Relative Humidity (%)' : []
                                   })
    for column in range(1,10):
        fixed_rel_humidity = pd.DataFrame({
                                 'attenuation (dB/km)' :
                                     attenuation_table[column][2:],
                                 'Frequency (kHz)' :
                                     attenuation_table[0][2:],
                                 'Relative Humidity (%)' :
                                     [attenuation_table[column][1]
                                      for j in range(21)]
                                 })
        attenuation_df = pd.concat([attenuation_df,
                                    fixed_rel_humidity]).astype(float)
        attenuation_df = attenuation_df.reset_index(drop = True)
        return attenuation_df

def test_air_molar_mass_compatible(): # pylint: disable=missing-function-docstring
    room = Environment()
    experimental_molar_mass = 0.02897
    assert abs(room.molar_mass - experimental_molar_mass) < PRECISION

@given(st.floats(250.,330.), st.floats(0,100))
def test_xw_in_appropriate_range(temperature, rel_humidity): # pylint: disable=missing-function-docstring
    xw_min = 0.
    xw_max = 0.19
    room = Environment(t_input = temperature, h_input = rel_humidity)
    assert room.rh_to_xw() >= xw_min
    assert room.rh_to_xw() <= xw_max

@given(st.floats(250.,330.), st.floats(0,100), st.floats(101325,111325))
def test_sound_speed_positive(temperature, rel_humidity, pressure): # pylint: disable=missing-function-docstring
    room = Environment(temperature, rel_humidity, pressure)
    assert room.sound_speed() > 0

def test_sound_speed_compatible_with_data(): # pylint: disable=missing-function-docstring
    data = read_kayelaby_speed()
    for i in range(len(data)):
        temperature_k = data['Temperature (°C)'][i] + 273.15
        rel_humidity = data['Relative Humidity (%)'][i]
        room = Environment(temperature_k,rel_humidity)
        assert abs(room.sound_speed() - data['Speed (m/s)'][i]) < ERROR
