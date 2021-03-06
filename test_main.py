from hypothesis import given
import hypothesis.strategies as st
from main import Environment
from main import Read_Kayelaby_Speed

def test_air_molar_mass_compatible():
    Room = Environment()
    assert abs(Room.Molar_Mass-0.02897) < 0.0001

@given(st.floats(250.,330.), st.floats(0,100))
def test_xw_in_appropriate_range(T, RH):
    Xw_min = 0.
    Xw_max = 0.19
    Room = Environment(T_input = T, H_input = RH)
    assert Room.RH_to_Xw() >= Xw_min
    assert Room.RH_to_Xw() <= Xw_max

@given(st.floats(250.,330.), st.floats(0,100), st.floats(101325,111325))
def test_sound_speed_positive(T, RH, P):
    Room = Environment(T, RH, P)
    assert Room.Sound_Speed() > 0

def test_sound_speed_compatible_with_data():
    Precision = 3 #will be reduced when a true Sound_Speed function will be implemented
    Data = Read_Kayelaby_Speed()
    for i in range(len(Data)):
        T_K = Data['Temperature (°C)'][i] + 273.15
        RH = Data['Relative Humidity (%)'][i]
        room = Environment(T_K,RH)
        assert abs(room.Sound_Speed() - Data['Speed (m/s)'][i]) < Precision
