from hypothesis import given
import hypothesis.strategies as st
from main import Environment

def test_air_molar_mass_positive():
    Room = Environment()
    assert Room.Air_Molar_Mass() > 0

def test_xw_lesser_than_one():
    Room = Environment()
    assert Room.RH_to_Xw() >= 0
    assert Room.RH_to_Xw() < 0.1

@given(st.floats(250,330),st.floats(0,100),st.floats(101325,111325))
def test_sound_speed_positive(T,H,P):
    Room = Environment(T,H,P)
    assert Room.Sound_Speed() > 0
