"""Tests for measure.py module."""
import collections
from hypothesis import given, settings
import hypothesis.strategies as st
import numpy as np
# from env import Environment
import measure
try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections

@given(st.floats(250.,330.), st.floats(0,100))
def test_positive_speed(temperature, humidity):
    """check: measured speed is positive
    """
    speed = measure.measure(method = 'pyroom')
    assert speed.any()>0
    
@given(st.data())
@settings(deadline=8000)
def test_regression_self_consistent(data):
    """check: regression provided by knn_regressor gives correct results.
    """
    precision = 1E-3
    fingerprint_lenght = 10#data.draw(st.integers(10,25), label='Fingerprint_n')
    n_humidity_samples = data.draw(st.integers(11,41), label='Humidity_n')
    n_temperature_samples = data.draw(st.integers(11,41), label='Temperature_n')
    index = data.draw(st.integers(n_temperature_samples+1,
                                   (n_temperature_samples-1)*n_humidity_samples
                                   ), label = 'Sample_index')
    reference_data = measure.generate_fingerprint(fingerprint_lenght,
                                              n_humidity_samples,
                                              n_temperature_samples)
    target = np.array(reference_data[['Temperature','Humidity']])[index,:]
    test = reference_data.drop(['Temperature','Humidity'], axis=1
                               )[reference_data.index==index]
    fitted_target = measure.knn_regressor(reference_data, test)[0,:]
    assert all((abs(target[i] - fitted_target[i]) < precision
                    for i in range(2)))
