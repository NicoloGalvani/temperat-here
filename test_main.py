"""Tests for main class analysis."""
import collections
from hypothesis import given, settings
import hypothesis.strategies as st
import numpy as np
# from env import Environment
import main
try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections

@given(st.data())
@settings(deadline=8000)
def test_regression_self_consistent(data):
    """check: regression provided by knn_regressor provides correct results.
    """
    precision = 1E-3
    fingerprint_lenght = data.draw(st.integers(5,50), label='Fingerprint_n')
    n_humidity_samples = data.draw(st.integers(11,41), label='Humidity_n')
    n_temperature_samples = data.draw(st.integers(11,41), label='Temperature_n')
    index = data.draw(st.integers(n_temperature_samples+1,
                                   (n_temperature_samples-1)*n_humidity_samples
                                   ), label = 'Sample_index')
    reference_data = main.generate_fingerprint(fingerprint_lenght,
                                              n_humidity_samples,
                                              n_temperature_samples)
    target = np.array(reference_data[['Temperature','Humidity']])[index,:]
    test = reference_data.drop(['Temperature','Humidity'], axis=1
                               )[reference_data.index==index]
    fitted_target = main.knn_regressor(reference_data, test)[0,:]
    assert all((abs(target[i] - fitted_target[i]) < precision
                    for i in range(2)))
