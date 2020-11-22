# Temperat-here
A python project for Unibo's course in Software and Computing in Applied Physics 
This project aims to study how sound's speed c depends on the temperature, humidity and pressure of air, with a focus on the propagation of uncertainties; taking as main reference the article "The variation of the specific heat ratio and the speed of sound in air with temperature, pressure, humidity, and CO2 concentration", by O.Cramer. 
The final goal is to simulate an experimental apparatus which measures c and employs it to evaluate T. 


# What: the theory

**Speed of sound**
The equation of state for a real gas with reference to density ρ, truncated to the second virial term, appears as:

<a href="https://www.codecogs.com/eqnedit.php?latex=P&space;=&space;R&space;T&space;\rho&space;&plus;&space;R&space;T&space;B(T)&space;\rho^2" target="_blank"><img src="https://latex.codecogs.com/svg.latex?P&space;=&space;R&space;T&space;\rho&space;&plus;&space;R&space;T&space;B(T)&space;\rho^2" title="P = R T \rho + R T B(T) \rho^2" /></a>


From this we can find heat capacity at constant P and at constant V:


<a href="https://www.codecogs.com/eqnedit.php?latex=c_p^1&space;(T)&space;=&space;c_p^0&space;-&space;\frac{R}{M}&space;\frac{p}{RT}&space;(T^2&space;\frac{d^2&space;B&space;(T)&space;}{d&space;T^2})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c_p^1&space;(T)&space;=&space;c_p^0&space;-&space;\frac{R}{M}&space;\frac{p}{RT}&space;(T^2&space;\frac{d^2&space;B&space;(T)&space;}{d&space;T^2})" title="c_p^1 (T) = c_p^0 - \frac{R}{M} \frac{p}{RT} (T^2 \frac{d^2 B (T) }{d T^2})" /></a>


<a href="https://www.codecogs.com/eqnedit.php?latex=c_v^1(T)&space;=&space;c_p^1(T)&space;-&space;\frac{R}{M}[1&plus;\frac{2p}{RT}(T&space;\frac{d&space;B(T)}{d&space;T})]" target="_blank"><img src="https://latex.codecogs.com/svg.latex?c_v^1(T)&space;=&space;c_p^1(T)&space;-&space;\frac{R}{M}[1&plus;\frac{2p}{RT}(T&space;\frac{d&space;B(T)}{d&space;T})]" title="c_v^1(T) = c_p^1(T) - \frac{R}{M}[1+\frac{2p}{RT}(T \frac{d B(T)}{d T})]" /></a>

where density has been substituted (valid in low density-high volume regime) according to:

<a href="https://www.codecogs.com/eqnedit.php?latex=\rho&space;=&space;\frac{1}{V}&space;=&space;\frac{P}{R&space;T}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\rho&space;=&space;\frac{1}{V}&space;=&space;\frac{P}{R&space;T}" title="\rho = \frac{1}{V} = \frac{P}{R T}" /></a>

c at 0 frequency is defined as:


<a href="https://www.codecogs.com/eqnedit.php?latex=c_0^2&space;=&space;\gamma(\frac{\delta&space;p}{\delta&space;\rho})&space;=&space;\gamma&space;(T)&space;\frac{RT}{M}(1&space;&plus;&space;\frac{2&space;p&space;B(T)}{RT})" target="_blank"><img src="https://latex.codecogs.com/svg.latex?c_0^2&space;=&space;\gamma(\frac{\delta&space;p}{\delta&space;\rho})&space;=&space;\gamma&space;(T)&space;\frac{RT}{M}(1&space;&plus;&space;\frac{2&space;p&space;B(T)}{RT})" title="c_0^2 = \gamma(\frac{\delta p}{\delta \rho}) = \gamma (T) \frac{RT}{M}(1 + \frac{2 p B(T)}{RT})" /></a>

while the dispersive effects produce a Δc frequency dependent, bringing to:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{1}{c_\phi}&space;=&space;\frac{1}{c_0}&space;-&space;\frac{\alpha_{\nu&space;N}}{2&space;\pi&space;f_{r&space;N}}&space;&plus;&space;\frac{\alpha_{\nu&space;O}}{2&space;\pi&space;f_{r&space;O}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{1}{c_\phi}&space;=&space;\frac{1}{c_0}&space;-&space;\frac{\alpha_{\nu&space;N}}{2&space;\pi&space;f_{r&space;N}}&space;&plus;&space;\frac{\alpha_{\nu&space;O}}{2&space;\pi&space;f_{r&space;O}}" title="\frac{1}{c_\phi} = \frac{1}{c_0} - \frac{\alpha_{\nu N}}{2 \pi f_{r N}} + \frac{\alpha_{\nu O}}{2 \pi f_{r O}}" /></a>

The main difficulties in c calculations are related to the virial coefficient.

**Virial Coefficient**
B(T) is the contribute to Pressure due to 2-particle interactions:

<a href="https://www.codecogs.com/eqnedit.php?latex=B(T)=\int_0^\infty(e^{\frac{-\Phi(r)}{kT}}-1)dr" target="_blank"><img src="https://latex.codecogs.com/gif.latex?B(T)=\int_0^\infty(e^{\frac{-\Phi(r)}{kT}}-1)dr" title="B(T)=\int_0^\infty(e^{\frac{-\Phi(r)}{kT}}-1)dr" /></a>

The equation is complicated by the necessity of modeling the intermolecular force Φ(r_{12}), so B is usually measured through experimental meaning and tabulated. There are many fitting equations to describe B as a function of T, and the code has four functions built-in:

- Lennard: using a Lennard-Jones-like potential to describe intermolecular force, it fits B faithfully to the analytical description.

- Exponential: taken by ...
<a href="https://www.codecogs.com/eqnedit.php?latex=B(T)&space;=&space;a&space;-&space;b&space;e^{\frac{c}{T}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?B(T)&space;=&space;a&space;-&space;b&space;e^{\frac{c}{T}}" title="B(T) = a - b e^{\frac{c}{T}}" /></a>

- Hyland: suggested by Hyland in "" to determine 
<a href="https://www.codecogs.com/eqnedit.php?latex=B(T)&space;=&space;a&space;-&space;\frac{b}{T}&space;10^{\frac{c}{T^2}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?B(T)&space;=&space;a&space;-&space;\frac{b}{T}&space;10^{\frac{c}{T^2}}" title="B(T) = a - \frac{b}{T} 10^{\frac{c}{T^2}}" /></a>

- Simple parabola
<a href="https://www.codecogs.com/eqnedit.php?latex=B(T)&space;=&space;a&space;T^2&space;&plus;&space;b&space;T&space;&plus;&space;c" target="_blank"><img src="https://latex.codecogs.com/svg.latex?B(T)&space;=&space;a&space;T^2&space;&plus;&space;b&space;T&space;&plus;&space;c" title="B(T) = a T^2 + b T + c" /></a>



The total second virial coefficient is given by:

<a href="https://www.codecogs.com/eqnedit.php?latex=B_{mix}&space;=&space;x_a^2&space;B_{aa}&space;&plus;&space;x_c^2&space;B_{cc}&space;&plus;&space;2&space;x_a&space;x_w&space;B_{aw}&space;&plus;&space;x_w^2&space;B_{ww}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?B_{mix}&space;=&space;x_a^2&space;B_{aa}&space;&plus;&space;x_c^2&space;B_{cc}&space;&plus;&space;2&space;x_a&space;x_w&space;B_{aw}&space;&plus;&space;x_w^2&space;B_{ww}" title="B_{mix} = x_a^2 B_{aa} + x_c^2 B_{cc} + 2 x_a x_w B_{aw} + x_w^2 B_{ww}" /></a>


**The simulation**
The basic concept of the experiment is the following:

- a speaker is set in place, and emits a sound wave of a specific frequency against a 'wall' at distance L;

- the 'wall' absorbes a fraction of the wave and reflects the remaining, again in direction of a microphone at distance l<<L from the speaker;

- the microphone registers the reflected sound wave after a time t, discerning it from the background noise.

The speed cphi can then be evaluated as simply as cphi = 2L/t (as long as l<<L, other configurations will be tested).


# How: the code
Class Environment has instances which resembles a space where temperature, humidity and pressure are homogeneous and constant. It should:
    -read experimental datas, and process them in order to estimate the dependance c(T,H,P)
    -study how uncertainties propagate, handling them with package 'uncertainties' https://pythonhosted.org/uncertainties/
    -simulate how a real ambience would affect such a measurements, including considerations regarding frequency, absorption and sound analysis with package 'pyroomacoustics' https://github.com/LCAV/pyroomacoustics
    -combine the previous considerations

# Who: contributions?
This project will be used to test my capabilites in writing python, both in the form and in the content: this implicates that the work should be done only by myself. Nonetheless, the project can be followed by anyone in the making, and after the evaluation it will be opened to community's contributions.
 
