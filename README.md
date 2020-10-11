# Temperat-here
a project for Unibo's course in Software and Computing in Applied Physics 

# What: the theory
This project aims to study how sound's speed c depends on the temperature, humidity and pressure of the gas in which propagates, with a focus on the propagation of uncertainties: the final goal is to understand how much they could affect an evaluation of temperature through experimental measurements of c.
...

# How: the code
Creation of a class, Environment, whose instances resembles a space where temperature, humidity and pressure are homogeneous and constant. It should:
    -read experimental datas, and process them in order to estimate the dependance c(T,H,P)
    -study how uncertainties propagate, handling them with package 'uncertainties' https://pythonhosted.org/uncertainties/
    -simulate how a real ambience would affect such a measurements, including considerations regarding frequency, absorption and sound analysis with package 'pyroomacoustics' https://github.com/LCAV/pyroomacoustics
    -combine the previous considerations

# Who: about contributions
This project will be used to test my capabilites in writing python, both in the form and in the content: this implicates that the work should be done only by myself. Nonetheless, the project can be followed by anyone in the making, and after the evaluation it will be opened to community's contributions.
