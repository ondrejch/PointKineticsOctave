# Point Kinetics Equation (PKE) solver for education and research purposes.

Adapted from https://github.com/ondrejch/MSBR-ORNL-4528/tree/master/dynamic_model/Octave

Author: Vikram Singh, viikraam@gmail.com and Ondrej Chvala, ochvala@utk.edu

This program takes an input file './reactivity.dat' with time (in s), reactivity,
and external neutron source as three columns, formatted as
< time reactivity ext-source > per line. 

The code plots for n(t) and C_i(t) for delayed neutron groups i=1,2,...6

