# Comp-labs-U-235-sim
Monte Carlo simulations of various topologies formed of pure U-235.

This project was focussed on establishing the critical dimensions of pure Uranium 235 in various topologies (not very realistic, I know). 

Gradually being updated to allow for `multiprocessing` support, will utilise `cpu_count() - 2` processes to give your CPU some room to do other things.

An example output of the `run_sphere.py` file has been included: `./length_oscillation_multi_plot_sphere_2.png`, the lines are individual Monte Carlo simulations that are being run in such a way as to converge on the critical mass and by extension, critical mass of a U-235 sphere. Note the fuzziness in the lines' final destinations, this is where the randomness required for MC simulations manifests. See the commit details of `./length_oscillation_multi_plot_sphere_2.png` for mean, range, and stdev associated with the critical radius.
