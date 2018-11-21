The code samples configurations for a single-particle harmonic oscillator using a standard Metropolis Monte Carlo scheme.
After the sampling is complete a histogram of the sampled configurations is produced and it is compared with the analytic solution and with the shape of the potential.
The "dynamics" in configuration space and energy space is also shown at the end of the execution for the first 100 steps (editable via #define instance).
The code also prints the acceptance ratio, i.e. the fraction of the accepted moves compared to the proposed ones. 

It is possible to modify from the input file (mc_sampling.in) the number of MC iterations, the maximum displacement for a single move (which is what affects the acceptance ratio) and the initial position of the particle.
It is also possible to modify the potential parameters of the harmonic oscillator and the inverse temperature of the system (this has effects on the shape of the potential through the parameters M and W and on the distribution through the inverse temperature BETA).
The parameters for the histogram are also editable. In particular it is possible to define the domain of the histogram through the parameter L (the histogram will be from -.5L and +.5L) and the number of bins through the parameter NBINS.

It could be interested to investigate:
- how the number of MC iteration affects the sampled distribution.
- how the maximum displacement affects the acceptance ratio and the sampled distribution for a fixed number of MC iterations.
- why, for an initial state not at the minimum of the potential, it is important to throw away the equilibration steps. What happens with an initial position far away from the minimum and with a small maximum displacement?

It could be also nice to implement a different potential or even do the same code in a 2-dimension configuration space.

Histogram parameters can be left untouched, as long as one does not want to investigate the sampling of the tails of the distribution for potential parameters which deviate a lot from unity.

It requires GNUplot for plot generation.