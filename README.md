# Generate figures used in the paper: The impact of high frequency-based stability on the onset of action potentials in neuron models
Python and MYOKIT files used to generate the graphs in the article.: A partially averaged system to model neuron responses to interferential current stimulation
Authors: Cerpa E., Corrales N., Courdurier M., Medina L., Paduro E.

## Download 
From a command line:
```
git clone https://github.com/estebanpaduro/qs-simulations
```

Current version only works on Linux, windows MYOKIT release seems to be bugged for now.
for MYOKIT installation see: http://myokit.org/linux


## Contents
The scripts are organized in the following folders:

* `numerical_experiments.ipynb`: Jupiter Notebook with the simulation files.
* `FHN.mmt`: MYOKIT source for the experiment for the FitzHigh-Nagumo model.
* `HH.mmt`: MYOKIT source for the experiment for the Hodgkin-Huxley model.

## References

[1] Cerpa E., Corrales N., Courdurier M., Medina L., Paduro E. (2023) The impact of high frequency-based stability on the onset of action potentials in neuron models.

## License

[MIT](LICENSE)
