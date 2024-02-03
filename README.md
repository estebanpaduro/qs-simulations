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

* `FHN_exp1.mmt`: MYOKIT source for the experiment 1 for the FitzHigh-Nagumo model.
* `FHN_exp2.mmt`: MYOKIT source for the experiment 2 for the FitzHugh-Nagumo model.
* `FHN_generate_data_exp1.py`: Source files to generate the data for the FHN model for experiment 1.
* `FHN_generate_data_exp2.py`: Source files to generate the data for the FHN model for experiment 2.
* `FHN_reports_from_data_exp1.py`: Generate figures from the csv files generated by `FHN_generate_data_exp1.py`.
* `FHN_reports_from_data_exp2.py`: Generate figures from the csv files generated by `FHN_generate_data_exp2.py`.
* `main_data.py`: Generates data by combining information from `FHN_generate_data_exp1.py` and `FHN_generate_data_exp2.py`.
* `main_figures.py`: Generates figures this article by combining information from `FHN_reports_from_data_exp1.py` and `FHN_reports_from_data_exp2.py`.


## References

[1] Cerpa E., Corrales N., Courdurier M., Medina L., Paduro E. (2023) The impact of high frequency-based stability on the onset of action potentials in neuron models.

## License

[MIT](LICENSE)
