# Generate figures used in the paper: The impact of high frequency-based stability on the onset of action potentials in neuron models
Python and MYOKIT files used to generate the graphs in the article.: The impact of high frequency-based stability on the onset of action potentials in neuron models
Authors: Cerpa E., Corrales N., Courdurier M., Medina L., Paduro E.




## Download 
From a command line:
```
git clone https://github.com/estebanpaduro/qs-simulations
```

Current version if the scripts for data generation seems to only works on Linux, windows MYOKIT release seems to be bugged for now.
for MYOKIT installation see: http://myokit.org/linux


## Requirements
The following code ran with the following library versions:

matplotlib  ==  3.8.2

Myokit      ==  1.35.4

numpy       ==  1.24.3

pandas      ==  2.1.1

seaborn     ==  0.13.2

scipy       ==  1.11.3


## Contents
The scripts are organized in the following folders:

* `FHN_exp1.mmt`: MYOKIT source for the experiment 1 for the FitzHigh-Nagumo model.
* `FHN_exp2.mmt`: MYOKIT source for the experiment 2 for the FitzHugh-Nagumo model.
* `fig1_example_input`: Generate figure 1 in the paper with an example of the input current under consideration
* `fig2_generate_data.py`: Source files to generate the data for the FHN model for experiment 1.
* `fig2_generate_figures.py`: Generate images in figure 2 from the csv files generated by `fig2_generate_data.py`.
* `fig3_generate_data.py`: Source files to generate the data for the FHN model for experiment 2.
* `fig3_generate_figures.py`: enerate images in figure 2 from the csv files generated by `fig3_generate_data.py`.


## References

[1] Cerpa E., Corrales N., Courdurier M., Medina L., Paduro E. (2023) The impact of high frequency-based stability on the onset of action potentials in neuron models.

## Funding
This work has been partially supported by ANID Millennium Science Initiative Program through Millennium Nucleus for Applied Control and Inverse Problems NCN19-161

## License

[MIT](LICENSE)
