# ad-scal-heterogeneous

This repository contains the open-source MATLAB code for the research paper "On the Interpretation of Unsteady State Experiments in Heterogeneous Rock by Stochastic Methods" by Omidreza Amrollahinasab, Boris Jammernegg, Siroos Azizmohammadi, and Holger Ott. The code, built upon the MATLAB Reservoir Simulation Toolbox (MRST), simulates immiscible displacement in heterogeneous porous media. It allows for the interpretation of unsteady-state core flooding experiments, considering rock heterogeneity.

## How to Cite

If you use this code in your research, please cite the following paper:

Amrollahinasab, O., Jammernegg, B., Azizmohammadi, S., & Ott, H. (2025). On the Interpretation of Unsteady State Experiments in Heterogeneous Rock by Stochastic Methods. *InterPore Journal, 2*(2), IPJ040625-5. https://doi.org/10.69631/ipj.v2i2nr44

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Getting Started

To use this code, you need a working installation of MATLAB and the MATLAB Reservoir Simulation Toolbox (MRST).

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/omidreza-amrollahi/ad-scal-heterogeneous.git
    ```
2.  **Add the module to MRST:**
    Run the `startup_user.m` script to add the `ad-scal-heterogeneous` module to your MRST path.
3.  **Run an example:**
    Open the `examples/Example_run.m` file in MATLAB to see an example of how to configure and run a simulation. This example demonstrates a 3D simulation of a heterogeneous core.

## Repository Structure

The repository is organized as follows:

* `docs/`: Contains user manuals, lists of examples, and templates for history matching.
* `examples/`: Includes settings files and required data for running example simulations.
* `models/`: Contains the main modules for running the simulations.
* `utils/`: Includes utility modules for outputting and plotting results, as well as input/output functions.
* `startup_user.m`: Script to load the package into MRST.

## Examples

The `examples` directory contains various files to run simulations, including:
* `Example_run.m`: A script to run a 3D heterogeneous core simulation.
* `SettingsDecaneBrine3D.txt`: An example settings file.
* Data files such as `pressure_CO2_brine.txt` and `pressure_decane_brine.txt` which contain experimental data.

For a comprehensive list of examples and keywords for the settings files, please refer to the documents in the `docs` folder.
