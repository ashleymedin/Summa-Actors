# SUMMA-Actors: Structure for Unifying Multiple Modeling Alternatives with Actors
SUMMA-Actors is a powerful extension of the existing [SUMMA](https://github.com/CH-Earth/summa#readme) hydrological modeling framework, designed to leverage the Actor Model for enhanced scalability and fault-tolerance. SUMMA-Actors is built using the [C++ Actor Framework](https://github.com/actor-framework/actor-framework) and the key highlights include:
  * Scalability: Actors process messages concurrenty, effortlessly scaling to thousands of HRUs/GRUs.
  * Fault-Tolerance: Individual HRUs/GRUs can fail without affecting the rest of the simulation.

## Resources
  * SUMMA-Actors Wiki: https://github.com/uofs-simlab/Summa-Actors/wiki
  * SUMMA Documentation: https://summa.readthedocs.io/en/latest/
  * Laugh-test framework: https://git.cs.usask.ca/numerical_simulations_lab/hydrology/laugh_tests

## Bug Reports, Feature Requests, and Questions
For bug reports, feature requests, and questions, please open an issue on the GitHub at https://github.com/uofs-simlab/Summa-Actors/issues

## Quick Start
SUMMA-Actors seamlessly integrates with two versions of the SUMMA hydrological modeling framework:
  * SUMMA version 4.x.x:
      *  Built on top of the most up-to-data SUMMA codebase: https://github.com/ashleymedin/summa/tree/develop
      * Includes latest features and bug fixes, including the option to use the Sundials numerical solver.
  * SUMMA version 3.x.x:
      * Built on top of the original SUMMA codebase: https://github.com/CH-Earth/summa
      * Suitable for existing projects and familiar workflows.

The build process for both verisions is largely the same, below are the steps common to both versions. For version specific instructions, please refer to the relevant sections below.
### Directory Structure
The directory structure for SUMMA-Actors is as follows:
```
Summa-Actors/
├── bin/
├── build/
│   ├── includes/
│   ├── source/
│   ├── summa/ (Versions Specific - 3.x.x or 4.x.x)
│   ├── v3_build_scripts/
│   └── v4_build_scripts/
├── containers/
│   ├── apptainer.def
│   └── Dockerfile
├── .gitignore
└── README.md
```

### Dependencies
  * g++
  * gfortran
  * [OpenBLAS](https://github.com/xianyi/OpenBLAS)
  * [NetCDF-Fortran](https://github.com/Unidata/netcdf-fortran)
  * NetCDF-C 
  * [C++ Actor Framework (0.18.6)](https://github.com/actor-framework/actor-framework/releases/tag/0.18.6)
  * [Sundials v7.0.0 (Only for SUMMA version 4.x.x)](https://github.com/LLNL/sundials/releases/tag/v7.0.0)

Install most dependencies using your preferred package manager. If you’re using Ubuntu, check our Dockerfile for specific installation examples. 

### Version 4.x.x Build Instructions
  1) git clone https://github.com/uofs-simlab/Summa-Actors.git
  2) cd Summa-Actors/build
  3) git clone -b develop https://github.com/ashleymedin/summa.git
  4) cd v4_build_scripts
  5) Modify the environment variables in the `build_v4_local.sh` script to match your system.
  6) ./build_v4_local.sh

Note: Modify and use the `build_v4_cluster.sh` script for building on the [Digitial Research Alliance of Canada](https://docs.alliancecan.ca/wiki/Getting_started) clusters.

### Version 3.x.x Build Instructions
  1) git clone https://github.com/uofs-simlab/Summa-Actors.git
  2) cd Summa-Actors/build
  3) git clone https://github.com/CH-Earth/summa.git
  4) cd v3_build_scripts
  5) Modify the environment variables in the `build_v3_local.sh` script to match your system.
  6) ./build_v3_local.sh

Note: Modify and use the `build_v3_cluster.sh` script for building on the [Digitial Research Alliance of Canada](https://docs.alliancecan.ca/wiki/Getting_started) clusters.

## Running SUMMA-Actors
Running SUMMA-Actors is similar to running the original version of SUMMA. **Input and configuration files remain identical** alowing exising projects and `fileManager.txt` files to be used seamlessly with SUMMA-Actors. Please refer to the [SUMMA documentation](https://summa.readthedocs.io/en/latest/) regarding input files and simulation configuration. The only difference, if desired, is the option to use a `config.json` file to fine tune how SUMMA-Actors will perform. Please refer to the [relevant section](###Config-File-and-Advanced-Features) for more information on the `config.json` file and the more advanced features of SUMMA-Actors.

Below is the help message for SUMMA-Actors, which provides a brief overview of both the avialable options and the currently unimplemented options.
```   
Usage: summa_actors -m master_file [-g startGRU countGRU] [-c config_file] [-b backup_server] [-s server_mode]
  Available options:
    -m, --master:         Define path/name of master file (can be specified in config)
    -g, --gru:            Run a subset of countGRU GRUs starting from index startGRU 
    -c, --config:         Path name of the Summa-Actors config file (optional but recommended)
    -s, --suffix          Add fileSuffix to the output files
        --gen-config:     Generate a config file
    -b, --backup-server:  Start backup server, requires a server and config_file
        --server-mode:    Enable server mode
    -h, --help:           Print this help message
  Unimplemented Options:
    -n, --newFile         Define frequency [noNewFiles,newFileEveryOct1] of new output files
    -h, --hru             Run a single HRU with index of iHRU
    -r, --restart         Define frequency [y,m,d,e,never] to write restart files
    -p, --progress        Define frequency [m,d,h,never] to print progress
    -v, --version         Display version information of the current build
```

### Example Usage

```bash
./summa_actors -g 1 10 -m /path/to/master_file.txt
```

### Config File and Advanced Features

#### Config File
The `config.json` file is a JSON file that is used to configure SUMMA-Actors. It can be generated by running `./summa_actors --gen-config`, and allows some fine tunning of the SUMMA-Actors program including operating SUMMA-Actors in additional modes. The details of the config file can be found on our wiki page [here](https://github.com/uofs-simlab/Summa-Actors/wiki/Config-File). 

Example usage of the `config.json` file is as follows. Note that the `config.json` file has a field for the `file_master.txt` file, so the `-m` flag is not required when using the `config.json` file.

```bash
./summa_actors -g 1 10 -c /path/to/config.json
```

#### Advanced Features
SUMMA-Actors has additional feature that are not covered in this README. For more information on these features, please refer to the [SUMMA-Actors Advanced Features Wiki Page](https://github.com/uofs-simlab/Summa-Actors/wiki/Advanced-Features). Here is a short summary of some of the optional features:
 
 * Distributed Mode: Run SUMMA-Actors across nodes, or create your own ad-hoc cluster.
 * Data Assimilation Mode: Use SUMMA-Actors to perform data assimilation, restricting all HRUs to complete a timestep before moving to the next.
 * Asynchronous Mode: Default mode of SUMMA-Actors, where HRUs can complete timesteps concurrently and independently.

## Scientific Use:
Please feel free to contribute to our project by submitting pull requests or opening issues. We only ask that if you use SUMMA-Actors that you kindly cite one of our publications:
```bibtex
@article{klenk2024improving,
  title={Improving resource utilization and fault tolerance in large simulations via actors},
  author={Klenk, Kyle and Spiteri, Raymond J},
  journal={Cluster Computing},
  pages={1--18},
  year={2024},
  publisher={Springer}
}

@inproceedings{klenk2024high,
  title={High-Throughput Scientific Computation with Heterogeneous Clusters: A Kitchen-Sink Approach using the Actor Model},
  author={Klenk, Kyle and Moayeri, Mohammad Mahdi and Spiteri, Raymond J},
  booktitle={Proceedings of the 2024 SIAM Conference on Parallel Processing for Scientific Computing (PP)},
  pages={78--89},
  year={2024},
  organization={SIAM}
}
```


## Credits
The initial implementation of SUMMA is credited to the initial publications below. These 
publications can be found in [Water Resources Research](http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1944-7973).

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, L. D. Brekke, J. R. Arnold, D. J. Gochis, R. M. Rasmussen, 2015a: A unified approach for process-based hydrologic modeling: Part 1. Modeling concept. _Water Resources Research_, [doi:10.1002/2015WR017198](http://dx.doi.org/10.1002/2015WR017198).<a id="clark_2015a"></a>

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, D. J. Gochis, R. M. Rasmussen, D. G. Tarboton, V. Mahat, G. N. Flerchinger, D. G. Marks, 2015b: A unified approach for process-based hydrologic modeling: Part 2. Model implementation and case studies. _Water Resources Research_, [doi:10.1002/2015WR017200](http://dx.doi.org/10.1002/2015WR017200).<a id="clark_2015b"></a>

We also credit the original creators of the C++ Actor Framework which allowed us to implement the actor model into SUMMA-Actors. Links to their research work can be found 
below.

 * Charousset, D., Schmidt, T. C., Hiesgen, R., Wählisch, M., 2013: Native actors: 
 a scalable software platform for distributed, heterogeneous environments. _AGERE!_, 
 [doi:10.1145/2541329.2541336](http://dx.doi.org/10.1145/2541329.2541336).

 * Charousset, D., Schmidt, T. C., Hiesgen, R., 2016: Revisiting actor programming in 
 C++. _Computer Languages, Systems & Structures_, [doi:10.1016/j.cl.2016.01.002](http://
 dx.doi.org/10.1016/j.cl.2016.01.002)



