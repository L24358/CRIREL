<!--Source code: https://github.com/othneildrew/Best-README-Template/edit/master/README.md -->

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/L24358/CRIREL">
    <img src="https://github.com/L24358/CRIREL/blob/main/graphs/CRIREL.PNG" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">CRIREL</h3>

  <p align="center">
    Code for <strong><em>Augmenting Flexibility: Mutual Inhibition Between Inhibitory Neurons Expands Functional Diversity</em></strong> by Belle (Pei-Hsien) Liu, Alexander James White, and Chung-Chuan Lo. <br/> bioRxiv link: https://www.biorxiv.org/content/10.1101/2020.11.08.371179v1
    <br />
    <a href="https://github.com/L24358/CRIREL"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/L24358/CRIREL">View Demo</a>
    ·
    <a href="https://github.com/L24358/CRIREL/issues">Report Bug</a>
    ·
    <a href="https://github.com/L24358/CRIREL/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li>
      <a href="#usage">Usage</a>
      <ul>
        <li><a href="#flysim">flysim</a></li>
        <li><a href="#dynalysis">dynalysis</a></li>
        <li><a href="#equilibrium_points">equilibrium_points</a></li>
        <li><a href="#CPG">CPG</a></li>
        <li><a href="#bistable_decision">bistable_decision</a></li>
        <li><a href="#functions">functions</a></li>
        <li><a href="#large_network">large_network</a></li>
      </ul>
    </li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

In the face of a complex environment, it is necessary for neural circuits to develop a “Swiss army knife” toolkit, in which neural networks flexibly switch between different functionalities to meet the demands placed on them. Here, we develop a theoretical framework for how recurrent circuits give rise to this flexibility. Moreover, mutual inhibition, an abundant yet less studied structure, is critical in expanding the functionality of the network, far beyond what feedback inhibition alone can accomplish. Using dynamical systems theory, we show that mutual inhibition doubles the number of cusp bifurcations in small neural circuits. As a concrete example, we build a class of functional motifs we call Coupled Recurrent Inhibitory and Recurrent Excitatory Loops (CRIRELs). These CRIRELs have the advantage of being multi-functional, performing a plethora of functions, including decisions, memory, toggle, etc. Finally, we demonstrate how this trend holds for larger networks.

<!-- USAGE EXAMPLES -->
## Usage

### flysim

Flysim is an open-sourced neural network simulator that runs on linux. It takes two files as input: (1) a configuration file, where the connectome of the neural network is specified, and (2) a protocol file, where the stimulus given is specified. Further specifications, such as neuonal models or numerical methods, can be found by ``./flysim.out -h``. In this project, we use:\
``./flysim.out -pro <protocol-file> -conf <configuration-file> -s moderate - nmodel LIF``

Link to its published paper: https://www.frontiersin.org/10.3389/conf.fninf.2014.18.00043/event_abstract \
(Note: currently, only version 6 is open-sourced. The execution file in this folder is version 7.21.)

Other files in this folder are shell scripts for executing ``flysim.out`` recursively.\
``ss_flysim_iter=0``: runs flysim in current directory.\
``ss_flysim_iter=1``: runs flysim in all first-level sub-directories.\
``ss_flysim_iter=2``: runs flysim in all second-level sub-directories.

### dynalysis

Contains useful modules for simulation and analysis.

``gen_pro``: generates protocol file.\
``gen_conf``: generates configuration files.\
``classes``: includes some commonly used classes. In particular, class ``motif`` generates small neuronal circuits when the ID is given. See how to specify circuit ID here:

### equilibrium_points

To access the complexity of the dynamical system, we counted the number of equilibrium points across a wide range of parameters. For a given neural circuit, it iterates through the parameter space, and for each parameter set it explores the phase space by stimulating the circuit differently in time. 

To generate the .pro (protocol file) and .conf (configuration file) files required for simulation, adjust the circuit ID in ``mkfiles.py``, and run the file.\
Go to the generated file, and run ``./ss_flysim.sh``.\
To obtain the analysis of the results, run ``python analysis.py``.

### CPG

To access the complexity of the dynamical system, we counted the number of central pattern generators (CPG) across a wide range of parameters. For a given neural circuit, it iterates through the parameter space, and determines whether a parameter set is capable of oscillation by analyzing its inter-spike interval (ISI).

### bistable_decision

Determines whether a circuit is capable of performing the functions switch (called bistable here) and decision-making.

### functions

Simulates the functions of a CRIREL circuit.

### large_network

Generates a neural network containing 100 neurons.

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/L24358/CRIREL/issues) for a list of proposed features (and known issues).


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

Belle (Pei-Hsien) Liu - belle.l24358@gmail.com

Project Link: [https://github.com/L24358/CRIREL](https://github.com/L24358/CRIREL)

