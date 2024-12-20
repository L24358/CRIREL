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
        <li><a href="#simulation-tool-flysim">Simulation Tool: Flysim</a></li>
        <li><a href="#neural-circuit-and-simulation-protocol">Neural Circuit and Simulation Protocol</a></li>
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

### Simulation Tool: Flysim

Flysim is an open-sourced neural network simulator that runs on linux. It takes two files as input: (1) a configuration file (``.conf``), where the connectome of the neural circuit is specified, and (2) a protocol file (``.pro``), where the stimulus given is specified. Further specifications, such as neuonal models or numerical methods, can be found by ``./flysim.out -h``. In this project, we use:\
``./flysim.out -pro <protocol-file> -conf <configuration-file> -s moderate - nmodel LIF``.

Link to its published paper: https://www.frontiersin.org/10.3389/conf.fninf.2014.18.00043/event_abstract \
(Note: currently, only version 6 is open-sourced. The execution file in this folder is version 7.21.)

### Neural Circuit and Simulation Protocol

In this project, generation of the .conf file can be done through the ``exec_conf`` function in ``code/dynalysis/gen_conf.py`` by specifying the neural circuit ID (and other optional adjustments). Generation of the .pro file is similarly executed through ``exec_pro`` in ``code/dynalysis/gen_pro.py``. For detailed description of the parameters involved, please see the documentation for Flysim.

The scripts that call the execution functions for each of the results (and thus specifies the parameters) are given in the ``code/fig*/`` folders. Raw simulation outputs are provided in ``results/``. 

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/L24358/CRIREL/issues) for a list of proposed features (and known issues).


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

Belle (Pei-Hsien) Liu - belleliu@uw.edu

Project Link: [https://github.com/L24358/CRIREL](https://github.com/L24358/CRIREL)

