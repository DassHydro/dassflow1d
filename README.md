# DassFlow1D
This git repository contains the development version of DassFlow1D.

## Table of Contents
1. [ General Information. ](#geninfo)
2. [ Repository organisation. ](#reporg)
3. [ Requirements. ](#req)
4. [ Install project. ](#ins)
5. [ Compile and install the adjoint code. ](#adj)
6. [ First steps and Website. ](#next)

<a name="geninfo"></a>
## 1. General information
<strong>WARNING: </strong>For the moment, dassflow1d and other versions work on <strong><ins>Linux only</ins></strong> *(might work on some MacOS installations, no warranty)*.

You can clone this project to your own machine using the following command:
````
git clone http://github.com/DassHydro-dev/dassflow1d
````

<a name="reporg"></a>
## 2. Repository organisation
<ul>
  <li><strong>build</strong>:
  <ul>
    <li><strong>obj</strong>: contains the binaries after installation / compilation. </li>
    <li><strong>api</strong>: contains the application programming interface. </li>
    <li><strong>tap</strong>: contains the Tapenade files used to generate and execute the adjoint code. </li>
  </ul>
  </li>
  <li><strong>cases</strong>: contains up-to-date and validated test cases. </li>
  <li><strong>doc</strong>: contains several documentations for DassFlow1D. </li>
  <li><strong>libs</strong>: contains external libraries that the software uses. </li>
  <li><strong>src</strong>: contains all Fortran source codes. </li>
  <li><strong>unit_tests</strong>: contains unitary tests used to validate DassFlow1D features. </li>
</ul>

<a name="req"></a>
## 3. Requirements
### 3.1 Operating System
The **DassFlow1D** software is dedicated for **Linux or Unix** systems.
### 3.2 Compilers
One of these two compilers must be installed:  
- The GNU Fortran compiler: <code>gfortran</code>
- The INTEL Fortran compiler: <code>ifort</code>
### 3.3 Python environment
- <code>Python</code> >= 3.6
- <code>matplotlib</code> >= 2.0
- <code>numpy</code> >= 1.10
- <code>scipy</code> >= 1.3
- <code>f90wrap</code> = 0.1.4 ***(mandatory version)***

> <strong>Note:</strong> you can use a conda environment to setup this environment. In your terminal, type the following command (in your repository directory):
> ````
> conda env create -f doc/conda_env_dassflow-1d.yml
> ````
> To **activate** or **deactivate** the conda environment, you can use respectively <code>conda activate dassflow-1d</code> or <code>conda deactivate</code>.

<a name="ins"></a>
## 4. Install project
<ol>
  <li> Make sure all the <a href="#req">requirements</a> are met.</li>
  <li> Make sure the conda environment is activated using <code>conda activate dassflow-1d</code>.</li>
  <li> Open a terminal at the root of your repository.
  <li> Compile and install the required libraries: <br>
        
    make alllibs
  </li>
  <li> Install DassFlow1D <sup><a href="#fn1">[1]</a></sup>:
    
    make
  </li>
</ol>

<a id="fn1"><sup>[1]</sup> if you wish to use the adjoint code, see <a href="#adj">Compile and install the adjoint code</a> before you install using the <code>make</code> command.</a>

<a name="adj"></a>
## 5. Compile and install the adjoint code
### 5.1 Install Tapenade
To use the adjoint code for DassFlow1D, you need an **Automatic Differentiation (AD)** tool. Follow this <a href="https://tapenade.gitlabpages.inria.fr/tapenade/distrib/README.html">tutorial</a> to download and install <code>Tapenade</code> and add the following lines to your ````~/.bashrc```` to add tapenade to your PATH:
````
alias tapenade="tapenade_dir/bin/tapenade"
TAPENADE_HOME=tapenade_dir/bin
export PATH=$PATH:$TAPENADE_HOME
export PATH=$PATH:$"tapenade_dir"
````
<em>**Note:** tapenade_dir is the absolute path to the directory containing the tapenade files you just downloaded.</em>
### 5.2 Compile and install DassFlow1D with the adjoint code
<ol>
  <li> First, make sure all the <a href="#req">requirements</a> are met.</li>
  <li> Modify the <code>Makefile.inc</code> file to have <code>ADJOINT = 1</code> <em>(to tell the Makefile it needs to install the adjoint code for DassFlow1D)</em>.
  </li>
  <li> Then, generate the adjoint code by entering the following command (in your repository directory):

    make generate_adjoint
  </li>
  <li> Install DassFlow1D:
    
    make
  </li>
</ol>

<a name="next"></a>
## 6. First steps and Website
### 6.1 Test your DassFlow1D installation
To try and test your DassFlow1D installation, you can try to run a test case by entering the following commands in your terminal (in your repository directory):
````
conda activate dassflow-1d
. build/env.sh
cd cases/osse/strickler_inference/channel_powerlaw_spatial_fields
python run_osse_without_noise.py
````
> <strong>Note:</strong> If it works, last printed line should contain <code>CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH</code>.
### 6.2 Website
You can visit the <a href="" target="_blank">DassHydro</a> Website for more information on DassFlow1D and other versions.
