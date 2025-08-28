# JGU-HIM-summer-internship
This is the collection of the work carried out during the stay at Helmholtz-Institut Mainz at Johannes Gutenberg University in the summer of 2025. 

It consists of the literature found, the performed data analysis, the codes created for fitting the models, the files from JGU Data Center, and the presentations from the meetings.

## 1. Literature and Datasets
### 📊📄 ./Datasets 📊📄
Contains the collected papers on experiments involving empirical measurements of 7Li cross sections. (see [Pending tasks and future extensions](##-pending-tasks-and-future-extensions))
Also, contains the papers of ANC and Phase-shift used data. 

> Cross section data extracted from https://www-nds.iaea.org/exfor can be found in [Racam_sca Jupyter notebook](Sfactor/Racam_sca.ipynb) in **Manual S factor extraction** section.

### 📑 ./Support papers 📑
It consists of papers related to the research topic that were found in the initial literature search. (see [Pending tasks and future extensions](##-pending-tasks-and-future-extensions))

## 2. Analysis of results with Jupiter notebooks

### 📈 ./Graphic outfile 📈
These are simply the functions for reading Boscos.f output files in .dep and .dfo format.

### 📝🧮 ./Sfactor 📝🧮
These are Jupyter notebooks containing the data analysis of the Boscos.f and Racam5.f outputs, for the final refits of the models and the S-factor contribution calculations.
#### ▶ ./Sfactor/Racam_sca.ipynb
The initial calculations of S-factor, first time using the Racam5.f code. The results are used for the [13/08/2025 meeting](<Meeting presentations/13-08 meet.pptx>).
The results are not are not definite, because later was performed a refit for central potential depth of p-waves. 
#### ▶ ./Sfactor/Complete-Sfactor.ipynb
This is the more extensive and last analysis. It starts again from the beginning by evaluating the correspondence of Tokimoto's potentials, applying what it was learned earlier when attempting to adjust the binding energy. It focuses solely on fitting with empirical radius and masses. He explains how the process of building the latest models went.

Furthermore, as both Fortran codes are used, it has been necessary to ensure that both give the same results. (see [Complete HPC directory](##-complete-hpc-directory))

The results are the final ones presented in the last presentation [25/08/2025 meeting](<Meeting presentations/Corrected 25-08 final meet.pptx>). 

## 3. Refitting models and optimization problems

### ⚙️ ./Refit dev ⚙️
These are the codes for visualising the relationship between potential depth and binding energy for bound states (p3/2 and p1/2), and solving the optimisation problems of finding the optimal depth to reproduce this observable. **These codes are NOT USED to construct the final potential models.**

/LossField1D, /LossField2D and /VcVso_opt are folders with executable scripts. /Phaseshifts_Ind contains just data analysis.
#### ▶ ./Refit dev/LossField1D and ./Refit dev/LossField2D
**LossField1D**: Generates the tuples (Vc, Eb) of the corresponding binding energy for each central potential depth Vc given a state (p3/2 or p1/2) and convention. These last parameters are fixed in the input file. Calculations are performed in (TODO Check refit_ind LossVisul mesh_res)

**LossField2D**: Generates the tuples (Vc, Vso, Eb) of the corresponding binding energy for each pair of central potential depth Vc and spin-orbit strenght Vso given a state (p3/2 or p1/2) and convention. These last parameters are fixed in the input file. Calculations are performed in (TODO Check refit_ind LossVisulSurf )

These two folder have a very similar execution line:
##### + config.yaml : It is where the HPC access is stated (⚠️ This should be changed to your own configuration) and the parameters for sampling from the linear fine mesh.
##### create_inpfile.py : Each input file (e.g. 042_input.dat) for Bocos.f is created. State and conventions is specified here by changing the template.
##### create_outcsv.py : Read all previously generated outfiles (e.g. 033_output.txt), extract the binding energy for each of them and collects them into a final CSV file (e.g. res_run-20250803T192222.csv).
##### local_runner.py : This is the main code, where execution is initialised. You must execute only this script in your own machine. It is recommended to read the comments to understand the inner workings.
##### graph_results.py : Here, the previously generated CSV file is read and plotted. 

* /LossField1D/graph_results.py : It allows to find the optimal Vc simply by zooming in on the figure, provided a sufficient resolution in the fine mesh.
* /LossField2D/graph_results.py : It uses the pre-calculated CSV files in the folder to visualise the surfaces (Vc, Vso, EbGround) (Vc, Vso, EbExcited) and construct the loss function for the gradient descent of ./VcVso_opt

#### ▶ ./Refit dev/VcVso_opt
It contains the executables which perform the refit of central potential Vc and spin-orbital strength Vso to reproduce simultaneously binding energy for both ground and excited states. The optimization method consists in a gradient descent with an adecuate choice of constante learning rate.

##### config.yaml : It is where the HPC access is stated (⚠️ This should be changed to your own configuration) and the initial point for gradient descent.
##### exploration csv/ : Contains several execution logs for different values of learning rate
##### graph_results.ipynb : Analysis studying the properties of gradient descent. Used to find the optimal learning rate. It uses the data from /exploration csv/.
##### local_runner.py : This is the main code, where execution is initialised. You must execute only this script in your own machine.
##### remote_runner.py : File that is copied to the remote machine and manages execution on that machine. It contains the functions for input creation, execution with Boscos.f, output reading, csv writing, and the entire implementation of the gradient descent.
##### resultant dep_dfo/ : Contains the output dep and dfo files from Boscos, for each convention with the optimum potential parameters.
##### PhaseRefit_Dep.ipynb : Plots the phase-shifts and wavefunction for refited potential models. It uses the data from /resultant dep_dfo/.


## 4. Complete HPC directory
### 💻🗃️ ./gonzalee 💻🗃️
This folder contains the content that was in my *kphth* user account. ***It is especially interesting for the .dat input files***.
#### ▶ ./gonzalee/manual
The scripts and input files from the Boscos.f and Racam5.f manual examples. Used initially to familiarize oneself with the Fortran codes.
#### ▶ ./gonzalee/BindE
Exploration of Boscos.f script changing the radius and mass unit to reproduce the binding energy for both states. Results presented in [13/08/2025 meeting](<Meeting presentations/13-08 meet.pptx>). Here, the conventions were defined as a pair of coulomb radius and mass unit. 
#### ▶ ./gonzalee/refitInd


## Pending tasks and future extensions

> Falta pdf de Griffiths

> Particularmente interesantes los papers de "apuntado"

> No se pudieron extraer los errores para Sfactor data Griffiths

> Exfor no estaba los datos de Sfactor para el paper de Burzynski, solo  angular differential cross section.

> Optimización de gradinet descend en caso de utilizarse
