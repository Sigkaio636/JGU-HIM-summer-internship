# JGU-HIM-summer-internship
This is the collection of the work carried out during the stay at Helmholtz-Institut Mainz at Johannes Gutenberg University in the summer of 2025. 

It consists of a literature review in the current situation of 7Li cluster theory, the codes created for fitting the models and the their analysis, the files from JGU Data Center, and the presentations from the meetings.

Highly recommended reading [Pending tasks and future extensions](##-pending-tasks-and-future-extensions) to find out at what point this work ends.

## 1. Literature and Datasets
### 📊📄 ./Datasets 📊📄
Contains the collected papers on experiments involving empirical measurements of 7Li cross sections. (see [Pending tasks and future extensions](##-pending-tasks-and-future-extensions))
Also, contains the papers of ANC and Phase-shift used data. 

> Cross section data extracted from https://www-nds.iaea.org/exfor can be found in [Racam_sca Jupyter notebook](Sfactor/Racam_sca.ipynb) in **Manual S factor extraction** section.

### 📑 ./Support papers 📑
It consists of papers related to the research topic that were found in the initial literature search. They were helpful in understanding the state of the issue at the beginning. (see [Pending tasks and future extensions](##-pending-tasks-and-future-extensions))

## 2. Analysis of results with Jupiter notebooks

### 📈 ./Graphic outfile 📈
These are simply the functions for reading Boscos.f output files in .dep and .dfo format.

### 📝🧮 ./Sfactor 📝🧮
These are Jupyter notebooks containing the data analysis of the Boscos.f and Racam5.f outputs, for the final model refits and the S-factor contribution calculations.
#### ▶ ./Sfactor/Racam_sca.ipynb
The initial calculations of S-factor contributions, first time using the Racam5.f code. The results are used for the [13/08/2025 meeting](<Meeting presentations/13-08 meet.pptx>).
The results are not are not definite, because later was performed a refit for central potential depth of p-waves. 
#### ▶ ./Sfactor/Complete-Sfactor.ipynb
This is the more extensive and last analysis. It starts again from the beginning by evaluating the correspondence of Tokimoto's potentials, applying what it was learned earlier when attempting to adjust the binding energy (the convention of mass unit and Coulomb radius choice and central vs spin-orbit models don't have a significant advantage among them). It focuses solely on fitting with empirical radius and masses. It is explained how the process of building the latest models went.

Furthermore, as both Fortran codes are used, it has been necessary to ensure that both give the same results. (see [Complete HPC directory](##-complete-hpc-directory))

The results are the final ones presented in the last presentation [25/08/2025 meeting](<Meeting presentations/Corrected 25-08 final meet.pptx>). 

## 3. Refitting models and optimization problems

### ⚙️ ./Refit dev ⚙️
These are the codes for visualising the relationship between potential depth and binding energy for bound states (p3/2 and p1/2), and solving the optimisation problem of finding the potential depth *Vc* to reproduce this observable. **These codes are NOT USED to construct the final potential models.**

/LossField1D, /LossField2D and /VcVso_opt are folders with executable scripts. /Phaseshifts_Ind contains just data analysis.
#### ▶ ./Refit dev/LossField1D and ./Refit dev/LossField2D
**LossField1D**: Generates the tuples (*Vc*, *Eb*) of the corresponding binding energy for each central potential depth Vc given a state (p3/2 or p1/2) and convention. These last parameters are fixed in the input file. Calculations are performed in [gonzalee/refitInd/LossVisul](gonzalee/refitInd/LossVisul) and [gonzalee/refitInd/mesh_res](gonzalee/refitInd/mesh_res).

**LossField2D**: Generates the tuples (*Vc*, *Vso*, *Eb*) of the corresponding binding energy for each pair of central potential depth Vc and spin-orbit strenght Vso given a state (p3/2 or p1/2) and convention. These last parameters are fixed in the input file. Calculations are performed in [gonzalee/refitInd/LossVisulSurf](gonzalee/refitInd/LossVisulSurf)

These two folder have a very similar execution line:

  **+ config.yaml** : It is where the HPC access is stated (⚠️ This should be changed to your own configuration) and the parameters for sampling from the linear fine mesh.

  **+ create_inpfile.py** : Each input file (e.g. 042_input.dat) for Bocos.f is created. State and conventions is specified here by changing the template.

  **+ create_outcsv.py** : Reads all previously generated outfiles (e.g. 033_output.txt), extracts the binding energy for each of them and collects them into a final CSV file (e.g. res_run-20250803T192222.csv).

  **+ local_runner.py** : This is the main code, where execution is initialised. **You must execute only this script in your own machine**. It is recommended to read the comments to understand the inner workings.

  **+ graph_results.py** : The previously generated CSV file is read and plotted:

* /LossField1D/graph_results.py : It allows to find the optimal *Vc* simply by zooming in on the figure, provided a sufficient resolution in the fine mesh.
* /LossField2D/graph_results.py : It uses the pre-calculated CSV files in the folder to visualise the surfaces (*Vc*, *Vso*, *EbGround*) (*Vc*, *Vso*, *EbExcited*) and construct the loss function for the gradient descent of ./VcVso_opt

#### ▶ ./Refit dev/VcVso_opt
It contains the executables which perform the refit of central potential *Vc* and spin-orbital strength *Vso* to reproduce simultaneously binding energy for both ground and excited states. The optimization method consists in a gradient descent with an adecuate choice of constant learning rate. Visually explained in [13/08/2025 meeting](<Meeting presentations/13-08 meet.pptx>). 

  **+ config.yaml** : It is where the HPC access is stated (⚠️ This should be changed to your own configuration) and the initial point for gradient descent.

  **+ exploration csv/** : Contains several execution logs for different values of learning rate.

  **+ graph_results.ipynb** : Analysis studying the properties of gradient descent. Used to find the optimal learning rate. It uses the data from /exploration csv/.

  **+ local_runner.py** : This is the main code, where execution is initialised. **You must execute only this script in your own machine.**

  **+ remote_runner.py** : File that is copied to the remote machine and manages execution on that machine. It contains the functions for input creation, execution with Boscos.f, output reading, csv writing, and the entire implementation of the gradient descent.

  **+ resultant dep_dfo/** : Contains the output dep and dfo files from Boscos, for each convention with the optimum potential parameters.

  **+ PhaseRefit_Dep.ipynb** : Plots the phase-shifts and wavefunction for refited potential models. It uses the data from /resultant dep_dfo/.

(see [Pending tasks and future extensions](##-pending-tasks-and-future-extensions))

## 4. Complete HPC directory
### 💻🗃️ ./gonzalee 💻🗃️
This folder contains the content that was in my *kphth* user account. **It is especially interesting for the .dat input files 🔼**.
#### ▶ ./gonzalee/manual
The scripts and input files from the Boscos.f and Racam5.f manual examples. Used initially to familiarize oneself with the Fortran codes.
#### ▶ ./gonzalee/BindE
Exploration of Boscos.f script changing the radius and mass unit to reproduce the binding energy for both states. Results presented in [13/08/2025 meeting](<Meeting presentations/13-08 meet.pptx>). Here, the conventions were defined as a pair of Coulomb radius and mass unit. 
#### ▶ ./gonzalee/refitInd
Contains the directories where the calculations to visualize the relation between model parameters (*Vc*, *Vso*) and binding energy (*EbGround*, *EbExcited*) are performed.

  **/LossVisul** : Calculates the tuples (*Vc*, *Eb*) for fitting models with just central potential. 

  **/LossVisulSurf** : Calculates the tuples (*Vc*, *Vso*, *Eb*) to construct the loss function of the gradient descent.

  **/mesh_res** : Contains the input and output files for the fitted central potential models for each convention.

The name of the file indicates the state (*G* for ground = p3/2 and *E* for excited = p1/2) and later references the convention. 

  🔴 *eki* : My suggestion for convention, since it reproduced the binding energy of 8B in [30/07/2025 meeting](<Meeting presentations/30-07 meet.pptx>). **Estimated** Coloumb radius 2.33951 fm and **empirical** masses.

  🟡 *ekX* : An extension of my suggested convention. Instead of using the estimated radius, it uses the empirically measured radius of 7Li (from Dubovichenko 2009). **Empirical** Coloumb radius 2.35 fm and **empirical** masses.  

  🔵 *pieRM* : Convention suggested by Pierre. It uses integer multiples of the average between proton and neutron masses and for the Coloumb radius the same as for Woods-Saxon geometry in Tokimoto's paper. **Woods-Saxon** Coloumb radius 2.39 fm and **integer** masses.

  🟢 *pieM* : Modification of Pierre's convention. It uses the empirical radius of 7Li instead of Woods-Saxon's. **Empirical** Coloumb radius 2.35 fm and **integer** masses.

In [13/08/2025 meeting](<Meeting presentations/13-08 meet.pptx>) there is a summary table. 

#### ▶ ./gonzalee/OptimVcVso
This is the directory where the execution of gradient descent is performed. 

**/phsh_wf** : Contains the input and output files for the fitted central and spin-orbit potential models for each convention. Same file names as in ./gonzalee/refitInd/mesh_res are used.

#### ▶ ./gonzalee/radiContri
Various tests to check for consistency between Boscos and Racam code outputs. The Fortran scripts had to be modified to change the value of the constants and then recompiled.
+ Boscos : Original Boscos.f code
+ BoscosE : Boscos.f modified to show more decimals in binding energy and ANC.
+ Racam5 : Original Racam5.f code
+ Racam5eps : Racam5.f modified to have the same energy diferential as in Boscos.f of eps=1e-6, originally eps=1e-3 .

Important to use empirical masses as hm = hbar^2/(2*mu) in Racam is calculated in function of other constants.

#### ▶ ./gonzalee/CheckTokphsh
All input and output files used in [Complete-Sfactor.ipynb](Sfactor/Complete-Sfactor.ipynb). As the convention is unique here, empirical Coloumb radius and empirical masses, different file names correspond to different partial waves (e.g. Tok_dp.dat = Tokimoto's potential for d5/2 wave) (e.g. refit_dm.dat = Refitted model for d3/2 wave)
If file name starts with *WS*, it is for Racam5eps.out. Any questions regarding the names can be resolved by consulting the *Vc*, *Vso* values in the input file.

## 5. Pending tasks and future extensions

+ I couldn't get Griffiths's paper. Johannes Gutenberg University do not have permission to access it, but other german universities do have them. See https://cdnsciencepub.com/doi/10.1139/p61-167

+ From [Support papers](<Support papers>) they are expecially interesting the following papers. First three correspond to the same team and explain that the breakup of 6Li and 7Li is significantly mediated by a transfer reaction which create different channels. Last one is the most modern paper founded about the topic, also considers differents channels. It might be interesting to compare their results with the hypotheses used in the data analysis in the various empirical measurements articles.
  + [Insights into mechanisms and time-scales of breakup of 6,7Li; 2011](<Support papers/Insights into mechanisms.pdf>)
  + [Asymptotic and near-target direct breakup of 6Li and 7Li; 2016](<Support papers/Asymptotic near-target.pdf>)
  + [Challenges in describing nuclear reaction outcomes at near-barrier energies; 2017](<Support papers/Challenges describing outcomes.pdf>)
  + [Triton-alpha radiactive capture reaction at astrophysical energies; 2023](<Support papers/Faddeev eq.pdf>)

+ Sfactor errors could not be extracted for Griffiths' data. See [Racam_sca.ipynb](Sfactor/Racam_sca.ipynb)
  > Nor the errors of Ivanovich phase shift measurements. But as in Dubovichenko paper Fig 1 don't appear errorbars, we assume that they can't be recovered.

+ In [Exfor](https://www-nds.iaea.org/exfor), they aren't the Sfactor data of Burzynski paper, just the angular differential cross section.

+ If you wish to use the gradient descend code, it would be beneficial to make some simple optimizations. In particular, instead of generating all the input files for each iteration, to save space, overwrite a file. The generation of all was useful for debugging and ensuring that the code worked properly, but now it is no longer necessary.
  > A better learning rate could also be sought, perhaps by making it dynamic.
