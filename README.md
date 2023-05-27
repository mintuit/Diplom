 # Coalescent Model and its Approximations for Population Genetics.

## Introduction.
Coalescent theory is a probabilistic model for genealogical trees. It models ge- nealogies as a Markovian process with backwards time (ancestral lineages go from descendants to their ancestors). Structured coalescence is a generalisation of a ba- sic model to multiple populations with migration. Though efficient for simulations, the inference is computationally challenging under structured coalescent. MASCO is an approximation of coalescent which allows to calculate tree likelihoods efficiently for complicated population scenarios. In our work we develop an approach for com- puting Ancestral Recombination Graphs’ likelihood using Importance Sampling and MASCO.

 The main part of this work is to analyze how we can use MASCO in estimation of some parameters, like migration rates, population size using ARGs, additionally using Importance sampling. The first part of work is to simulate such situation (and the sequences of length 3 · 107) that we have 2 populations which population sizes are 10000. We suppose the symmetric migration between populations with the migration rate 9.6e-1. We observe 4 samples from each population. To do this we use a python software msprime. 

 The second part is to simulate a collection of ARGs of these 6 sequences. In this step we use the ARGweaver [4]. We get 1000 ARGs and get from them every 50th graph. In the ARG we do not take every tree, because of the correlation between neighbours trees, so choose every 100th tree. 

 The last part is recalculating the ARGs’ Likelihood using Structured Coalescent theory and Importance sampling. To calculate ARG’s Likelihood in MASCO model we use a code that solve the MASCO equations numerically using Runge-Kutta method. Then we use Importance Sampling to calculate the Likelihood of ARGs, where 𝑄(𝐺𝑖) is the likelihood of tree in Kingman Coalescent model, 𝑃Θ(𝐺𝑖) is the likelihood in MASCO model and 𝑃 𝑟Θ(𝐷|𝐺𝑖) is the probability of observed sequences.

For a detailed explanation you can read the file "Coalescent_Model_and_its_Approximations_for_Population_Genetics.pdf"

## Code navigation.

In the file "Simulation_seq.py" I simulate the situation with 2 populations with symmetric migration rate 9.7e-5. I get sequences, which are in the branch 'data'.

In the file 'From_Newick_to_MASCO.py' there is a code for 'translation' from Newick representation of tree into some representation that i need.

In the file 'Likelihood_and_migr_estimation.py' there is a code for calculating Likelihood of ARGs in the MASCO model using Importance Sampling; and there is a code
for finding Likelihood with different migration rates.

## Results.

This heatmap represant the quality of migration rates estimation. The axes show the migration rates from one population to another. 

<img width="282" alt="image" src="https://user-images.githubusercontent.com/82202781/172449998-122e2cdb-6bab-4904-9276-88ed3479a403.png">

We can see a square in the centre of the picture, that signals us that the Likelihood is maximum with symmetric migration rates about 0.8-1.1.

