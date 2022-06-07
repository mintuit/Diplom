Here I present the code of my work.

In the file "Simulation_seq.py" I simulate the situation with 2 populations with symmetric migration rate 9.7e-5. I get sequences, which are in the branch 'data'.

In the file 'From_Newick_to_MASCO.py' there is a code for 'translation' from Newick representation of tree into some representation that i need.

In the file 'Likelihood_and_migr_estimation.py' there is a code for calculating Likelihood of ARGs in the MASCO model using Importance Sampling; and there is a code
for finding Likelihood with different migration rates.

As a result in this situation we get such Heatmap:

<img width="282" alt="image" src="https://user-images.githubusercontent.com/82202781/172449998-122e2cdb-6bab-4904-9276-88ed3479a403.png">

Here we can see a square in the centre of the picture, that signals us that the Likelihood is maximum with migration rates about 0.8-1.1.
