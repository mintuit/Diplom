import numpy as np
import msprime
import sys

i = sys.argv[1] # first argument

l = 30000000
mutation_rate = 1.25e-08
demography = msprime.Demography()
demography.add_population(name="A", initial_size=20000)
demography.add_population(name="B", initial_size=20000)
demography.set_symmetric_migration_rate(["A", "B"], 9.6e-5)


ts = msprime.sim_ancestry({'A': 1, 'B':2},sequence_length=l, demography=demography, recombination_rate=1.6e-9, ploidy=1,
                          random_seed=42+173)
mts = msprime.sim_mutations(ts,
                        rate=mutation_rate,
                        random_seed=147+173)
import copy
np.random.seed(173)
A = np.random.randint(4, size=(l,))
# print(A, len(A))
A = "".join(map(str, A))
# print(A, len(A))
A = A.replace('0', 'A').replace('1', 'T').replace('2', 'G').replace('3', 'C').replace('[', '').replace(']', '')
# print(A, len(A))
# A = A.split()
# print(A, len(A))
N = []
M = []
my_file = open(str(int(i)-1)+'_diplom_migr.fasta', 'w')
N.append(A)
for chr in mts.variants() :
    j = int(chr.position)
    if list(chr.genotypes)[int(i)-1] == 0 :
        N[0] = N[0][:j] + str(chr.alleles[0]) + N[0][j + 1 :]
    else :
        N[0] = N[0][:j] + str(chr.alleles[1]) + N[0][j + 1 :]
M.append("".join(map(str, N[0])))
my_file.write('>haplotype' + str(int(i)-1) + '\n')
my_file.write(M[0] + '\n')
my_file.close()
