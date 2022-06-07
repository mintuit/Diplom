import msprime
import numpy as np
import random
import argparse
import copy
from struct import pack, unpack
from sys import getsizeof
from scipy.integrate import RK45, solve_ivp
import matplotlib.pyplot as plt
from IPython.display import SVG
import copy
import sys
from scipy.optimize import minimize
import json
from scipy import special  # для биномиальных коэффициентов
import struct


class Tree :

    def __init__(self, migration_rates: np.array, number_of_samples: np.array, number_of_populations: int,
                 coalescence_rates: np.array, coalescence_time: float, Q: float, N: int) :
        # Извлечение данных из файла и создание объекта
        # Данные в файле должны быть в таком порядке N, k, n, t, m, q, Q
        # для отладки
        self.__N = int(N)  # эталонный размер популяции
        self.__number_of_populations = int(number_of_populations)  # количество популяций (m)
        self.__number_of_samples = np.array(number_of_samples)
        self.__samples_amount = int(np.sum(self.__number_of_samples))  # всего образцов дано (n)
        self.__T = float(coalescence_time)  # время слияния всех популяций
        self.__migration_probability = np.array(migration_rates)
        self.__coalescence_probability = np.array(coalescence_rates)
        self.__Q = float(Q)
        self.__cur_samples_amount = self.__samples_amount

    def show(self) :
        print(self.__number_of_populations, self.__T, '\n')
        print(self.__number_of_samples, self.__samples_amount, '\n')
        print(self.__migration_probability, '\n\n', self.__coalescence_probability, '\n')
        # s = self.get_initial_states
        # print(s)

    @property
    def original_size(self) :
        return self.__N

    @property
    def number_of_populations(self) :
        return self.__number_of_populations

    @property
    def number_of_samples(self) :
        return self.__number_of_samples

    @property
    def samples_amount(self) :
        return self.__samples_amount

    @property
    def T(self) :
        return self.__T

    @property
    def migration_probability(self) :
        return self.__migration_probability

    @property
    def coalescence_probability(self) :
        return self.__coalescence_probability

    @property
    def Q(self) :
        return self.__Q

    @property
    def tree_newick(self) :
        return self.__tree_newick

    @property
    def cur_samples_amount(self) :
        return self.__cur_samples_amount

    @cur_samples_amount.setter
    def cur_samples_amount(self, cur_samples_amount) :
        if cur_samples_amount >= 1 :
            self.__cur_samples_amount = cur_samples_amount
        else :
            raise ValueError("the minimum number of samples has been reached")

    def get_initial_states(self) :
        """
        return result 1-D array of initial states
        """
        result = []
        current_samples_sum = 0
        # по каждой популяции
        for i in range(self.__number_of_populations) :
            # по каждому образцу
            for j in range(self.__samples_amount) :
                # если образец в начальный момент принадлежит этой популяции
                if (current_samples_sum <= j < current_samples_sum + self.__number_of_samples[i]) :
                    result.append(1)
                else :
                    result.append(0)
            current_samples_sum += self.__number_of_samples[i]
        return result


def system_of_DE_for_lines(data: Tree, p: np.ndarray) -> np.ndarray :
    """
    :param data:
    :param p: p means P(L_i = l_i | T)
    :return: (result) derivative of probability function for each line in any population,
             1-D array with length = m * n
    """
    result = []
    value_of_cur_fun = 0
    # по каждой популяции
    for pop in range(data.number_of_populations) :
        # по всем образцам
        for i in range(data.cur_samples_amount) :
            # *** обращение к нужному p через индекс ( pop * data.cur_samples_amount + i )
            ################################################
            # первое(1-я часть семмы) слагаемое
            for a in range(data.number_of_populations) :
                value_of_cur_fun += data.migration_probability[a][pop] * p[a * data.cur_samples_amount + i]
            ################################################
            # первое(2-я часть суммы) слагаемое
            for a in range(data.number_of_populations) :
                value_of_cur_fun -= data.migration_probability[pop][a] * p[pop * data.cur_samples_amount + i]
            ################################################
            # второе слагаемое
            sum_of_mult = 0
            for a in range(data.number_of_populations) :
                for j in range(data.cur_samples_amount) :
                    if j != i :
                        for k in range(data.cur_samples_amount) :
                            if (k != i) and (k != j) :
                                t_sum1 = 0
                                for l in range(data.number_of_populations) :
                                    t_sum1 += p[l * data.cur_samples_amount + k]
                                t_sum2 = 0
                                for l in range(data.number_of_populations) :
                                    t_sum2 += p[l * data.cur_samples_amount + j]
                                sum_of_mult += (p[pop * data.cur_samples_amount + k] / t_sum1) * (
                                        p[pop * data.cur_samples_amount + j] / t_sum2)

            value_of_cur_fun -= sum_of_mult * 0.5 * data.coalescence_probability[a] * p[
                pop * data.cur_samples_amount + i]
            ################################################
            # третье слагаемое
            tmp_sum = 0
            for k in range(data.cur_samples_amount) :
                if k != i :
                    t_sum = 0
                    for l in range(data.number_of_populations) :
                        t_sum += p[l * data.cur_samples_amount + k]
                    tmp_sum += p[pop * data.cur_samples_amount + k] / t_sum
            value_of_cur_fun -= p[pop * data.cur_samples_amount + i] * data.coalescence_probability[pop] * tmp_sum
            ################################################
            result.append(value_of_cur_fun)
            value_of_cur_fun = 0

    return result


def create_initial0(data: Tree,
                    period: int,
                    previous_states: np.ndarray = None,
                    lineage: np.ndarray = None) -> np.ndarray :
    """
    :param data:
    :param period:
    :param previous_states: probabilities before coalescence
    :param lineage: contains two samples that participated in coalescence
    :return: result: vector of initial states for current "period"
    """
    # *** обращение к нужному p через индекс ( pop * data.cur_samples_amount + i )
    if int(period) == 0 :
        return data.get_initial_states()
    else :
        lineage = np.sort(lineage)
        result = []
        conditional_prob = []
        t_sum0 = 0
        for l in range(data.number_of_populations) :
            t_sum0 += previous_states[l * (data.cur_samples_amount + 1) + lineage[0]]
        t_sum1 = 0
        for l in range(data.number_of_populations) :
            t_sum1 += previous_states[l * (data.cur_samples_amount + 1) + lineage[1]]
        for pop in range(data.number_of_populations) :
            conditional_prob.append((previous_states[pop * (data.cur_samples_amount + 1) + lineage[0]]) *
                                    (previous_states[pop * (data.cur_samples_amount + 1) + lineage[1]] / t_sum1) *
                                    data.coalescence_probability[pop])
        sum_conditional_prob = np.sum(conditional_prob)  # ??? нужно ли переводить в массив ???
        # по каждой популяции
        for pop in range(data.number_of_populations) :
            # по каждому образцу
            for i in range(data.cur_samples_amount + 1) :
                # если образец участвовал в коалесценции и у него меньшее 'id', то есть он остался
                if i == lineage[0] :
                    # !!!!!!!!!!!!!! можно использовать уже посчитанные выше
                    result.append((previous_states[pop * (data.cur_samples_amount + 1) + lineage[0]]) *
                                  (previous_states[pop * (data.cur_samples_amount + 1) + lineage[1]] / t_sum1) *
                                  data.coalescence_probability[pop])
                # если образец участвовал в коалесценции и у него большее 'id', то есть он не остался
                elif i == lineage[1] :  # добавила для лучшего понимания кода, но вообще этот if можно убрать
                    continue
                # если образец не участвовал в коалесценции
                elif i != lineage[0] and i != lineage[1] :
                    t_sum2 = 0
                    for j in range(data.number_of_populations) :
                        t_sum2 += previous_states[j * (data.cur_samples_amount + 1) + i]
                    result.append(
                        (previous_states[pop * (data.cur_samples_amount + 1) + i] / t_sum2) * sum_conditional_prob)

            # cur_samples_sum += tree.number_of_samples[pop]  ??нужно??
        return result


def prop_1(limits_list: np.array, lineage_list: np.array, tree: Tree) :
    previous_states = None
    sol_lines_list = []
    sol_init_states = []
    for period in range(tree.samples_amount - 1) :
        print(period, 'period')
        lineage = lineage_list[period]
        print(limits_list, 'limlist')
        limits = limits_list[period]  # ...
        if period != 0 :
            previous_states = np.array(sol_lines_list[-1])[:, -1]
        initial_states = create_initial0(tree, period, previous_states=previous_states, lineage=lineage)
        t_span = np.linspace(limits[0], limits[1], 10001)
        q = np.array(solve_ivp(lambda t, y : system_of_DE_for_lines(tree, y),
                               t_span=limits,
                               y0=initial_states,
                               t_eval=t_span).y)
        q[q < 0] = np.exp(-15)
        sol_lines_list.append(q)
        tree.cur_samples_amount -= 1
    result = np.sum(
        create_initial0(tree, period=1, previous_states=np.array(sol_lines_list[-1])[:, -1], lineage=np.array([0, 1])))
    return result

def tree_eq(time: np.array, N: int, k: int):,
    lam = 1/N,
    sol = [],
    for i in range(k, 1, -1):",
        sol.append(lam  * np.exp(-lam * 0.5 * i * (i-1) * time[k - i]))
    return np.prod(sol)

def migr_est(migr_mat: np.ndarray, limits_list: np.ndarray, lineage_list: np.ndarray, num_replicates: int) :
    m12 = migr_mat[0]
    m21 = migr_mat[1]
    migration_matrix1 = [[0, m12],
                         [m21, 0]]
    M12 = m12 / 10000
    M21 = m21 / 10000
    migration_matrix2 = [[0, M12],
                         [M21, 0]]
    P1 = np.zeros(num_replicates)
    lineage_list1 = np.copy(lineage_list)
    limits_list1 = np.copy(limits_list)
    for i in range(num_replicates) :
        lean_time = np.array(limits_list1[i]).reshape(len(limits_list1[i]), 1)
        lean_time2 = np.insert(lean_time[:len(lean_time) - 1], 0, 0).reshape(len(lean_time), 1)
        lin_time = np.hstack((lean_time2, lean_time+0.00001))
        print(lin_time, lean_time)
        tree_1 = Tree(migration_rates=migration_matrix2, number_of_samples=[3, 3], number_of_populations=2,
                      coalescence_rates=[1, 1], coalescence_time=1000000, Q=1, N=10000)
        S = tree_eq(time = limits_list, N = 10000, k = 6)
        P1[i] = prop_1(limits_list=lin_time, lineage_list=lineage_list1[i], tree=tree_1) / (S * num_replicates)
    D1 = np.log10(np.sum(P1))
    return (- D1)


def count_lines(filename, chunk_size=1 << 13) :
    with open(filename) as file :
        return sum(chunk.count('\n')
                   for chunk in iter(lambda : file.read(chunk_size), ''))


lenenen = count_lines('masco.txt')

with open("masco.txt") as file :
    lines = file.readlines()

linlist = []
timelist = []

for i in range(lenenen) :
    line = ''.join(lines[i])
    lin = line.split(';')
    l1 = json.loads(lin[0])
    l2 = json.loads(lin[1])
    linlist.append(l1)
    timelist.append(l2)


x = np.linspace(0.01, 2, num=10)
Q = np.empty([10, 10])
for i in range(10):
    for j in range(10):
        x0 = np.array([x[i], x[j]], dtype=np.double)
        print(x0)
        Q[i, j] = migr_est(migr_mat = x0, 
                    limits_list=np.array(timelist[: :100]), 
                    lineage_list=np.array(linlist[: :100]),
                    num_replicates=lenenen // 100)
        
    
my_file = open('migr1.txt', 'w')
my_file.write(str(Q) + '\n')
my_file.close()
