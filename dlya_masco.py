import numpy as np
import copy
import re

with open("out.30.smc", "r") as file:
    lines = file.readlines()
lines = lines[2:]
lines = lines[::2]

def newick_to_masco1(tree, n):
    lean_time = []
    lean_nd = []
    a = str(tree)
    y = []
    y.append(a)
    for o in range(n-1):
        b = y[o].count('(')
        s1 = re.sub("[0-9]", "", y[o]).replace('.', '').replace(':', '').replace(',', '')
        s2 = s1.split('()')
        l = len(s2)
        t = []
        for i in range(l - 1):
            t.append(len(s2[i]) + 2)
        t.append(0)
        for i in range(l-1):
            t[i] = t[i-1] + t[i]
        del t[-1]
        q = y[o].replace('(', ' ').replace(")", " ")
        q = q.split(' ')
        m=[]
        b = []
        for j in range(len(t)):
            k = int(t[j]) - 1
            b.append(q[k])
        for i in range(len(q)):
            if len(q[i]) > 0 and q[i][0] == ':' and q[i].count(',') > 0:
                abc = q[i].split(',')
                q[i] = abc[0]
                q.insert(i + 1, abc[1])
                i = i + 1
        for i in range(len(q)):
            if q[i].count(':') == 2:
                l = q[i].replace('|', '')
                m.append(l)
        s = []
        sec = []
        sc = []
        for i in range(len(b)):
            sec.append(b[i].split(':'))
            s.append(sec[i][-1])
            sc.append(float(sec[i][-1]))
        mi1 = min(sc)
        j1 = sc.index(mi1)
        mi2 = sec[j1][1].split(',')[0]
        mi1 = sec[j1][-1]
        x1, x2 = sec[j1][0], sec[j1][1].split(',')[-1]
        red2 = str(x1)+':'+str(mi2)+','+str(x2)+':'+str(mi1)
        lean_time.append(float(mi1))
        lean_nd.append([int(x1), int(x2)])
        j2 = q.index(red2)
        if o < (n - 2):
            d1 = q[j2 + 1].replace(':', '?').replace(',', '?').split('?')
            d = d1[1].replace(':', '').replace(',', '').replace('|', '')
            red1 = '('+str(x1)+':'+str(mi2)+','+str(x2)+':'+str(mi1)+')'+str(d1[0])+':'+str(d)
            su = float(mi1) + float(d)
            z = str(n + o)+':' + str(su)
            y.append(y[o].replace(red1, z))
        else:
            red1 = 0
            su = str(mi1)
    return lean_time, lean_nd

def newick_to_msprime(lin_nodes, n):
    res = []
    res1 = []
    for i in range(len(lin_nodes)):
        res.append([lin_nodes[i][0], int(i + n)])
        res.append([lin_nodes[i][1], int(i + n)])
    res.sort()
    for i in range(len(res)):
        if i == 0:
            res1.append('{' + str(res[i][0]) + ': ' + str(res[i][1]) + ',')
        elif i == len(res) - 1:
            res1.append(' ' + str(res[i][0]) + ': ' + str(res[i][1]) +  '}')
        else:
            res1.append(' ' + str(res[i][0]) + ': ' + str(res[i][1])+  ',')
    res1 = ''.join(res1)
    return res1

def msprime_to_masco(tree):
    t = []
    a = str(tree)
    b = a.replace(":", "").replace("{", "").replace("}", "").replace(" ", "").split(",")
    d = []
    for i in range(len(b)):
        d.append(int(b[i][len(str(i)):]))
    a = np.array(d)
    n = int(len(a)/2 + 1)
    m = copy.deepcopy(n)
    z = []
    for i in range(n-1):
        p, q = np.where(a == m)[0][0], np.where(a == m)[0][1]
        a = a - 2
        b = np.zeros(len(a) - 2)
        for j in range(len(a) - 2):
            if j < q and j !=p:
                b[j] = a[j]
            elif j == p:
                b[j] = a[m]
            elif q <= j < m - 1:
                b[j] = a[j + 1]
            elif m - 1 <= j:
                b[j] = a[j+2]
        m -= 1
        a = np.copy(b)
        z.append([p, q])
    return z

my_file = open('masco.txt', 'w')
lenen = len(lines)
for i in range(lenen):
    s = lines[i]
    k = s.find("(")
    lines[i] = lines[i][k:]
    s1 = lines[i].replace('[', '|').replace(']', '|')
    s1 = s1.split('|')
    s2 = s1[0:len(s1) - 3]
    s2 = ''.join(s2[: : 2]) + ')'
    lean_time, lean_nd = newick_to_masco1(s2, 6)
    lin_node = msprime_to_masco(newick_to_msprime(lean_nd, 6))
    my_file.write(str(lin_node) + ';' + str(lean_time) + '\n')
my_file.close()