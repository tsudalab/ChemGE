#! -*- coding: utf-8 -*-

import csv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()

best_all = []
mean_all = []
worst_all = []
std_all = []
for i in range(1, 11):
    with open('log-J-1000-2000-{0:02d}'.format(i), 'r') as f:
        best = []
        mean = []
        worst = []
        std = []
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if row[0][0] == '$':
                continue
            best.append(float(row[1]))
            mean.append(float(row[2]))
            worst.append(float(row[3]))
            std.append(float(row[4]))
        best_all.append(best)
        mean_all.append(mean)
        worst_all.append(worst)
        std_all.append(std)

y = []
y1 = []
y2 = []
for i in range(500):
    best_generation = [best_all[t][i] for t in range(len(best_all))]
    m = np.mean(best_generation)
    s = np.std(best_generation)
    y.append(m)
    y1.append(m-s)
    y2.append(m+s)

plt.plot(np.arange(len(y)), y, label="best score")
plt.fill_between(np.arange(len(y)), y1, y2, alpha=0.5)

y = []
y1 = []
y2 = []
for i in range(500):
    mean_generation = [mean_all[t][i] for t in range(len(best_all))]
    m = np.mean(mean_generation)
    s = np.std(mean_generation)
    y.append(m)
    y1.append(m-s)
    y2.append(m+s)

plt.plot(np.arange(len(y)), y, label="mean score")
plt.fill_between(np.arange(len(y)), y1, y2, alpha=0.5)

plt.xlabel("generation")
plt.ylabel("J score")
plt.legend(loc="lower right")

plt.savefig("J-score-generation.pdf")
