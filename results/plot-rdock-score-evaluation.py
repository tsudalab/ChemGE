#! -*- coding: utf-8 -*-

import csv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()

evaluation = []
score = []
with open('log-rdock-score-evaluation.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        evaluation.append(int(row[0]))
        score.append(float(row[1]))

plt.plot(evaluation, score, label="ChemGE (32, 64)")

scores = []
with open('log-zinc-40000', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        scores.append(float(row[0]))

y = []
y1 = []
y2 = []
mins = scores[0::10000]
for i in range(1, 10000):
    mins2 = scores[i::10000]
    mins = [min(mins[j], mins2[j]) for j in range(len(mins))]
    m = np.mean(mins)
    s = np.std(mins)
    y.append(m)
    y1.append(m-s)
    y2.append(m+s)

plt.plot(np.arange(len(y)), y, label="ZINC")
plt.fill_between(np.arange(len(y)), y1, y2, alpha=0.5)

plt.xlabel("evaluation")
plt.ylabel("rdock score")
plt.legend(loc="upper right")

plt.savefig("rdock-score-evaluation.pdf")
