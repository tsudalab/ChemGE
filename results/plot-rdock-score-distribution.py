#! -*- coding: utf-8 -*-

import csv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()

scores = []
with open('log-rdock-zinc', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        score = float(row[0])
        if score < 100:
            scores.append(score)
sns.distplot(scores, label="ZINC")
plt.xlim(-50, 100)

scores = []
with open('mol-rdock-known-active.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        score = float(row[0])
        if score < 100:
            scores.append(score)
sns.distplot(scores, label="known active")

scores = []
with open('mol-rdock-chemge-all.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        score = float(row[0])
        if score < 100:
            scores.append(score)
sns.distplot(scores, label="ChemGE")

plt.xlim(-60, 100)
plt.legend()
plt.xlabel("docking score")
plt.ylabel("ratio")
plt.savefig("rdock-score-dist.pdf")
