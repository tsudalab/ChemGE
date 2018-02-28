#! -*- coding: utf-8 -*-

import csv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()

score = []
time = []
cnt = 0
with open('log-rdock-qsub-2', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0][0] == '%':
            time.append(float(row[2])/60.0)
            cnt = 0
        elif row[0][0] == 'l':
            break
        else:
            if cnt == 0:
                score.append(float(row[0]))
            cnt += 1

plt.plot(time, score, linestyle='-.', label="2 instances")

score = []
time = []
cnt = 0
with open('log-rdock-qsub-4', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0][0] == '%':
            time.append(float(row[2])/60.0)
            cnt = 0
        elif row[0][0] == 'l':
            break
        else:
            if cnt == 0:
                score.append(float(row[0]))
            cnt += 1

plt.plot(time, score, linestyle='--', label="4 instances")

score = []
time = []
cnt = 0
with open('log-rdock-qsub-8', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0][0] == '%':
            time.append(float(row[2])/60.0)
            cnt = 0
        elif row[0][0] == 'l':
            break
        else:
            if cnt == 0:
                score.append(float(row[0]))
            cnt += 1

plt.plot(time, score, linestyle='-', label="8 instances")

plt.xlabel("elapsed time (minute)")
plt.ylabel("rdock score")
plt.legend(loc="upper right")
plt.savefig("rdock-score-parallel.pdf")
