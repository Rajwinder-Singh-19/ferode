import csv
import matplotlib.pyplot as plt

t_undamped = []
y_undamped = []
with open("undamped_result.csv") as file:
    reader = csv.reader(file)
    for row in reader:
        t_undamped.append(float(row[0]))
        y_undamped.append(float(row[1]))

t_damped = []
y_damped = []
with open("damped_result.csv") as file:
    reader = csv.reader(file)
    for row in reader:
        t_damped.append(float(row[0]))
        y_damped.append(float(row[1]))

fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].plot(t_undamped, y_undamped, linewidth=4)
ax[0].set_title("y'' + y = 0 [Undamped Harmonic Oscillations]", fontweight="bold", fontsize=16)
ax[0].set_xlabel("t", fontweight="bold", fontsize=16)
ax[0].set_ylabel("y", fontweight="bold", fontsize=16)
ax[1].set_xlabel("t", fontweight="bold", fontsize=16)
ax[1].set_ylabel("y", fontweight="bold", fontsize=16)
ax[0].tick_params(axis='x', labelsize=12)
ax[0].tick_params(axis='y', labelsize=12)
ax[1].tick_params(axis='x', labelsize=12)
ax[1].tick_params(axis='y', labelsize=12)
ax[1].plot(t_damped, y_damped, linewidth=4)
ax[1].set_title("y'' + cy' + y = 0 [Damped Harmonic Oscillations]", fontweight="bold", fontsize=16)
ax[0].grid(True)
ax[1].grid(True)
plt.show()
