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
ax[0].plot(t_undamped, y_undamped)
ax[0].set_title("y'' + y = 0 [Undamped Harmonic Oscillations]")
ax[1].plot(t_damped, y_damped)
ax[1].set_title("y'' + cy' + y = 0 [Damped Harmonic Oscillations]")
ax[0].grid(True)
ax[1].grid(True)
plt.show()
