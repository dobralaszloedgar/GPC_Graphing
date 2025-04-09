import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np

# Initial variables
txt_file = "G:/Edgar Dobra/GPC Samples/Spring 2025/03.04.2025_GB-OT2_PS-b.txt"
color_palette = 2
RI_calibration = "G:/Edgar Dobra/GPC Samples/Calibration Curves/RI Calibration Curve 2025 January.txt"
x_lim = [1e3, 1e8]  # default set to [1e3, 1e8]
y_lim = [-0.05, 1.2]  # default set to [-0.05, 1.2]
normalize_peak_range = [1e4, 1e8]  # default set to [1e4, 1e8]


matplotlib.rcParams['font.family'] = 'Avenir Next LT Pro'
matplotlib.rcParams['font.size'] = 18
plt.figure(figsize=(7, 5))
plt.subplots_adjust(bottom=0.19)
plt.subplots_adjust(left=0.19)

file = open(txt_file, "r")
content = file.readlines()
file.close()
list_content = []
for row in content:
    list_content.append(row.strip("\n").split("\t"))
data_array = np.array(list_content)

file_RI = open(RI_calibration, "r")
content_RI = file_RI.readlines()
file_RI.close()
list_content_RI = []
for row in content_RI:
    list_content_RI.append(row.strip("\n").split("\t"))
data_array_RI = np.array(list_content_RI)


def color(number):
    if color_palette == 1:
        if number == 0:
            return "#FA7921"
        elif number == 1:
            return "#5BC0EB"
        elif number == 2:
            return "#E55934"
        elif number == 3:
            return "#9BC53D"
        elif number == 4:
            return "#4059AD"
        elif number == 5:
            return "#037171"
        elif number == 6:
            return "#FF3562"

    elif color_palette == 2:
        if number == 0:
            return "#A4036F"
        elif number == 1:
            return "#048BA8"
        elif number == 2:
            return "#16DB93"
        elif number == 3:
            return "#F29E4C"
        elif number == 4:
            return "#EFEA5A"
        elif number == 5:
            return "#A18276"

    elif color_palette == 3:
        if number == 0:
            return "#FF5F0F"
        elif number == 1:
            return "#0455A4"

    elif color_palette == 4:
        if number == 0:
            return "#DF2E38"
        elif number == 1:
            return "#5D9C59"

    elif color_palette == 5:
        if number == 0:
            return "#DF2E38"
        elif number == 1:
            return "#5D9C59"

    elif color_palette == 11:
        if number == 0:
            return "#FA7921"
        elif number == 1:
            return "#5BC0EB"
        elif number == 2:
            return "#037171"
    elif color_palette == 12:
        if number == 0:
            return "#E55934"
        elif number == 1:
            return "#9BC53D"
        elif number == 2:
            return "#037171"
    elif color_palette == 13:
        if number == 0:
            return "#4059AD"
        elif number == 1:
            return "#037171"
    elif color_palette == 14:
        if number == 0:
            return "#A4036F"
        elif number == 1:
            return "#048BA8"
    elif color_palette == 15:
        if number == 0:
            return "#16DB93"
        elif number == 1:
            return "#F29E4C"
    elif color_palette == 16:
        if number == 0:
            return "#EFEA5A"
        elif number == 1:
            return "#A18276"
    elif color_palette == 17:
        if number == 0:
            return "#C5D86D"
        elif number == 1:
            return "#0D1321"
    elif color_palette == 18:
        if number == 0:
            return "#42113C"
        elif number == 1:
            return "#6BD425"
    elif color_palette == 19:
        if number == 0:
            return "#663F46"
        elif number == 1:
            return "#3C362A"
    elif color_palette == 20:
        if number == 0:
            return "#6874E8"
        elif number == 1:
            return "#F7ACCF"

    elif color_palette == 101:
        if number == 0:
            return "#ea5545"
        elif number == 1:
            return "#f46a9b"
        elif number == 2:
            return "#ef9b20"
        elif number == 3:
            return "#edbf33"
        elif number == 4:
            return "#ede15b"
        elif number == 5:
            return "#bdcf32"
        elif number == 6:
            return "#87bc45"
        elif number == 7:
            return "#27aeef"
        elif number == 8:
            return "#b33dc6"
        elif number == 9:
            return "#964B00"

    elif color_palette == 102:
        if number == 0:
            return "#e60049"
        elif number == 1:
            return "#0bb4ff"
        elif number == 2:
            return "#50e991"
        elif number == 3:
            return "#e6d800"
        elif number == 4:
            return "#9b19f5"
        elif number == 5:
            return "#ffa300"
        elif number == 6:
            return "#dc0ab4"
        elif number == 7:
            return "#b3d4ff"
        elif number == 8:
            return "#00bfa0"
        elif number == 9:
            return "#80461B"


def max_of_y_within_range(x_array, y_array):
    mask = np.logical_and(x_array > normalize_peak_range[0], x_array < normalize_peak_range[1])
    return np.max(y_array[mask])


def min_of_y_within_range(x_array, y_array):
    mask = np.logical_and(x_array > 1e3, x_array < 1e4)
    return np.min(y_array[mask])


def get_data(index):
    x_a = np.delete(data_array_RI[:, 1], [0, 1], 0)
    x_a = x_a.astype(float)
    x_a = 10**x_a
    y_a = np.delete(data_array[:, index * 2 + 1], [0, 1], 0).astype(float)

    min_y = min_of_y_within_range(x_a, y_a)
    for value in range(len(y_a)):
        y_a[value] = y_a[value] - min_y

    max_y = max_of_y_within_range(x_a, y_a)
    for value in range(len(y_a)):
        y_a[value] = y_a[value] / max_y

    return x_a, y_a


for i in range(int(len(data_array[0, :])/2)):
    x = get_data(i)[0]
    y = get_data(i)[1]
    plt.plot(x, y, color(i), label=data_array[0, 2*i])


plt.xscale("log")
ax = plt.gca()
ax.set_xlim(x_lim)
ax.set_ylim(y_lim)
plt.xlabel("Molecular weight (g/mol)", fontsize=20, name='Avenir Next LT Pro', fontstyle='italic', fontweight='demi')
plt.ylabel("Normalization", fontsize=20, name='Avenir Next LT Pro', fontstyle='italic', fontweight='demi')
plt.legend()
plt.savefig("Figure1.png", dpi=1200)
plt.show()
