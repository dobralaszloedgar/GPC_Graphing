import os

import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np

# Initial variables
txt_file = f"G:/Edgar Dobra/GPC Samples/2025 Spring/JBB Campaign 1/6_GT-GF_RI.txt"
output_dir = r"G:/Edgar Dobra/Polymer GPC Graphs/2025 Spring/JBB Campaign 1"
os.makedirs(output_dir, exist_ok=True)  # create folder if it doesn’t exist
color_palette = 98
RI_calibration = "G:/Edgar Dobra/GPC Samples/Calibration Curves/RI Calibration Curve 2025 April.txt"
x_lim = [1e3, 1e8]  # default set to [1e3, 1e8]
y_lim = [-0.05, 1.1]  # default set to [-0.05, 1.2]
normalize_peak_range = [1e3, 1e4]  # default set to [1e4, 1e8]

baseline_method = 'flat'  # can change to 'flat' (1 range), 'linear'(2 ranges), or 'quadratic' (3 ranges)
baseline_ranges = [[1e7, 5e7]]  # MW ranges for baseline calculation

matplotlib.rcParams['font.family'] = 'Avenir Next LT Pro'
matplotlib.rcParams['font.size'] = 18
plt.figure(figsize=(16, 8))
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

    elif color_palette == 98:
        if number == 0:
            return "#298c8c"
        elif number == 1:
            return "#f1a226"

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
    # load x (MW) and raw y for this trace
    x_a = np.delete(data_array_RI[:, 1], [0, 1], 0).astype(float)
    x_a = 10 ** x_a
    y_a = np.delete(data_array[:, index * 2 + 1], [0, 1], 0).astype(float)


    # 2) subtract local minimum in the normalization window
    min_y = min_of_y_within_range(x_a, y_a)
    y_a = y_a - min_y

    # 3) divide by local maximum in the peak window
    max_y = max_of_y_within_range(x_a, y_a)
    y_a = y_a / max_y

    return x_a, y_a


# your plotting setup …
for i in range(int(data_array.shape[1] / 2)):
    x, y = get_data(i)
    plt.plot(x, y, color(i), label=data_array[0, 2 * i])
    # optional: to see the fitted baseline itself, uncomment:
    # plt.plot(x, baseline, '--', color(i), alpha=0.5)

plt.xscale("log")

plt.xscale("log")
ax = plt.gca()
ax.set_xlim(x_lim)
ax.set_ylim(y_lim)
plt.xlabel("Molecular weight (g/mol)", fontsize=20, name='Avenir Next LT Pro', fontstyle='italic', fontweight='demi')
plt.ylabel("Normalization", fontsize=20, name='Avenir Next LT Pro', fontstyle='italic', fontweight='demi')
plt.legend()
plt.savefig(os.path.join(output_dir, os.path.splitext(os.path.basename(txt_file))[0] + '.png'), dpi=300)
#plt.show()
