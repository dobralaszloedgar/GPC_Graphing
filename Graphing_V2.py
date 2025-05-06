import os
import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np

# Initial variables
txt_file = "G:/Edgar Dobra/GPC Samples/2025 Spring/JBB Campaign 2/UV_Samples_All.txt"
output_dir = r"G:/Edgar Dobra/Polymer GPC Graphs/2025 Spring/JBB Campaign 2/UV Graphs"
os.makedirs(output_dir, exist_ok=True)
color_palette = 99
RI_calibration = "G:/Edgar Dobra/GPC Samples/Calibration Curves/UV Calibration Curve 2025 April.txt"
x_lim = [1e3, 1e8]
y_lim = [-0.05, 1.1]
normalize_peak_range = [1e3, 1e4]

baseline_method = 'automatic'
baseline_ranges = [[1e3, 1.01e3], [5e6, 5e7], [8e7, 1e8]]

# Setup plot styling
matplotlib.rcParams['font.family'] = 'Avenir Next LT Pro'
matplotlib.rcParams['font.size'] = 18

# Load data files
def load_data(file_path):
    with open(file_path, "r") as f:
        return np.array([line.strip("\n").split("\t") for line in f])

data_array = load_data(txt_file)
data_array_RI = load_data(RI_calibration)

# Color palette function
def color(number):
    palettes = {98: {0: "#298c8c", 1: "#f1a226"}, 99: {0: "#ff000d", 1: "#29A829"}}
    return palettes.get(color_palette, {}).get(number % 2, "#000000")



def fit_baseline(x, y, bl_ranges, method):
    # build (x̄,ȳ) for each range
    pts = []
    for lo, hi in bl_ranges:
        m = (x>=lo)&(x<=hi)
        if not m.any():
            raise ValueError(f"No data in baseline range [{lo},{hi}]")
        pts.append((x[m].mean(), y[m].mean()))
    xs, ys = map(np.array, zip(*pts))

    if method == 'flat':
        base = np.full_like(y, ys.mean())

    elif method in ('linear','quadratic'):
        deg = 1 if method=='linear' else 2
        c = np.polyfit(xs, ys, deg)
        base = np.polyval(c, x)

    elif method == 'logarithmic':
        # fit ȳ ≈ a·log₁₀(x) + b
        a, b = np.polyfit(np.log10(xs), ys, 1)
        base = a*np.log10(x) + b

    elif method == 'exponential':
        # fit ln(ȳ) ≈ B·x̄ + ln(A)  => base = A·e^{B x}
        if (ys<=0).any():
            raise ValueError("need positive ys for exponential")
        B, lnA = np.polyfit(xs, np.log(ys), 1)
        A = np.exp(lnA)
        base = A*np.exp(B*x)

    else:
        raise ValueError(f"unknown method {method}")

    return y-base, base


def auto_baseline_correction(x, y, bl_ranges):
    methods = ['flat','linear','quadratic','logarithmic','exponential']
    mask = np.zeros_like(x, bool)
    for lo,hi in bl_ranges:
        mask |= (x>=lo)&(x<=hi)

    best = None
    best_score = np.inf
    best_yc = None
    best_base = None

    for m in methods:
        try:
            yc, base = fit_baseline(x, y, bl_ranges, m)
            score = np.mean(np.abs(yc[mask]))
            if score < best_score:
                best_score = score
                best = m
                best_yc = yc
                best_base = base
        except:
            continue

    if best is None:
        raise RuntimeError("all baseline fits failed")
    print(f"→ automatic chose {best} (mean|res|={best_score:.2e})")
    return best_yc, best_base
# Enhanced baseline correction with multiple methods
def baseline_correction(x_mw, y, method='automatic'):
    # Create local copy of baseline ranges for safety
    local_ranges = [list(bl) for bl in baseline_ranges]
    ref_points = []
    method_requirements = {
        'none': 0,
        'flat': 1,
        'linear': 2,
        'quadratic': 3,
        'logarithmic': 3,
        'exponential': 2,
        'automatic': 3
    }

    required_ranges = method_requirements.get(method, 0)
    if len(baseline_ranges) != required_ranges:
        raise ValueError(f"{method} requires {required_ranges} baseline ranges")

    # Special handling for automatic flat correction
    if method == 'automatic':
        return auto_baseline_correction(x_mw, y, baseline_ranges)

    # Handle other methods
    for bl in baseline_ranges:
        mask = (x_mw >= bl[0]) & (x_mw <= bl[1])
        if not mask.any():
            raise ValueError(f"No data in baseline range {bl}")
        x_val = x_mw[mask].mean()
        y_val = y[mask].mean()
        if method == 'logarithmic':
            ref_points.append((np.log10(x_val), y_val))
        else:
            ref_points.append((x_val, y_val))

    xs, ys = zip(*ref_points)

    if method == 'none':
        return y, np.zeros_like(y)
    elif method == 'flat':
        return y - np.mean(ys), np.full_like(y, np.mean(ys))
    elif method == 'exponential':
        if any(y <= 0 for y in ys):
            raise ValueError("Exponential baseline requires positive range values")
        coeffs = np.polyfit(xs, np.log(ys), 1)
        base = np.exp(np.polyval(coeffs, x_mw))
    elif method == 'logarithmic':
        coeffs = np.polyfit(xs, ys, 2)
        base = np.polyval(coeffs, np.log10(x_mw))
    else:  # linear/quadratic
        deg = 1 if method == 'linear' else 2
        coeffs = np.polyfit(xs, ys, deg)
        base = np.polyval(coeffs, x_mw)

    return fit_baseline(x_mw, y, baseline_ranges, method)

# Normalization function
def normalize_peaks(x_mw, y):
    mask = (x_mw >= normalize_peak_range[0]) & (x_mw <= normalize_peak_range[1])
    peak_max = np.max(y[mask])
    return y / peak_max

# Data processing function
def process_data(index):
    # Get molecular weight axis (same for all traces)
    x_mw = 10 ** np.delete(data_array_RI[:, 1], [0, 1], 0).astype(float)

    # Get raw y-data for SPECIFIC TRACE (column index*2+1)
    y_raw = np.delete(data_array[:, index * 2 + 1], [0, 1], 0).astype(float)

    # INDEPENDENT baseline correction for this trace
    y_corrected, baseline = baseline_correction(x_mw, y_raw, method=baseline_method)

    # Normalize using THIS TRACE'S peak maximum
    y_normalized = normalize_peaks(x_mw, y_corrected)

    return x_mw, y_normalized

# Plot all traces
def create_individual_plots():
    total_graphs = int(data_array.shape[1] / 4)  # 4 columns per graph = 2 traces
    for graph_num in range(total_graphs):
        plt.figure(figsize=(16, 8))
        plt.subplots_adjust(bottom=0.19, left=0.19)
        all_peaks = []
        line_artists = []
        legend_labels = []

        # Process both traces in this 4-column group (columns 1 & 3)
        for group_position in [0, 2]:  # First and third columns in the 4-column group
            # Calculate absolute column index
            absolute_col = graph_num * 4 + group_position

            # Verify column exists in data
            if absolute_col >= data_array.shape[1]:
                continue

            # Get label from FIRST ROW of column 1 and 3 in the group
            label = data_array[0, absolute_col]

            # Get corresponding data columns (current col = label, next col = data)
            try:
                x_mw = 10 ** np.delete(data_array_RI[:, 1], [0, 1], 0).astype(float)
                y_raw = np.delete(data_array[:, absolute_col + 1], [0, 1], 0).astype(float)
            except IndexError:
                continue

            # Process data
            y_corrected, _ = baseline_correction(x_mw, y_raw, baseline_method)
            y_normalized = normalize_peaks(x_mw, y_corrected)

            # Plot and store legend info
            line, = plt.plot(x_mw, y_normalized,
                             color=color(graph_num * 2 + group_position // 2),
                             label=label)
            line_artists.append(line)
            legend_labels.append(label)

            # Track peaks for y-axis scaling
            mask = (x_mw >= x_lim[0]) & (x_mw <= x_lim[1])
            if mask.any():
                all_peaks.append(np.max(y_normalized[mask]))

        # Configure plot if we have data
        if line_artists:
            # Dynamic Y-axis
            y_max = (max(all_peaks) + 0.2) if all_peaks else 1.2
            plt.ylim(-0.05, y_max)

            # Labels and legend
            plt.xscale('log')
            plt.xlim(x_lim)
            plt.xlabel("Molecular weight (g/mol)",
                       fontsize=20, fontstyle='italic', fontweight='demi')
            plt.ylabel("Normalized Intensity",
                       fontsize=20, fontstyle='italic', fontweight='demi')
            plt.legend(line_artists, legend_labels,
                       loc='upper right',
                       frameon=False,
                       fontsize=14)

        # Save output
        output_path = os.path.join(output_dir, f"{graph_num + 1}_GT-GF_UV.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()


# Execute the plotting
create_individual_plots()
