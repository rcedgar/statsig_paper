import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

plotfn = "../figures/3di_freqs.svg"

# --- Data Initialization ---
nupper, nlower = 864550, 1084762
ntotal = nupper + nlower
freq_upper, freq_lower, letters = {}, {}, []

def add(letter, fu, fl):
    letters.append(letter)
    freq_upper[letter] = (fu * nupper) / ntotal
    freq_lower[letter] = (fl * nlower) / ntotal

add('A', 0.0302, 0.0169); add('C', 0.0121, 0.0233); add('D', 0.061, 0.0767)
add('E', 0.0135, 0.00945); add('F', 0.019, 0.00931); add('G', 0.0257, 0.00873)
add('H', 0.0274, 0.011); add('I', 0.0169, 0.00627); add('K', 0.0137, 0.00912)
add('L', 0.0165, 0.053); add('M', 0.00722, 0.0036); add('N', 0.0143, 0.0146)
add('P', 0.0382, 0.0399); add('Q', 0.0287, 0.0269); add('R', 0.0229, 0.0134)
add('S', 0.0155, 0.0482); add('T', 0.0163, 0.00461); add('V', 0.0299, 0.167)
add('W', 0.0174, 0.00743); add('Y', 0.0171, 0.00721)

freq_h, freq_s, freq_t, freq_l = {}, {}, {}, {}
def addss(letter, fh, fs, ft, fl):
    freq_h[letter], freq_s[letter], freq_t[letter], freq_l[letter] = fh, fs, ft, fl

addss('A', 0.0236, 0.3649, 0.0737, 0.5378); addss('C', 0.6611, 0.0084, 0.0835, 0.2470)
addss('D', 0.0062, 0.2883, 0.0522, 0.6534); addss('E', 0.0011, 0.6897, 0.0099, 0.2993)
addss('F', 0.0312, 0.3649, 0.0854, 0.5185); addss('G', 0.0217, 0.3487, 0.0623, 0.5673)
addss('H', 0.0209, 0.3020, 0.0724, 0.6047); addss('I', 0.0006, 0.8432, 0.0028, 0.1534)
addss('K', 0.0000, 0.7736, 0.0013, 0.2251); addss('L', 0.7388, 0.0021, 0.1186, 0.1405)
addss('M', 0.0000, 0.8531, 0.0007, 0.1462); addss('N', 0.5114, 0.0136, 0.1504, 0.3245)
addss('P', 0.1344, 0.0094, 0.3115, 0.5447); addss('Q', 0.3727, 0.0645, 0.1513, 0.4115)
addss('R', 0.3091, 0.0592, 0.1480, 0.4836); addss('S', 0.7561, 0.0012, 0.1089, 0.1338)
addss('T', 0.0078, 0.6252, 0.0257, 0.3413); addss('V', 0.7975, 0.0002, 0.0849, 0.1173)
addss('W', 0.0036, 0.6750, 0.0159, 0.3055); addss('Y', 0.0009, 0.7051, 0.0090, 0.2850)

# --- Logic & Scaling ---
total_upper = {k: freq_upper[k] + freq_lower[k] for k in letters}
sorted_keys = sorted(letters, key=lambda k: total_upper[k], reverse=True)

# Determine the actual height of the upper bars to remove white space
max_u = max(total_upper.values())
y_top_limit = max_u * 1.1  # 10% padding
# To make top twice as long as bottom, bottom limit must be half of top
y_bottom_limit = y_top_limit / 2.0
gap = y_top_limit * 0.15

fig, ax = plt.subplots(figsize=(8, 5))

color_hc = '#aec7e8'
color_lc = '#1f77b4'

# 1. Plot Upper Bars
u1 = [freq_upper[k] for k in sorted_keys]
u2 = [freq_lower[k] for k in sorted_keys]
ax.bar(sorted_keys, u2, bottom=[gap + v for v in u1], label='Low complexity', color=color_lc)
ax.bar(sorted_keys, u1, bottom=gap, label='High complexity', color=color_hc)

# 2. Plot Lower Bars (Scaled to fit in half the space)
# Since lower data sums to 1.0, we multiply by y_bottom_limit to fill that space
s_factor = y_bottom_limit 
l1 = [freq_h[k]*s_factor for k in sorted_keys]
l2 = [freq_s[k]*s_factor for k in sorted_keys]
l3 = [freq_t[k]*s_factor for k in sorted_keys]
l4 = [freq_l[k]*s_factor for k in sorted_keys]

color_h = "saddlebrown"
color_s = "chocolate"
color_t = "orange"
color_l = "peachpuff"
ax.bar(sorted_keys, [-v for v in l1], bottom=-gap, label='Helix', color=color_h)
ax.bar(sorted_keys, [-v for v in l2], bottom=[-gap-v for v in l1], label='Strand', color=color_s)
ax.bar(sorted_keys, [-v for v in l3], bottom=[-gap-v1-v2 for v1,v2 in zip(l1,l2)], label='Turn', color=color_t)
ax.bar(sorted_keys, [-v for v in l4], bottom=[-gap-v1-v2-v3 for v1,v2,v3 in zip(l1,l2,l3)], label='Loop', color=color_l)

# 3. Formatting
ax.set_ylim(-(y_bottom_limit + gap), (y_top_limit + gap))
ax.axhline(0, color='black', linewidth=1)
ax.set_xticks(range(len(sorted_keys)))
ax.set_xticklabels(sorted_keys, fontweight='bold')
ax.spines['bottom'].set_position('zero')

# 4. Correct Ticks
# Top Ticks: actual data values
u_ticks = [gap + i*(max_u/3) for i in range(4)]
ax.set_yticks(u_ticks, [f"{i*(max_u/3):.2f}" for i in range(4)], minor=False)

# Lower Ticks: Manually mapped 0 to 1
l_ticks = [-gap - (i*0.2 * s_factor) for i in range(6)]
ax.set_yticks(list(ax.get_yticks()) + l_ticks)
all_labels = [f"{i*(max_u/3):.2f}" for i in range(4)] + [f"{i*0.2:.1f}" for i in range(6)]
ax.set_yticklabels(all_labels)

# plt.title('Sequence Distribution (2x Scale) vs Structure Propensity (1x Scale)', pad=30)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
print(plotfn)
plt.savefig(plotfn)