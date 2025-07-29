import numpy as np

voltages = [560]

# Example data (non-normal distributions)
infileoff = open("size_no_sc_combined.txt")
l  = infileoff.readlines()
l = [i[:-1].split(",") for i in l]

gain_list_off = []
gain_list_on = []

for i in voltages:
    gain_list_voltage = []
    for entry in l:
        if entry[0] == str(i):
            gain_list_voltage.append(int(entry[1]))
    gain_list_off.append(gain_list_voltage)

infileon = open("size_sc_combined.txt")
l  = infileon.readlines()
l = [i[:-1].split(",") for i in l]

for i in voltages:
    gain_list_voltage = []
    for entry in l:
        if entry[0] == str(i):
            gain_list_voltage.append(int(entry[1]))
    gain_list_on.append(gain_list_voltage)

group_off = np.array(gain_list_off[0])
group_on = np.array(gain_list_on[0])

# Step 1: Observed test statistic (difference in means)
obs_diff = np.mean(group_off) - np.mean(group_on)

# Step 2: Combine data
combined = np.concatenate([group_on, group_off])
n_a = len(group_off)
n_permutations = 10000
perm_diffs = []

# Step 3: Permutation loop
for _ in range(n_permutations):
    np.random.shuffle(combined)
    perm_group_a = combined[:n_a]
    perm_group_b = combined[n_a:]
    perm_diff = np.mean(perm_group_a) - np.mean(perm_group_b)
    perm_diffs.append(perm_diff)

# Step 4: Compute two-sided p-value
perm_diffs = np.array(perm_diffs)
p_value = np.mean(np.abs(perm_diffs) >= np.abs(obs_diff))

print(f"Observed difference in means: {obs_diff:.3f}")
print(f"P-value from permutation test: {p_value:.4f}")

import matplotlib.pyplot as plt

plt.hist(perm_diffs, bins=30, edgecolor='k', alpha=0.7)
plt.axvline(obs_diff, color='red', linestyle='--', label='Observed diff')
plt.axvline(-obs_diff, color='red', linestyle='--')
plt.title('Permutation Null Distribution')
plt.xlabel('Difference in Means')
plt.ylabel('Frequency')
plt.legend()
plt.show()
