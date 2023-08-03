import numpy as np
import matplotlib.pyplot as plt

bandwidth = np.array([18953.0, 19002.4, 18943.9, 18899.1, 18993.3])
bandwidth32 = np.array([51052.8, 54312.8, 47293.1, 47917.8, 49658.8])
bandwidth64 = np.array([50141.1, 46512.9, 45273.5, 59126.8, 43907.9])
bandwidth128 = np.array([40354.1, 41087.9, 39789.4, 40890.1, 41407.3])

# threads 128
bandwidth_guided = np.array([38798.0, 37148.6, 40148.9, 50632.9, 33791.0])
bandwidth_static = np.array([1185669.0, 1231355.3, 1080658.0, 1082401.0, 1120348.3])
bandwidth_dynamic = np.array([213.6, 230.6, 222.3, 231.8, 223.7])


dft = np.array([4.477748, 4.480523, 4.458912, 4.472387, 4.452438])
dft32 = np.array([0.143021, 0.144944, 0.144730, 0.145131, 0.144137])
dft64 = np.array([0.076717, 0.078026, 0.076500, 0.076694, 0.076794])
dft128 = np.array([0.069757, 0.066805, 0.070983, 0.070422, 0.074288])

bandwidth_mean = np.mean(bandwidth)
bandwidth32_mean = np.mean(bandwidth32)
bandwidth64_mean = np.mean(bandwidth64)
bandwidth128_mean = np.mean(bandwidth128)

bandwidth_static_mean = np.mean(bandwidth_static)
bandwidth_dynamic_mean = np.mean(bandwidth_dynamic)
bandwidth_guided_mean = np.mean(bandwidth_guided)

dft_mean = np.mean(dft)
dft32_mean = np.mean(dft32)
dft64_mean = np.mean(dft64)
dft128_mean = np.mean(dft128)

bandwidth_std = np.std(bandwidth)
bandwidth32_std = np.std(bandwidth32)
bandwidth64_std = np.std(bandwidth64)
bandwidth128_std = np.std(bandwidth128)

bandwidth_static_std = np.std(bandwidth_static)
bandwidth_dynamic_std = np.std(bandwidth_dynamic)
bandwidth_guided_std = np.std(bandwidth_guided)

dft_std = np.std(dft)
dft32_std = np.std(dft32)
dft64_std = np.std(dft64)
dft128_std = np.std(dft128)

schedule = ['static', 'dynamic', 'guided']
x_poss = np.arange(len(schedule))
CTEss = [bandwidth_static_mean, bandwidth_dynamic_mean, bandwidth_guided_mean]
errors = [bandwidth_static_std, bandwidth_dynamic_std, bandwidth_guided_std]

figs, axs = plt.subplots()
bar_colorss = ['tab:red', 'tab:blue', 'tab:green', 'tab:orange']
axs.bar(x_poss, CTEss, yerr=errors, align='center', alpha=0.5, color = bar_colors, ecolor='black', capsize=10)
axs.set_ylabel('Bandwidth (MB/s), log scale')
axs.set_xticks(x_poss)
axs.set_xticklabels(schedule)
axs.set_title('Bandwidth for different schedule ')
axs.yaxis.grid(True)
axs.set_yscale("log")

# Save the figure and show
plt.tight_layout()
plt.savefig('/home/yifei/Pictures/schedule.pdf')
plt.show()