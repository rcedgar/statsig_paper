import matplotlib.pyplot as plt

def viridis_color(i, n=None):
	if n is None:
		assert i >= 0 and i <= 1
		return plt.get_cmap('viridis')(i)
	else:
		assert i >= 0 and i <= n
		return plt.get_cmap('viridis')(i/(n-1))
