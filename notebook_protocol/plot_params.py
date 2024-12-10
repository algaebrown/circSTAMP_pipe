import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
import seaborn as sns
from pathlib import Path
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams["image.cmap"] = "Dark2"
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Dark2.colors)

figdir = Path('/tscc/nfs/home/hsher/scratch/circular_fig/')