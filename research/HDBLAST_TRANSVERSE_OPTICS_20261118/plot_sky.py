"""Render the two normalized angular patterns from the analytic quadrupole."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plot(path):
    lon = np.linspace(-np.pi, np.pi, 481)
    lat = np.linspace(-np.pi/2, np.pi/2, 241)
    az, alt = np.meshgrid(lon, lat)
    quad = np.cos(alt)**2 * np.cos(2*az)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), subplot_kw={'projection': 'mollweide'},
                             layout='constrained')
    for ax, values, cmap, lo, title in zip(axes, [quad, quad**2], ['RdBu_r','viridis'],
            [-1, 0], ['Shell area-distance change', 'Bulk volume focusing']):
        im = ax.pcolormesh(az, alt, values, cmap=cmap, vmin=lo, vmax=1, shading='auto', rasterized=True)
        ax.set_title(title, fontsize=14, pad=12)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(alpha=.22)
        bar = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=.067, pad=.09,
                           ticks=[-1,0,1] if lo < 0 else [0,.5,1])
        bar.set_label('First-order angular factor Q' if lo < 0 else 'Quadratic angular factor Q²')
    fig.savefig(path, dpi=160, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    plot(Path(__file__).with_name('directional_focusing.png'))
