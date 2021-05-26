import cProfile
import io
import logging
import os
import pathlib
import pstats
import sys
import numpy as np
import yaml
import pandas as pd
from pprint import pprint


def mkdir(init_directory):
    dirs = list(sorted(list(set(["/".join(init_directory.split("/")[: j + 1]) for j in range(len(init_directory.split("/")))]))))
    for directory in dirs:
        if not os.path.exists(directory):
            try:
                pathlib.Path(directory).mkdir(exist_ok=True)
            except FileNotFoundError:
                logging.error("Unable to find or create directory {}".format(directory))
                sys.exit()


def load_config_file(config_file="src/config.yaml"):
    return yaml.load(open(config_file), Loader=yaml.FullLoader)


def convert_bins_into_labels(bins):
    return [" - ".join([str(round(bins[j], 2)) + " - " + str(round(bins[j + 1], 2))]) for j, e in enumerate(bins) if j < len(bins) - 1]


def output_figure(f, name):
    print(name)
    for ext in ["png", "jpg"]:
        for res in [150, 300]:
            f.savefig(name + "_{}_DPI.{}".format(res, ext), dpi=res)


def show_values_on_bars(axs, i=0, fontsize=13, rotation=0):
    def _show_on_single_plot(ax):
        for p in ax.patches:
            print(p)
            _x = p.get_x() + p.get_width() / 2
            _y = p.get_y() + (p.get_height()) + 20
            if i == 0:
                value = "{:.0f}".format(p.get_height())
            if i == 2:
                value = "{:.2f}".format(p.get_height())

            if i == 3:
                value = "{:.3f}".format(p.get_height())
            ax.text(_x, _y, str(round(int(value), 0)), ha="center", fontsize=fontsize, rotation=rotation, color="black")

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)
