# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Functions to plot panels for the paper about estimating_uncertainties_in_simulated_ENSO submitted to JAMES
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from inspect import stack as inspect__stack
from typing import Literal, Union
# cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# cmocean
import cmocean
# matplotlib
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Polygon
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
# numpy
from numpy import array as numpy__array
from numpy import isnan as numpy__isnan
from numpy import linspace as numpy__linspace
from numpy import ndarray as numpy__ndarray
from numpy import meshgrid as numpy__meshgrid
from numpy import zeros as numpy__zeros
# panda
import pandas
# scipy
from scipy.ndimage import gaussian_filter
# seaborn
import seaborn

# local package
from enso_change.plots.plot_tools import is_list_to_plot
from enso_change.tools.default_tools import none_to_default
from enso_change.tools.print_tools import print_fail
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default parameters for plots
# ---------------------------------------------------------------------------------------------------------------------#
# general values
default_plot: dict[str, list[float] | bool | float | int | str] = {
    # horizontal and vertical size of the plot
    "size_x": 4.,
    "size_y": 4.,
    # axes ticks
    "lab": [-1., 1.],
    "lab_side": "both",
    "tic": [-1., 1.],
    "nbr_minor": 1,
    # reference color / line / marker
    "color": "k",
    "linestyle": "-",
    "linewidth": 1.,
    "marker": "D",
    "size": 100.,
    # text
    "fontsize": 15.,
    "horizontalalignment": "left",
    "rotation": 0.,
    "text": "",
    "verticalalignment": "center",
    "weight": "normal",
    # plot_ax_title
    "title_column_h": "left",
    "title_column_v": "bottom",
    "title_column_x": 50.,
    "title_column_y": 10.,
    "title_row_h": "center",
    "title_row_v": "center",
    "title_row_x": 35.,
    "title_row_y": 50.,
    "title_x": 4.05,
    "title_y": 4.1,
    # plot_boxplot
    "boxplot_flier": True,
    "boxplot_flier_size": 5.,
    "boxplot_linestyle": "-",
    "boxplot_linewidth": 2,
    "boxplot_mean": True,
    "boxplot_mean_size": 9.,
    "boxplot_vert": True,
    "boxplot_zorder": 1,
    # plot_curve
    "curve_zorder": 2,
    # plot_heatmap
    "heatmap_annotate": False,
    "heatmap_format": "{x:.1f}",
    "heatmap_colorbar": "cmo.balance",
    "heatmap_colorbar_bool": True,
    "heatmap_colorbar_extend": "neither",
    "heatmap_colorbar_fraction": 0.15,
    "heatmap_colorbar_location": "right",
    "heatmap_colorbar_shrink": 1.0,
    # plot_marker
    "marker_position": "bottom",
    "marker_zorder": 3,
    # plot_shading
    "shading_alpha": 1.,
    "shading_color": "r",
    "shading_zorder": 1,
    # plot_violinplot
    "violinplot_alpha": 0.8,
    "violinplot_bw_method": "scott",
    "violinplot_flier": True,
    "violinplot_mean": True,
    "violinplot_median": False,
    "violinplot_orientation": "vertical",
    "violinplot_points": 100,
    "violinplot_quantiles": None,
    "violinplot_side": "both",
    "violinplot_zorder": 1,
}
# legend
default_legend = {
    "line": {
        "color": default_plot["color"],
        "linestyle": None,
        "linewidth": 2},
    "marker": {
        "edgecolor": default_plot["color"],
        "facecolor": default_plot["color"],
        "linestyle": default_plot["linestyle"],
        "linewidth": default_plot["linewidth"],
        "marker": None,
        "s": default_plot["size"]},
    "position": {
        "x": 5,
        "y": 95},
    "text": {
        "color": default_plot["color"],
        "fontsize": default_plot["fontsize"],
        "horizontalalignment": default_plot["horizontalalignment"],
        "verticalalignment": default_plot["verticalalignment"],
        "rotation": default_plot["rotation"]},
}
# map
default_map: dict[str, float | list[float] | ccrs.PlateCarree | str] = {
    "colorbar": "cmo.balance",
    "color": "gainsboro",
    "lat_lim": [-20., 20.],
    "lat_tic": [-15, 0, 15],
    "lon_lim": [130., 290.],
    "lon_tic": [140, 180, 220, 260],
    "projection": ccrs.PlateCarree(central_longitude=0),
    "size_x": 7.2,
    "size_y": 3.6,
}
default_region = {
    # "n30e": {"alpha": 1, "color": "k", "facecolor": "none", "linestyle": "-", "linewidth": 2,
    #          "latitude": [-5., 5.], "longitude": [210., 270.], "short_name": "N3", "zorder": 5},
    # "n34e": {"alpha": 1, "color": "k", "facecolor": "none", "linestyle": "-", "linewidth": 2,
    #          "latitude": [-5., 5.], "longitude": [190., 240.], "short_name": "N3.4", "zorder": 5},
    # "n40e": {"alpha": 1, "color": "k", "facecolor": "none", "linestyle": "-", "linewidth": 2,
    #          "latitude": [-5., 5.], "longitude": [160., 210.], "short_name": "N4", "zorder": 5},
    # "plew": {"alpha": 1., "color": "k", "facecolor": "none", "linestyle": "--", "linewidth": 2,
    #          "latitude": [-5., 5.], "longitude": [110., 180.], "short_name": "", "zorder": 5},
    # "plee": {"alpha": 1., "color": "k", "facecolor": "none", "linestyle": "-", "linewidth": 2,
    #          "latitude": [-5., 5.], "longitude": [180., 280.], "short_name": "", "zorder": 4},
    # --> test
    "n30e": {"alpha": 1, "color": "k", "facecolor": "none", "linestyle": "-", "linewidth": 2,
             "latitude": [-5., 5.], "longitude": [210., 270.], "short_name": "N3", "zorder": 2},
    "n34e": {"alpha": 1, "color": "limegreen", "facecolor": "none", "linestyle": ":", "linewidth": 2,
             "latitude": [-5., 5.], "longitude": [190., 240.], "short_name": "N3.4", "zorder": 5},
    "n40e": {"alpha": 1, "color": "k", "facecolor": "none", "linestyle": "-", "linewidth": 2,
             "latitude": [-5., 5.], "longitude": [160., 210.], "short_name": "N4", "zorder": 2},
    "n30s": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
             "latitude": [-7., 3.], "longitude": [210., 270.], "short_name": "N3S", "zorder": 3},
    "n43s": {"alpha": 1, "color": "limegreen", "facecolor": "none", "linestyle": ":", "linewidth": 2,
             "latitude": [-7., 3.], "longitude": [160., 270.], "short_name": "N3&4S", "zorder": 4},
    "peen": {"alpha": 1, "color": "limegreen", "facecolor": "none", "linestyle": "--", "linewidth": 2,
             "latitude": [2., 8.], "longitude": [230., 285.], "short_name": "EEP$_{N}$", "zorder": 3},
    "pees": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
             "latitude": [-4., 2.], "longitude": [230., 285.], "short_name": "EEP$_{S}$", "zorder": 2},
    "peqntx": {"alpha": 1, "color": "limegreen", "facecolor": "none", "linestyle": "-", "linewidth": 2,
               "latitude": [5., 15.], "longitude": [215., 275.], "short_name": "N3S", "zorder": 3},
    "peqspr": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
               "latitude": [-7., 3.], "longitude": [200., 280.], "short_name": "N3S", "zorder": 3},
    "peqstx": {"alpha": 1, "color": "limegreen", "facecolor": "none", "linestyle": "-", "linewidth": 2,
               "latitude": [-14., -4.], "longitude": [240., 290.], "short_name": "N3S", "zorder": 3},
    "p063pr": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
               "latitude": [-8., 4.], "longitude": [190., 280.], "short_name": "PEQS", "zorder": 3},
    "p074pr": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
               "latitude": [-6., 4.], "longitude": [210., 280.], "short_name": "PEQS", "zorder": 3},
    "p117pr": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
               "latitude": [-6., 4.], "longitude": [210., 275.], "short_name": "PEQS", "zorder": 3},
    "p120pr": {"alpha": 1, "color": "magenta", "facecolor": "none", "linestyle": "-", "linewidth": 3,
               "latitude": [-8., 3.], "longitude": [190., 280.], "short_name": "PEQS", "zorder": 3},
    "p043n": {"alpha": 1, "color": "limegreen", "facecolor": "none", "linestyle": "--", "linewidth": 2,
              "latitude": [2., 8.], "longitude": [228., 278.], "short_name": "EEP$_{N}$", "zorder": 3},
    "p043s": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
              "latitude": [-4., 2.], "longitude": [228., 278.], "short_name": "EEP$_{S}$", "zorder": 2},
    "p068n": {"alpha": 1, "color": "magenta", "facecolor": "none", "linestyle": "--", "linewidth": 4,
              "latitude": [2., 8.], "longitude": [229., 279.], "short_name": "EEP$_{N}$", "zorder": 3},
    "p068s": {"alpha": 1, "color": "magenta", "facecolor": "none", "linestyle": "-", "linewidth": 3,
              "latitude": [-4., 2.], "longitude": [229., 279.], "short_name": "EEP$_{S}$", "zorder": 2},
    "pcain": {"alpha": 1, "color": "limegreen", "facecolor": "none", "linestyle": "--", "linewidth": 2,
              "latitude": [5., 10.], "longitude": [210., 270.], "short_name": "EEP$_{N}$", "zorder": 3},
    "pcais": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
              "latitude": [-2.5, 2.5], "longitude": [210., 270.], "short_name": "EEP$_{S}$", "zorder": 2},
    "pgxw": {"alpha": 1, "color": "r", "facecolor": "none", "linestyle": "-", "linewidth": 2,
              "latitude": [-5, 5], "longitude": [140., 190.], "short_name": "WEP", "zorder": 2},
    "pgxe": {"alpha": 1, "color": "limegreen", "facecolor": "none", "linestyle": "--", "linewidth": 2,
             "latitude": [-5, 5], "longitude": [220., 280.], "short_name": "EEP", "zorder": 3},
    "plaw": {"alpha": 1, "color": "magenta", "facecolor": "none", "linestyle": "-", "linewidth": 3,
             "latitude": [-4, 4], "longitude": [135., 230.], "short_name": "WEP", "zorder": 2},
    "plae": {"alpha": 1, "color": "magenta", "facecolor": "none", "linestyle": "--", "linewidth": 4,
             "latitude": [-10, 10], "longitude": [230., 270.], "short_name": "EEP", "zorder": 3},
    "pl2w": {"alpha": 1, "color": "magenta", "facecolor": "none", "linestyle": "-", "linewidth": 3,
             "latitude": [-5, 5], "longitude": [135., 185.], "short_name": "WEP", "zorder": 2},
    "pl2e": {"alpha": 1, "color": "magenta", "facecolor": "none", "linestyle": "--", "linewidth": 4,
             "latitude": [-6, 6], "longitude": [230., 270.], "short_name": "EEP", "zorder": 3},
    # <--
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def _ax_axis_x(
        ax,
        fontsize: float = None,
        x_lab: list = None,
        x_lab_rot: float = None,
        x_lim: list = None,
        x_nam: str = None,
        x_nbr_minor: int = None,
        x_tic: list = None,
        **kwargs):
    """
    Define x-axis

    Inputs:
    -------
    :param ax: matplotlib Axes object
    :param fontsize: float, optional
        Size of the font; e.g., fontsize = 15
    :param x_lab: list, optional
        List of labels to write on x-axis
    :param x_lab_rot: float, optional
        Label rotation on x-axis
    :param x_lim: list, optional
        Bounds of x-axis
    :param x_nam: str, optional
        Name of the x-axis
    :param x_nbr_minor: int, optional
        Number of minor x tick(s) between major ticks
    :param x_tic: list, optional
        Tick positions on x-axis
    """
    # set default values if not given
    fontsize = none_to_default(fontsize, default_plot["fontsize"])
    x_lab = none_to_default(x_lab, default_plot["lab"])
    x_lab_rot = none_to_default(x_lab_rot, default_plot["rotation"])
    x_nam = none_to_default(x_nam, default_plot["text"])
    x_nbr_minor = none_to_default(x_nbr_minor, default_plot["nbr_minor"])
    x_tic = none_to_default(x_tic, default_plot["tic"])
    x_lim = none_to_default(x_lim, [min(x_tic), max(x_tic)])
    # major x ticks
    ax.set_xticks(x_tic, minor=False)
    # minor x ticks
    if isinstance(x_nbr_minor, int) is True and x_nbr_minor > 0:
        # put lx_nbr_minor minor x tick between major ticks
        dx = abs(x_tic[1] - x_tic[0]) / (x_nbr_minor + 1)
        minor = list()
        xx = min(x_tic)
        while xx > min(x_lim):
            xx -= dx
        while xx < max(x_lim):
            if xx not in x_tic:
                minor.append(xx)
            xx += dx
        ax.set_xticks(minor, minor=True)
    # tick labels
    if 0 < x_lab_rot < 90:
        ax.set_xticklabels(x_lab, ha="right")
        pad = -3
    else:
        ax.set_xticklabels(x_lab)
        pad = 2
    # bounds of x-axis
    ax.set_xlim(x_lim)
    ax.set_xlabel(x_nam, fontsize=fontsize)
    # params of x-axis
    ax.tick_params(axis="x", direction="out", which="both", labelsize=fontsize, bottom=True, top=True, left=False,
                   right=False, labelrotation=x_lab_rot, pad=pad)
    ax.xaxis.set_tick_params(length=6, width=1.5, which="major")
    ax.xaxis.set_tick_params(length=3, width=1.0, which="minor")


def _ax_axis_y(
        ax,
        fontsize: float = None,
        y_axe: str = None,
        y_col: str = None,
        y_lab: list = None,
        y_lim: list = None,
        y_nam: str = None,
        y_nbr_minor: int = None,
        y_tic: list = None,
        **kwargs):
    """
    Define y-axis

    Inputs:
    -------
    :param ax: matplotlib Axes object
    :param fontsize: float, optional
        Size of the font; e.g., fontsize = 15
    :param y_axe: str, optional
        Side where to plot the y-axis ('both', 'left', 'right')
    :param y_col: str optional
        Color of y_lab and y_nam
    :param y_lab: list, optional
        List of labels to write on y-axis
    :param y_lim: list, optional
        Bounds of y-axis
    :param y_nam: str, optional
        Name of the y-axis
    :param y_nbr_minor: int, optional
        Number of minor y tick(s) between major ticks
    :param y_tic: list, optional
        Tick positions on y-axis
    """
    # set default values if not given
    fontsize = none_to_default(fontsize, default_plot["fontsize"])
    y_axe = none_to_default(y_axe, default_plot["lab_side"])
    y_col = none_to_default(y_col, default_plot["color"])
    y_lab = none_to_default(y_lab, default_plot["lab"])
    y_nam = none_to_default(y_nam, default_plot["text"])
    y_nbr_minor = none_to_default(y_nbr_minor, default_plot["nbr_minor"])
    y_tic = none_to_default(y_tic, default_plot["tic"])
    y_lim = none_to_default(y_lim, [min(y_tic), max(y_tic)])
    # major y ticks
    ax.set_yticks(y_tic, minor=False)
    # minor y ticks
    if isinstance(y_nbr_minor, int) is True and y_nbr_minor > 0:
        # put y_nbr_minor minor y tick between major ticks
        dy = abs(y_tic[1] - y_tic[0]) / (y_nbr_minor + 1)
        minor = list()
        yy = min(y_tic)
        while yy > min(y_lim):
            yy -= dy
        while yy < max(y_lim):
            if yy not in y_tic:
                minor.append(yy)
            yy += dy
        ax.set_yticks(minor, minor=True)
    # tick labels
    ax.set_yticklabels(y_lab)
    # bounds of y-axis
    ax.set_ylim(y_lim)
    ax.set_ylabel(y_nam, fontsize=fontsize)
    # params of y-axis
    if y_axe == "both":
        ytic_left, ytic_right = True, True
        ylab_left, ylab_right = True, False
    elif y_axe == "left":
        ytic_left, ytic_right = True, False
        ylab_left, ylab_right = True, False
        ax.yaxis.set_label_position("left")
    else:
        # ytic_left, ytic_right = False, True
        # ylab_left, ylab_right = False, True
        # ax.yaxis.set_label_position("right")
        ytic_left, ytic_right = True, True
        ylab_left, ylab_right = False, True
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
    ax.tick_params(axis="y", colors=y_col, direction="out", which="both", labelsize=fontsize, bottom=False, top=False,
                   left=ytic_left, right=ytic_right, labelbottom=False, labeltop=False, labelleft=ylab_left,
                   labelright=ylab_right)
    ax.yaxis.set_tick_params(length=6, width=1.5, which="major")
    ax.yaxis.set_tick_params(length=3, width=1.0, which="minor")
    ax.yaxis.label.set_color(y_col)


def _ax_axis_y2(
        ax2,
        fontsize: float = None,
        y2_axe: str = None,
        y2_col: str = None,
        y2_lab: list = None,
        y2_lim: list = None,
        y2_nam: str = None,
        y2_nbr_minor: int = None,
        y2_tic: list = None,
        **kwargs):
    tmp_kwargs = copy__deepcopy(kwargs)
    tmp_kwargs.update({"fontsize": fontsize, "y_axe": y2_axe, "y_col": y2_col, "y_lab": y2_lab, "y_lim": y2_lim,
                      "y_nam": y2_nam, "y_nbr_minor": y2_nbr_minor, "y_tic": y2_tic})
    _ax_axis_y(ax2, **tmp_kwargs)


def _ax_legend(
        ax,
        legend_param: dict = None,
        legend_txt: list = None,
        x_size: float = None,
        **kwargs):
    """
    Plot the legend if a keyword is used

    Input:
    ------
    :param ax: matplotlib Axes object
    :param legend_param, dict, optional
        Dictionary with two nested levels [legend_name, legend_keys1, legend_keys2]
        The first level (legend_name) must correspond to the values in legend_txt
        The second level (legend_keys1) must be made of keywords supported in this code:
            line: definition of the line (color, linestyle, linewidth)
            marker: definition of the line (edgecolor, facecolor, linestyle, linewidth, s)
            position: horizontal / vertical distance (percentage) from the left / bottom axis (x, y)
            text: definition of the line (color, fontsize, horizontalalignment, verticalalignment, rotation)
    :param legend_txt: list, optional
        List of string to write on the figure
    :param x_size: float, optional
        Horizontal size of the plot
    """
    # set default values if not given
    x_size = none_to_default(x_size, default_plot["size_x"])
    # plot dimensions
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    dx = (x2 - x1) / 100.
    dy = (y2 - y1) / 100.
    # legend
    if isinstance(legend_txt, list):
        for txt in legend_txt:
            # get legend parameters
            if isinstance(legend_param, dict) is False or (
                    isinstance(legend_param, dict) is True and txt not in list(legend_param.keys())):
                dict_t = copy__deepcopy(default_legend)
            else:
                dict_t = dict((k, {}) for k in list(default_legend.keys()))
                for k1 in list(default_legend.keys()):
                    if k1 not in list(legend_param[txt].keys()):
                        dict_t[k1] = default_legend[k1]
                    else:
                        for k2 in list(default_legend[k1].keys()):
                            if k2 not in list(legend_param[txt][k1].keys()):
                                dict_t[k1][k2] = default_legend[k1][k2]
                            else:
                                dict_t[k1][k2] = legend_param[txt][k1][k2]
            # legend position
            xx1 = x1 + dict_t["position"]["x"] * dx
            xx2 = xx1 + 4 * dx * default_plot["size_x"] / x_size
            xx3 = xx2 + 4 * dx * default_plot["size_x"] / x_size
            xx4 = xx3 + 1 * dx * default_plot["size_x"] / x_size
            xx5 = copy__deepcopy(xx1)
            yy1 = y1 + dict_t["position"]["y"] * dy
            # plot line
            if dict_t["line"]["linestyle"] is not None:
                ax.plot([xx1, xx3], [yy1, yy1], clip_on=False, zorder=11, **dict_t["line"])
                xx5 = copy__deepcopy(xx4)
            # plot marker
            if dict_t["marker"]["marker"] is not None:
                ax.scatter([xx2], [yy1], clip_on=False, zorder=12, **dict_t["marker"])
                xx5 = copy__deepcopy(xx4)
            # plot txt
            ax.text(xx5, yy1, txt, **dict_t["text"])


def _ax_text(
        ax,
        text: list = None,
        text_bbox: list = None,
        text_c: list = None,
        text_fs: list = None,
        text_ha: list = None,
        text_r: list = None,
        text_va: list = None,
        text_w: list = None,
        text_x: list = None,
        text_y: list = None,
        **kwargs):
    """
    Plot text

    Input:
    ------
    :param ax: matplotlib Axes object
    :param text, list, optional
        List of text to write; e.g., text = ['a', 'b', 'c']
    :param text_bbox, list, optional
        List of background for the text; e.g., text = [dict(facecolor='red', alpha=0.5), None, None]
    :param text_c, list, optional
        List of color for text to write; e.g., text_c = ['k', 'r', 'limegreen']
    :param text_fs, list, optional
        List of fontsize for text to write; e.g., text_fs = [12, 15, 20]
    :param text_ha, list, optional
        List of horizontal alignment for text to write; e.g., text_ha = ['center', 'left', 'right']
    :param text_r, list, optional
        List of rotation for text to write; e.g., text_r = [0, 45, 90]
    :param text_va, list, optional
        List of vertical alignment for text to write; e.g., text_va = ['bottom', 'center', 'top']
    param text_w, list, optional
        List of font weight for text to write; e.g., text_x = ['regular', 'bold', 'regular']
        Choice: {a numeric value in range 0-1000, 'ultralight', 'light', 'normal', 'regular', 'book', 'medium', 'roman',
                 'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'}
    :param text_x, list, optional
        List of horizontal position for text to write (using horizontal axis); e.g., text_x = [1, 2, 3]
    :param text_y, list, optional
        List of vertical position for text to write (using left vertical axis); e.g., text_y = [1, 2, 3]
    """
    if isinstance(text, list) and len(text) > 0 and isinstance(text_x, list) and len(text_x) == len(text) and \
            isinstance(text_y, list) and len(text_y) == len(text):
        if not isinstance(text_bbox, list) or (isinstance(text_bbox, list) and len(text_bbox) != len(text)):
            text_bbox = [None] * len(text)
        if not isinstance(text_c, list) or (isinstance(text_c, list) and len(text_c) != len(text)):
            text_c = [default_plot["color"]] * len(text)
        if not isinstance(text_fs, list) or (isinstance(text_fs, list) and len(text_fs) != len(text)):
            text_fs = [default_plot["fontsize"]] * len(text)
        if not isinstance(text_ha, list) or (isinstance(text_ha, list) and len(text_ha) != len(text)):
            text_ha = [default_plot["horizontalalignment"]] * len(text)
        if not isinstance(text_r, list) or (isinstance(text_r, list) and len(text_r) != len(text)):
            text_r = [default_plot["rotation"]] * len(text)
        if not isinstance(text_va, list) or (isinstance(text_va, list) and len(text_va) != len(text)):
            text_va = [default_plot["verticalalignment"]] * len(text)
        if not isinstance(text_w, list) or (isinstance(text_w, list) and len(text_w) != len(text)):
            text_w = [default_plot["weight"]] * len(text)
        for t, b, c, f, h, r, v, w, x, y in zip(
                text, text_bbox, text_c, text_fs, text_ha, text_r, text_va, text_w, text_x, text_y):
            ax.text(x, y, t, bbox=b, color=c, fontsize=f, ha=h, rotation=r, va=v, weight=w)


def _ax_titles(
        ax,
        fontsize: float = None,
        title: str = None,
        title_col: str = None,
        title_col_fs: float = None,
        title_col_h: Literal["left", "right", "center"] = None,
        title_col_v: Literal["top", "bottom", "center"] = None,
        title_col_x: float = None,
        title_col_y: float = None,
        title_row: str = None,
        title_row_fs: float = None,
        title_row_h: Literal["left", "right", "center"] = None,
        title_row_v: Literal["top", "bottom", "center"] = None,
        title_row_x: float = None,
        title_row_y: float = None,
        title_x: float = None,
        title_y: float = None,
        x_size: float = None,
        y_size: float = None,
        **kwargs):
    """
    Plot column and row names

    Input:
    ------
    :param ax: matplotlib Axes object
    :param fontsize: float
        Size of the font
    :param title: str, optional
        Name of the panel; e.g. title = 'a'
    :param title_col: str, optional
        Name of the column; e.g. title_col = 'sst'
    :param title_col_fs: float
        Size of the font of the column name; e.g., title_col_fs = 15
    :param title_col_h: {"left", "right", "center"}, optional
        Horizontal alignment of the column name; e.g., title_col_h = 'left'
    :param title_col_v, {"top", "bottom", "center"}, optional
        Vertical alignment of the column name; e.g., title_col_v = 'bottom'
    :param title_col_x: float, optional
        Horizontal position of the column name (percentage of horizontal axis); e.g., title_col_x = 50
    :param title_col_y: float, optional
        Vertical position of the column name (percentage of vertical axis); e.g., title_col_y = 10
    :param title_row: str, optional
        Name of the row; e.g. title_row = 'sst'
    :param title_row_fs: float
        Size of the font of the row name; e.g., title_row_fs = 15
    :param title_row_h: {"left", "right", "center"}, optional
        Horizontal alignment of the row name; e.g., title_row_h = 'center'
    :param title_row_v, {"top", "bottom", "center"}, optional
        Vertical alignment of the row name; e.g., title_row_v = 'center'
    :param title_row_x: float, optional
        Horizontal position of the row name (percentage of horizontal axis); e.g., title_row_x = 30
    :param title_row_y: float, optional
        Vertical position of the row name (percentage of vertical axis); e.g., title_row_y = 50
    :param title_x: float, optional
        Horizontal position of the panel name (percentage of horizontal axis); e.g., title_x = 4
    :param title_y: float, optional
        Vertical position of the panel name (percentage of vertical axis); e.g., title_y = 4
    :param x_size: float, optional
        Horizontal size of the plot; e.g., x_size = 8
    :param y_size: float, optional; e.g., y_size = 8
        Vertical size of the plot
    """
    # set default values if not given
    fontsize = none_to_default(fontsize, default_plot["fontsize"])
    title = none_to_default(title, default_plot["text"])
    title_col = none_to_default(title_col, default_plot["text"])
    title_col_fs = none_to_default(title_col_fs, fontsize)
    title_col_h = none_to_default(title_col_h, default_plot["title_column_h"])
    title_col_v = none_to_default(title_col_v, default_plot["title_column_v"])
    title_col_x = none_to_default(title_col_x, default_plot["title_column_x"])
    title_col_y = none_to_default(title_col_y, default_plot["title_column_y"])
    title_row = none_to_default(title_row, default_plot["text"])
    title_row_fs = none_to_default(title_row_fs, fontsize)
    title_row_h = none_to_default(title_row_h, default_plot["title_row_h"])
    title_row_v = none_to_default(title_row_v, default_plot["title_row_v"])
    title_row_x = none_to_default(title_row_x, default_plot["title_row_x"])
    title_row_y = none_to_default(title_row_y, default_plot["title_row_y"])
    title_x = none_to_default(title_x, default_plot["title_x"])
    title_y = none_to_default(title_y, default_plot["title_y"])
    x_size = none_to_default(x_size, default_plot["size_x"])
    y_size = none_to_default(y_size, default_plot["size_y"])
    # plot dimensions
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    dx = (x2 - x1) / 100.
    dy = (y2 - y1) / 100.
    # panel title
    if len(title) in [1, 2]:
        yy = y2 - title_y * dy * default_plot["size_y"] / y_size
        for k in list(range(len(title))):
            xx = x1 + (1 + k) * title_x * dx * default_plot["size_x"] / x_size
            ax.plot([xx], [yy], markersize=17, color="darkgrey", marker="s", fillstyle="full",
                    markeredgecolor="darkgrey", markeredgewidth=0, zorder=11)
        xx = x1 + (1 + (len(title) - 1) / 2) * title_x * dx * default_plot["size_x"] / x_size
        ax.text(xx, yy, title, fontsize=fontsize, fontweight="bold", zorder=12, ha="center", va="center")
    else:
        ax.set_title(title, fontsize=fontsize, y=1 + 0.02 * default_plot["size_y"] / y_size, loc="right")
    # column title
    xx = x1 + title_col_x * dx
    yy = y2 + title_col_y * dy * default_plot["size_y"] / y_size
    ax.text(xx, yy, title_col, fontsize=title_col_fs, color="k", weight="bold", ha=title_col_h, va=title_col_v)
    # row title
    xx = x1 - title_row_x * dx * default_plot["size_x"] / x_size
    yy = y1 + title_row_y * dy
    ax.text(xx, yy, title_row, fontsize=title_row_fs, color="k", weight="bold", ha=title_row_h, va=title_row_v,
            rotation=90)


def _plot_annotate(
        ax,
        ann_arr: list = None,
        ann_coords1: list = None,
        ann_coords2: list = None,
        ann_clip: list = None,
        ann_t: list = None,
        ann_t_c: list = None,
        ann_t_fs: list = None,
        ann_t_ha: list = None,
        ann_t_r: list = None,
        ann_t_va: list = None,
        ann_xy1: list = None,
        ann_xy2: list = None,
        ann_z: list = None,
        **kwargs):
    if isinstance(ann_xy1, list) and isinstance(ann_xy2, list) and len(ann_xy1) == len(ann_xy2):
        # arrow
        if ann_coords1 is None or (isinstance(ann_coords1, list) and len(ann_coords1) != len(ann_xy1)):
            ann_coords1 = ["data"] * len(ann_xy1)
        if ann_coords2 is None or (isinstance(ann_coords2, list) and len(ann_coords2) != len(ann_xy1)):
            ann_coords2 = ["data"] * len(ann_xy1)
        if ann_arr is None or (isinstance(ann_arr, list) and len(ann_arr) != len(ann_xy1)):
            ann_arr = [dict(facecolor="k", shrink=0)] * len(ann_xy1)
        if ann_clip is None or (isinstance(ann_clip, list) and len(ann_clip) != len(ann_xy1)):
            ann_clip = [False] * len(ann_xy1)
        if ann_z is None or (isinstance(ann_z, list) and len(ann_z) != len(ann_xy1)):
            ann_z = [1] * len(ann_xy1)
        # text
        if not isinstance(ann_t, list) or (isinstance(ann_t, list) and len(ann_t) != len(ann_xy1)):
            ann_t = [""] * len(ann_xy1)
        if ann_t_c is None or (isinstance(ann_t_c, list) and len(ann_t_c) != len(ann_xy1)):
            ann_t_c = [default_plot["color"]] * len(ann_xy1)
        if ann_t_fs is None or (isinstance(ann_t_fs, list) and len(ann_t_fs) != len(ann_xy1)):
            ann_t_fs = [default_plot["fontsize"]] * len(ann_xy1)
        if ann_t_ha is None or (isinstance(ann_t_ha, list) and len(ann_t_ha) != len(ann_xy1)):
            ann_t_ha = [default_plot["horizontalalignment"]] * len(ann_xy1)
        if ann_t_r is None or (isinstance(ann_t_r, list) and len(ann_t_r) != len(ann_xy1)):
            ann_t_r = [default_plot["rotation"]] * len(ann_xy1)
        if ann_t_va is None or (isinstance(ann_t_va, list) and len(ann_t_va) != len(ann_xy1)):
            ann_t_va = [default_plot["verticalalignment"]] * len(ann_xy1)
        # -- plot
        l_to_plot = zip(ann_t, ann_t_c, ann_t_fs, ann_t_ha, ann_t_r, ann_t_va, ann_xy1, ann_xy2, ann_coords1,
                        ann_coords2, ann_arr, ann_clip, ann_z)
        for tt, tc, tf, th, tr, tv, xy1, xy2, co1, co2, prop, clip, zz in l_to_plot:
            ax.annotate(tt, color=tc, fontsize=tf, ha=th, rotation=tr, va=tv, xy=xy1, xycoords=co1, xytext=xy2,
                        textcoords=co2, arrowprops=prop, annotation_clip=clip, zorder=zz)


def _plot_boxplot(
        ax,
        box_c: list = None,
        box_flie: list = None,
        box_fs: list = None,
        box_ls: list = None,
        box_lw: list = None,
        box_mean: list = None,
        box_ms: list = None,
        box_vert: list = None,
        box_w: list = None,
        box_x: list = None,
        box_y: list = None,
        box_z: list = None,
        **kwargs):
    """
    Plot boxplots

    Input:
    ------
    :param ax: matplotlib Axes object
    :param box_c: list, optional
        List of color for boxplots; e.g., box_c = ['k', 'b']
    :param box_flie: list, optional
        List of boolean outliers for boxplots, if True, show the outliers; e.g., box_flie = [True, False]
    :param box_fs: list, optional
        List of size of outlier marker for boxplots; e.g., box_fs = [2, 3]
    :param box_mean: list, optional
        List of boolean means for boxplots, if True, show the means; e.g., box_mean = [True, False]
    :param box_ms: list, optional
        List of size of mean marker for boxplots; e.g., box_ms = [14, 15]
    :param box_vert: list, optional
        List of boolean orientations for boxplots, if True, draws vertical, else horizontal;
        e.g., box_vert = [True, False]
    :param box_w: list, optional
        List of width for boxplots; e.g., box_w = [0.5, 0.5]
    :param box_x: list, optional
        List of x values for boxplots; e.g., box_x = [1, 2]
    :param box_y: list, optional
        List of y values for boxplots; e.g., box_y = [ [1, 2, 3], [1, 2, 3] ]
    :param box_z: list, optional
        List of zorder for boxplots (drawing order); e.g., box_z = [2, 4]
    """
    if isinstance(box_x, list) is True and isinstance(box_y, list) is True and len(box_x) == len(box_y):
        # control boxplot parameters
        if isinstance(box_c, list) is False or (isinstance(box_c, list) is True and len(box_c) != len(box_x)):
            box_c = [default_plot["color"]] * len(box_x)
        if isinstance(box_flie, list) is False or (isinstance(box_flie, list) is True and len(box_flie) != len(box_x)):
            box_flie = [default_plot["boxplot_flier"]] * len(box_x)
        if isinstance(box_fs, list) is False or (isinstance(box_fs, list) is True and len(box_fs) != len(box_x)):
            box_fs = [default_plot["boxplot_flier_size"]] * len(box_x)
        if isinstance(box_ls, list) is False or (isinstance(box_ls, list) is True and len(box_ls) != len(box_x)):
            box_ls = [default_plot["boxplot_linestyle"]] * len(box_x)
        if isinstance(box_lw, list) is False or (isinstance(box_lw, list) is True and len(box_lw) != len(box_x)):
            box_lw = [default_plot["boxplot_linewidth"]] * len(box_x)
        if isinstance(box_mean, list) is False or (isinstance(box_mean, list) is True and len(box_mean) != len(box_x)):
            box_mean = [default_plot["boxplot_mean"]] * len(box_x)
        if isinstance(box_ms, list) is False or (isinstance(box_ms, list) is True and len(box_ms) != len(box_x)):
            box_ms = [default_plot["boxplot_mean_size"]] * len(box_x)
        if isinstance(box_vert, list) is False or (isinstance(box_vert, list) is True and len(box_vert) != len(box_x)):
            box_vert = [default_plot["boxplot_vert"]] * len(box_x)
        if isinstance(box_w, list) is False or (isinstance(box_w, list) is True and len(box_w) != len(box_x)):
            dx = 1e20
            for k1 in range(len(box_x)):
                for k2 in range(k1 + 1, len(box_x)):
                    if abs(box_x[k1] - box_x[k2]) < dx:
                        dx = abs(box_x[k1] - box_x[k2])
            if dx == 1e20:
                dx = 1
            box_w = [dx * 0.7] * len(box_vert)
        if isinstance(box_z, list) is False or (isinstance(box_z, list) is True and len(box_z) != len(box_x)):
            box_z = [default_plot["boxplot_zorder"]] * len(box_x)
        # plot boxplot
        for c, f, fs, ls, lw, m, ms, v, w, x, y, z in zip(
                box_c, box_flie, box_fs, box_ls, box_lw, box_mean, box_ms, box_vert, box_w, box_x, box_y, box_z):
            dict_t = {
                "boxprops": dict(linestyle=ls, linewidth=lw, color=c),
                "capprops": dict(linestyle=ls, linewidth=lw, color=c),
                "flierprops": dict(marker="o", markersize=fs, markeredgecolor=c, markerfacecolor=c, markeredgewidth=0),
                "meanprops": dict(marker="D", markersize=ms, markeredgecolor=c, markerfacecolor=c, markeredgewidth=0),
                "medianprops": dict(linestyle=ls, linewidth=lw, color=c),
                "whiskerprops": dict(linestyle=ls, linewidth=lw, color=c)}
            ax.boxplot(y, positions=[x], whis=[5, 95], widths=w, labels=[""], showmeans=m, showfliers=f, vert=v,
                       zorder=z, **dict_t)
        # redo axes, the boxplot function causes troubles
        _ax_axis_x(ax, **kwargs)
        _ax_axis_y(ax, **kwargs)


def _plot_curve(
        ax,
        cur_c: list = None,
        cur_ls: list = None,
        cur_lw: list = None,
        cur_on: list = None,
        cur_x: list = None,
        cur_y: list = None,
        cur_z: list = None,
        **kwargs):
    """
    Plot curves

    Input:
    ------
    :param ax: matplotlib Axes object
    :param cur_c: list, optional
        List of colors for curves; e.g., cur_c = ['k', 'b']
    :param cur_ls: list, optional
        List of linestyle for curves; e.g., cur_ls = ['-', ':']
    :param cur_lw: list, optional
        List of linewidth for curves; e.g., cur_lw = [0.5, 0.5]
    :param cur_on: list, optional
        List of clip_on for curves; e.g., cur_on = [True, False]
    :param cur_x: list, optional
        List of x values for curves; e.g., cur_x = [ [1, 2, 3], [1, 2, 3] ]
    :param cur_y: list, optional
        List of y values for curves; e.g., cur_y = [ [10, 2, 10], [3, 6, 9] ]
    :param cur_z: list, optional
        List of zorder for curves (drawing order); e.g., cur_z = [2, 4]
    """
    if isinstance(cur_x, list) is True and isinstance(cur_y, list) is True and len(cur_x) == len(cur_y):
        # control curve parameters
        if cur_c is None or (isinstance(cur_c, list) is True and len(cur_c) != len(cur_x)):
            cur_c = [default_plot["color"]] * len(cur_x)
        if cur_ls is None or (isinstance(cur_ls, list) is True and len(cur_ls) != len(cur_x)):
            cur_ls = [default_plot["linestyle"]] * len(cur_x)
        if cur_lw is None or (isinstance(cur_lw, list) is True and len(cur_lw) != len(cur_x)):
            cur_lw = [default_plot["linewidth"]] * len(cur_x)
        if cur_on is None or (isinstance(cur_on, list) is True and len(cur_on) != len(cur_x)):
            cur_on = [True] * len(cur_x)
        if cur_z is None or (isinstance(cur_z, list) is True and len(cur_z) != len(cur_x)):
            cur_z = [default_plot["curve_zorder"]] * len(cur_x)
        # plot curves
        for c, ls, lw, on, x, y, z in zip(cur_c, cur_ls, cur_lw, cur_on, cur_x, cur_y, cur_z):
            ax.plot(x, y, clip_on=on, color=c, ls=ls, lw=lw, zorder=z)


def _plot_curve2(
        ax2,
        cur2_c: list = None,
        cur2_ls: list = None,
        cur2_lw: list = None,
        cur2_on: list = None,
        cur2_x: list = None,
        cur2_y: list = None,
        cur2_z: list = None,
        **kwargs):
    tmp_kwargs = {}
    l1 = ["cur_c", "cur_ls", "cur_lw", "cur_on", "cur_x", "cur_y", "cur_z"]
    l2 = [cur2_c, cur2_ls, cur2_lw, cur2_on, cur2_x, cur2_y, cur2_z]
    for k1, k2 in zip(l1, l2):
        if k2 is not None:
            tmp_kwargs[k1] = k2
    _plot_curve(ax2, **tmp_kwargs)


def _plot_density(
        ax, den_c: list = None, den_x: list = None, den_y: list = None, den_z: list = None, **kwargs):
    if isinstance(den_x, list) is True and isinstance(den_y, list) is True and len(den_x) == len(den_y):
        # control marker parameters
        if den_c is None or (isinstance(den_c, list) is True and len(den_c) != len(den_x)):
            den_c = [default_plot["color"]] * len(den_x)
        if den_z is None or (isinstance(den_z, list) is True and len(den_z) != len(den_x)):
            den_z = [1] * len(den_x)
        # plot using seaborn
        for c, x, y, z in zip(den_c, den_x, den_y, den_z):
            dd = pandas.DataFrame(numpy__array([x, y]).T, columns=['x', 'y'])
            seaborn.kdeplot(data=dd, color=c, levels=[0.2, 0.4, 0.6, 0.8], ax=ax, x='x', y='y', zorder=z)


def _plot_heatmap(
        ax,
        fontsize: float | int = None,
        hea_a_bool: bool = None,
        hea_a_color: str = None,
        hea_a_font: float | int = None,
        hea_a_form: str = None,
        hea_c: str = None,
        hea_cb: bool = None,
        hea_cl: Literal["left", "right", "top", "bottom"] = None,
        hea_cw: list[float | int] = None,
        hea_x: numpy__ndarray = None,
        hea_y: numpy__ndarray = None,
        hea_z: numpy__ndarray = None,
        hea_zv: numpy__ndarray = None,
        z_lab: list = None,
        z_nam: str = "",
        z_tic: list = None,
        z_tic_bool: bool = True,
        **kwargs):
    """
    Plot curves

    Inputs:
    -------
    :param ax: matplotlib Axes object
    :param fontsize: float or int, optional
        Size of the label font; e.g., fontsize = 15.
    :param hea_a_bool: bool, optional
        True to annotate the heatmap; e.g., hea_a_bool = True
    :param hea_a_color: str, optional
        Colors for annotation; e.g., hea_a_color = 'k', 'b'
    :param hea_a_font: float or int, optional
        Size of the annotation font; e.g., hea_a_font = 15.
    :param hea_a_form: float or int, optional
        The format of the annotations inside the heatmap; e.g., hea_a_form = '{x:.2f}'
    :param hea_c: str, optional
        Name of a colorbar; e.g., hea_c = 'cmo.balance'
    :param hea_cb: bool, optional
        True to plot the colorbar; e.g., hea_cb = True
    :param hea_cl: str, optional
        The location, relative to the parent axes, where the colorbar axes is created; e.g., hea_cl = 'right'
        The following options are available: 'left', 'right', 'top', 'bottom'
    :param hea_cw: list[float | int], optional
        List of values between which the colobar will be set to white
    :param hea_z: numpy__ndarray, optional
        The heat values over which the contour is drawn, the color-mapping is controlled by cmap
    :param hea_zv: numpy__ndarray, optional
        The heat values for annotation.
        Might be different from hea_z if the color must be different from the value (e.g., color only if significant)
    :param z_lab: list, optional
        Tick labels for colorbar; e.g., z_lab = ['1', '2']
    :param z_nam: str, optional
        Label for the colorbar
        Default is ''
    :param x_size: float, optional
        Horizontal size of the plot
    :param y_size: float, optional
        Vertical size of the plot
    :param z_tic: list, optional
        Ticks for colorbar; e.g., z_lab = [1, 2]
    """
    if isinstance(hea_z, numpy__ndarray):
        # set default values if not given
        fontsize = none_to_default(fontsize, default_plot["fontsize"])
        hea_a_bool = none_to_default(hea_a_bool, default_plot["heatmap_annotate"])
        hea_a_color = none_to_default(hea_a_color, default_plot["color"])
        hea_a_font = none_to_default(hea_a_font, default_plot["fontsize"])
        hea_a_form = none_to_default(hea_a_form, default_plot["heatmap_format"])
        hea_c = none_to_default(hea_c, default_plot["heatmap_colorbar"])
        hea_cb = none_to_default(hea_cb, default_plot["heatmap_colorbar_bool"])
        hea_cl = none_to_default(hea_cl, default_plot["heatmap_colorbar_location"])
        z_lab = none_to_default(z_lab, default_plot["lab"])
        z_tic = none_to_default(z_tic, default_plot["tic"])
        # colorbar levels
        if isinstance(z_tic_bool, bool) and z_tic_bool:
            dif, mul = round(float(z_tic[1] - z_tic[0]), 5), 5
            delta = round(dif / mul, 5)
            lev = [round(k1 + k2 * delta, 5) for k1 in z_tic[:-1] for k2 in range(mul)] + [z_tic[-1]]
        else:
            lev = copy__deepcopy(z_tic)
            if isinstance(hea_c, list) is True and len(lev) == len(hea_c):
                lev = []
                for k1, (k2, k3) in enumerate(zip(z_tic[:-1], z_tic[1:])):
                    if k1 == 0:
                        lev += [z_tic[0] - (k3 - k2) / 2]
                    lev += [k2 + (k3 - k2) / 2]
                    if k1 == len(z_tic) - 2:
                        lev += [z_tic[-1] + (k3 - k2) / 2]
        if isinstance(hea_c, list) and len(lev) != len(hea_c) + 1:
            lev = ticker.MaxNLocator(nbins=len(lev) + 1).tick_values(min(lev), max(lev))
        # shading
        if hea_c == "correlation":
            cmap = pyplot.get_cmap("cmo.balance")
            cmap = cmap(numpy__linspace(0.15, 0.85, len(lev) - 1))
            cmap = ListedColormap(cmap)
        elif isinstance(hea_c, list):
            cmap = (matplotlib.colors.ListedColormap(hea_c))
        else:
            cmap = pyplot.get_cmap(hea_c)
        #     cmap = cmap(numpy__linspace(0, 1, len(lev) - 1))
        # if isinstance(hea_cw, list) is True and len(hea_cw) == 2:
        #     for ii in list(range(len(cmap.colors))):
        #         if lev[ii] >= hea_cw[0] and lev[ii+1] <= hea_cw[1]:
        #             cmap.colors[ii] = [0, 0, 0, 0]
        cmap.set_bad("white")
        norm = BoundaryNorm(lev, ncolors=cmap.N, clip=True)
        im = ax.pcolormesh(hea_x, hea_y, hea_z, cmap=cmap, norm=norm)
        if isinstance(hea_cb, bool) and hea_cb:
            x1, x2 = ax.get_position().x0, ax.get_position().x1
            y1, y2 = ax.get_position().y0, ax.get_position().y1
            if hea_cl in ["bottom", "top"]:
                cax = pyplot.axes((x1, y1 - (y2 - y1) / 3.5, x2 - x1, (y2 - y1) / 15))
            elif hea_cl in ["left", "right"]:
                cax = pyplot.axes((x2 + (x2 - x1) / 10, y1, (x2 - x1) / 30, y2 - y1))
            else:
                cax = pyplot.axes((x1 + 2 * (x2 - x1) / 3, y1 - 2.3 * (y2 - y1), (x2 - x1) / 3, (y2 - y1) / 5))
            orientation = "horizontal" if hea_cl in ["bottom", "custom", "top"] else "vertical"
            cbar = pyplot.colorbar(im, cax=cax, orientation=orientation, ticks=z_tic, label=z_lab, pad=0.3)
            cbar.ax.tick_params(labelsize=fontsize)
            label_pad = 2 if hea_cl in ["bottom", "custom", "top"] else 0
            label_rotation = 0 if hea_cl in ["bottom", "custom", "top"] else 90
            cbar.set_label(z_nam, fontsize=fontsize, labelpad=label_pad, rotation=label_rotation)
        # annotate heatmap
        if hea_a_bool:
            # set text formatter
            text_format = ticker.StrMethodFormatter(hea_a_form)
            # loop over the data and create a `Text` for each "pixel"
            values = copy__deepcopy(hea_z) if hea_zv is None else copy__deepcopy(hea_zv)
            for i in range(values.shape[0]):
                for j in range(values.shape[1]):
                    if isinstance(values[i, j], (float, int)) is True and values[i, j] == -9999:
                        ax.text(j + 0.5, i + 0.5, "N/A", color="k", horizontalalignment="center",
                                size=hea_a_font, verticalalignment="center")
                    elif isinstance(values[i, j], (float, int)) and not numpy__isnan(values[i, j]):
                        im.axes.text(j + 0.5, i + 0.5, text_format(round(float(values[i, j]), 1), None),
                                     color=hea_a_color, horizontalalignment="center", size=hea_a_font,
                                     verticalalignment="center")


def _plot_hist(
        ax,
        his_data: list[numpy__ndarray] =None,
        his_stat: list[str] = None,
        his_bins: list[int] = None,
        his_binwidth: list[float] = None,
        his_binrange: list[tuple[float, float]] = None,
        his_cumulative: list[bool] = None,
        his_element: list[str] = None,
        his_fill: list[bool] = None,
        his_shrink: list[float] = None,
        his_kde: list[bool] = None,
        his_color: list[str] = None,
        **kwargs):
    if isinstance(his_data, list) and len(his_data) > 0:
        # control curve parameters
        if his_stat is None or (isinstance(his_stat, list) and len(his_stat) != len(his_data)):
            his_stat = ["probability"] * len(his_data)
        if his_bins is None or (isinstance(his_bins, list) and len(his_bins) != len(his_data)):
            his_bins = ["auto"] * len(his_data)
        if his_binwidth is None or (isinstance(his_binwidth, list) and len(his_binwidth) != len(his_data)):
            his_binwidth = [None] * len(his_data)
        if his_binrange is None or (isinstance(his_binrange, list) and len(his_binrange) != len(his_data)):
            his_binrange = [None] * len(his_data)
        if his_cumulative is None or (isinstance(his_cumulative, list) and len(his_cumulative) != len(his_data)):
            his_cumulative = [False] * len(his_data)
        if his_element is None or (isinstance(his_element, list) and len(his_element) != len(his_data)):
            his_element = ["step"] * len(his_data)
        if his_fill is None or (isinstance(his_fill, list) and len(his_fill) != len(his_data)):
            his_fill = [False] * len(his_data)
        if his_shrink is None or (isinstance(his_shrink, list) and len(his_shrink) != len(his_data)):
            his_shrink = [1] * len(his_data)
        if his_kde is None or (isinstance(his_kde, list) and len(his_kde) != len(his_data)):
            his_kde = [False] * len(his_data)
        if his_color is None or (isinstance(his_color, list) and len(his_color) != len(his_data)):
            his_color = ["k"] * len(his_data)
        # plot histograms
        for bi, br, bw, co, cu, da, e, f, k, sh, st in zip(
                his_bins, his_binrange, his_binwidth, his_color, his_cumulative, his_data, his_element, his_fill,
                his_kde, his_shrink, his_stat):
            seaborn.histplot(data=da, bins=bi, binrange=br, binwidth=bw, color=co, cumulative=cu, element=e, fill=f,
                             kde=k, shrink=sh, stat=st, ax=ax)
        # redo axes, the boxplot function causes troubles
        _ax_axis_x(ax, **kwargs)
        _ax_axis_y(ax, **kwargs)


def _plot_marker(
        ax,
        mar_c: list = None,
        mar_lc: list = None,
        mar_ls: list = None,
        mar_lw: list = None,
        mar_m: list = None,
        mar_s: list = None,
        mar_x: list = None,
        mar_y: list = None,
        mar_z: list = None,
        mar_colorbar: str = None,
        mar_colorbar_lab: list = None,
        mar_colorbar_nam: str = None,
        mar_colorbar_position: str = None,
        mar_colorbar_tic: list = None,
        **kwargs):
    """
    Plot markers

    Input:
    ------
    :param ax: matplotlib Axes object
    :param mar_c: list, optional
        List of facecolor for markers; e.g., mar_c = ['k', 'b']
    :param mar_colorbar: str, optional
        Name of a mar_colorbar_tics; e.g., mar_colorbar = 'cmo.rain'
    :param mar_colorbar_lab: list, optional
        Tick labels for colorbar; e.g., mar_colorbar_lab = ['1', '2']
    :param mar_colorbar_nam: str, optional
        Label for the colorbar; e.g., mar_colorbar_nam = 'temperature'
        Default is ''
    :param mar_colorbar_position: str, optional
        Position of the colorbar; e.g., mar_colorbar_position = 'bottom'
        Two colorbar positions are accepted: 'bottom', 'right'
        Default is 'bottom'
    :param mar_colorbar_tic: list, optional
        Ticks for colorbar; e.g., mar_colorbar_tic = [1, 2]
    :param mar_ls: list, optional
        List of linestyle for markers; e.g., mar_ls = ['-', ':']
    :param mar_lc: list, optional
        List of edgecolor for markers; e.g., mar_lc = ['k', 'b']
    :param mar_lw: list, optional
        List of linewidth for markers; e.g., mar_lw = [0.5, 0.5]
    :param mar_m: list, optional
        List of marker for markers; e.g., mar_m = ['D', 'o']
    :param mar_s: list, optional
        List of size for markers; e.g., mar_s = [50, 100]
    :param mar_x: list, optional
        List of x values for markers; e.g., mar_x = [ [1, 2, 3], [1, 2, 3] ]
    :param mar_y: list, optional
        List of y values for markers; e.g., mar_y = [ [10, 2, 10], [3, 6, 9] ]
    :param mar_z: list, optional
        List of zorder for markers (drawing order); e.g., mar_z = [2, 4]
    """
    # set default values if not given
    mar_colorbar_lab = none_to_default(mar_colorbar_lab, default_plot["lab"])
    mar_colorbar_nam = none_to_default(mar_colorbar_nam, default_plot["text"])
    mar_colorbar_position = none_to_default(mar_colorbar_position, default_plot["marker_position"])
    mar_colorbar_tic = none_to_default(mar_colorbar_tic, default_plot["tic"])
    # left axis
    if isinstance(mar_x, list) is True and isinstance(mar_y, list) is True and len(mar_x) == len(mar_y):
        # control marker parameters
        if mar_c is None or (isinstance(mar_c, list) is True and len(mar_c) != len(mar_x)):
            mar_c = [default_plot["color"]] * len(mar_x)
        if mar_lc is None or (isinstance(mar_lc, list) is True and len(mar_lc) != len(mar_x)):
            mar_lc = [default_plot["color"]] * len(mar_x)
        if mar_ls is None or (isinstance(mar_ls, list) is True and len(mar_ls) != len(mar_x)):
            mar_ls = [default_plot["linestyle"]] * len(mar_x)
        if mar_lw is None or (isinstance(mar_lw, list) is True and len(mar_lw) != len(mar_x)):
            mar_lw = [default_plot["linewidth"]] * len(mar_x)
        if mar_m is None or (isinstance(mar_m, list) is True and len(mar_m) != len(mar_x)):
            mar_m = [default_plot["marker"]] * len(mar_x)
        if mar_s is None or (isinstance(mar_s, list) is True and len(mar_s) != len(mar_x)):
            mar_s = [default_plot["size"]] * len(mar_x)
        if mar_z is None or (isinstance(mar_z, list) is True and len(mar_z) != len(mar_x)):
            mar_z = [default_plot["marker_zorder"]] * len(mar_x)
        # put parameters that are unique in one dictionary, the ones that are not in another
        ll1 = [mar_c, mar_lc, mar_ls, mar_lw, mar_m, mar_s, mar_z]
        ll2 = ["facecolors", "edgecolors", "linestyle", "linewidths", "marker", "s", "zorder"]
        dict_sca, dict_lis = {}, {}
        for k1, k2 in zip(ll1, ll2):
            ll3 = list(set(k1))
            if len(ll3) == 1:
                dict_sca[k2] = ll3[0]
            else:
                dict_lis[k2] = {"full": k1, "unique": ll3}
        scat = None
        if isinstance(mar_colorbar, str) is True and all([isinstance(k, str) is False for k in mar_c]):
            dif, mul = round(float(mar_colorbar_tic[1] - mar_colorbar_tic[0]), 5), 5
            delta = round(dif / mul, 5)
            lev = [round(k1 + k2 * delta, 5) for k1 in mar_colorbar_tic[:-1]
                   for k2 in range(mul)] + [mar_colorbar_tic[-1]]
            if mar_colorbar_nam == "correlation":
                cmap = pyplot.get_cmap("cmo.balance")
                cmap = cmap(numpy__linspace(0.15, 0.85, len(lev) - 1))
                cmap = ListedColormap(cmap)
            else:
                cmap = pyplot.get_cmap(mar_colorbar)
                cmap = cmap(numpy__linspace(0, 1, len(lev) - 1))
            cmap.set_bad("black")
            norm = BoundaryNorm(lev, ncolors=cmap.N, clip=True)
            dict_sca.update({"cmap": cmap, "norm": norm})
        else:
            mar_colorbar = None
        # plot markers
        if len(list(dict_lis.keys())) == 0:
            # all markers have the same parameters
            scat = ax.scatter(mar_x, mar_y, clip_on=False, **dict_sca)
        elif (len(mar_x) < 101 or len(list(dict_lis.keys())) > 4) and mar_colorbar is None:
            # few markers to plot --> plot them one by one
            # or
            # too few unique parameters --> plot markers one by one
            for c, lc, ls, lw, m, s, x, y, z in zip(mar_c, mar_lc, mar_ls, mar_lw, mar_m, mar_s, mar_x, mar_y, mar_z):
                scat = ax.scatter(x, y, marker=m, s=s, facecolors=c, edgecolors=lc, linestyle=ls, linewidths=lw,
                                  clip_on=False, zorder=z)
        else:
            # there are several non-unique parameters
            list_keys = list(dict_lis.keys())
            n1 = list_keys[0]
            f1, u1 = dict_lis[n1]["full"], dict_lis[n1]["unique"]
            if len(list_keys) == 1:
                for k1 in u1:
                    dict_t = {**dict_sca, **{"clip_on": False, n1: k1}}
                    list_x, list_y = [], []
                    for i1, xx, yy in zip(f1, mar_x, mar_y):
                        if i1 == k1:
                            list_x.append(xx)
                            list_y.append(yy)
                    scat = ax.scatter(list_x, list_y, **dict_t)
            elif len(list_keys) == 2:
                n2 = list_keys[1]
                f2, u2 = dict_lis[n2]["full"], dict_lis[n2]["unique"]
                for k1 in u1:
                    for k2 in u2:
                        dict_t = {**dict_sca, **{"clip_on": False, n1: k1, n2: k2}}
                        list_x, list_y = [], []
                        for i1, i2, xx, yy in zip(f1, f2, mar_x, mar_y):
                            if i1 == k1 and i2 == k2:
                                list_x.append(xx)
                                list_y.append(yy)
                        scat = ax.scatter(list_x, list_y, **dict_t)
            elif len(list_keys) == 3:
                n2 = list_keys[1]
                f2, u2 = dict_lis[n2]["full"], dict_lis[n2]["unique"]
                n3 = list_keys[2]
                f3, u3 = dict_lis[n3]["full"], dict_lis[n3]["unique"]
                for k1 in u1:
                    for k2 in u2:
                        for k3 in u3:
                            dict_t = {**dict_sca, **{"clip_on": False, n1: k1, n2: k2, n3: k3}}
                            list_x, list_y = [], []
                            for i1, i2, i3, xx, yy in zip(f1, f2, f3, mar_x, mar_y):
                                if i1 == k1 and i2 == k2 and i3 == k3:
                                    list_x.append(xx)
                                    list_y.append(yy)
                            scat = ax.scatter(list_x, list_y, **dict_t)
            elif len(list_keys) == 4:
                n2 = list_keys[1]
                f2, u2 = dict_lis[n2]["full"], dict_lis[n2]["unique"]
                n3 = list_keys[2]
                f3, u3 = dict_lis[n3]["full"], dict_lis[n3]["unique"]
                n4 = list_keys[3]
                f4, u4 = dict_lis[n4]["full"], dict_lis[n4]["unique"]
                for k1 in u1:
                    for k2 in u2:
                        for k3 in u3:
                            for k4 in u4:
                                dict_t = {**dict_sca, **{"clip_on": False, n1: k1, n2: k2, n3: k3, n4: k4}}
                                list_x, list_y = [], []
                                for i1, i2, i3, i4, xx, yy in zip(f1, f2, f3, f4, mar_x, mar_y):
                                    if i1 == k1 and i2 == k2 and i3 == k3 and i4 == k4:
                                        list_x.append(xx)
                                        list_y.append(yy)
                                scat = ax.scatter(list_x, list_y, **dict_t)
            else:
                print_fail(inspect__stack(), "too many non unique parameters for scatter")
        # colorbar
        fontsize = kwargs["fontsize"] if "fontsize" in list(kwargs.keys()) else default_plot["fontsize"]
        if mar_colorbar is not None and mar_colorbar_position == "bottom":
            x1, x2 = ax.get_position().x0, ax.get_position().x1
            y1, y2 = ax.get_position().y0, ax.get_position().y1
            cax = pyplot.axes((x1, y1 - (y2 - y1) / 3.5, x2 - x1, (y2 - y1) / 15))
            cbar = pyplot.colorbar(
                scat, cax=cax, orientation="horizontal", ticks=mar_colorbar_tic, label=mar_colorbar_lab, pad=0.3)
            cbar.ax.tick_params(labelsize=fontsize)
            cbar.set_label(mar_colorbar_nam, fontsize=fontsize, labelpad=2)
        elif mar_colorbar is not None and mar_colorbar_position == "right":
            x1, x2 = ax.get_position().x0, ax.get_position().x1
            y1, y2 = ax.get_position().y0, ax.get_position().y1
            cax = pyplot.axes((x2 + (x2 - x1) / 20, y1, (x2 - x1) / 30, y2 - y1))
            cbar = pyplot.colorbar(
                scat, cax=cax, orientation="vertical", ticks=mar_colorbar_tic, label=mar_colorbar_lab, pad=0.35)
            cbar.ax.tick_params(labelsize=fontsize)
            cbar.set_label(mar_colorbar_nam, fontsize=fontsize, rotation=90)


def _plot_marker2(
        ax2,
        mar2_c: list = None,
        mar2_lc: list = None,
        mar2_ls: list = None,
        mar2_lw: list = None,
        mar2_m: list = None,
        mar2_s: list = None,
        mar2_x: list = None,
        mar2_y: list = None,
        mar2_z: list = None,
        mar2_colorbar: str = None,
        mar2_colorbar_lab: list = None,
        mar2_colorbar_nam: str = None,
        mar2_colorbar_position: str = None,
        mar2_colorbar_tic: list = None,
        **kwargs):
    tmp_kwargs = {}
    l1 = ["mar_c", "mar_lc", "mar_ls", "mar_lw", "mar_m", "mar_s", "mar_x", "mar_y", "mar_z"]
    l2 = [mar2_c, mar2_lc, mar2_ls, mar2_lw, mar2_m, mar2_s, mar2_x, mar2_y, mar2_z]
    for k1, k2 in zip(l1, l2):
        if k2 is not None:
            tmp_kwargs[k1] = k2
    _plot_marker(ax2, **tmp_kwargs)


def _plot_shading(
        ax,
        sha_a=None,
        sha_c=None,
        sha_x=None,
        sha_y1=None,
        sha_y2=None,
        sha_z=None,
        **kwargs):
    """
    Plot shading between y1 and y2

    Input:
    ------
    :param ax: matplotlib Axes object
    :param sha_a: list[float | int], optional
        List of alpha for shadings; e.g., sha_c = [0.5, 1]
    :param sha_c: list, optional
        List of color for shadings; e.g., sha_c = ['k', 'b']
    :param sha_x: list, optional
        List of x values for shadings; e.g., sha_x = [ [1, 2, 3], [1, 2, 3] ]
    :param sha_y1: list, optional
        List of lowest y values for shadings; e.g., sha_y1 = [ [1, 2, 3], [3, 4, 5] ]
    :param sha_y2: list, optional
        List of highest y values for shadings; e.g., sha_y2 = [ [2, 3, 4], [3, 6, 9] ]
    :param sha_z: list, optional
        List of zorder for shadings (drawing order); e.g., sha_z = [2, 4]
    """
    # left axis
    if isinstance(sha_x, list) is True and isinstance(sha_y1, list) is True and isinstance(sha_y2, list) is True \
            and len(sha_x) == len(sha_y1) and len(sha_x) == len(sha_y2):
        # control shading parameters
        if sha_a is None or (isinstance(sha_a, list) is True and len(sha_a) != len(sha_x)):
            sha_a = [default_plot["shading_alpha"]] * len(sha_x)
        if sha_c is None or (isinstance(sha_c, list) is True and len(sha_c) != len(sha_x)):
            sha_c = [default_plot["shading_color"]] * len(sha_x)
        if sha_z is None or (isinstance(sha_z, list) is True and len(sha_z) != len(sha_x)):
            sha_z = [default_plot["shading_zorder"]] * len(sha_x)
        # plot shadings
        for a, c, x, y1, y2, z in zip(sha_a, sha_c, sha_x, sha_y1, sha_y2, sha_z):
            ax.fill_between(x, y1, y2=y2, alpha=a, color=c, linewidth=0.0, zorder=z)


def _plot_violin(
        ax,
        vio_a: list[float] = None,
        vio_b: list = None,
        vio_bar_ec: list = None,
        vio_bar_ls: list = None,
        vio_bar_lw: list[float] = None,
        vio_bod_ec: list = None,
        vio_bod_fc: list = None,
        vio_bod_ls: list = None,
        vio_bod_lw: list[float] = None,
        vio_c: list = None,
        vio_flie: list[bool] = None,
        vio_ls: list = None,
        vio_lw: list[float] = None,
        vio_max_ec: list = None,
        vio_max_ls: list = None,
        vio_max_lw: list[float] = None,
        vio_mea: list[bool] = None,
        vio_mea_ec: list = None,
        vio_mea_ls: list = None,
        vio_mea_lw: list[float] = None,
        vio_med: list[bool] = None,
        vio_med_ec: list = None,
        vio_med_ls: list = None,
        vio_med_lw: list[float] = None,
        vio_min_ec: list = None,
        vio_min_ls: list = None,
        vio_min_lw: list[float] = None,
        vio_o: list = None,
        vio_p: list[int] = None,
        vio_q: list[list[float]] = None,
        vio_qua_ec: list = None,
        vio_qua_ls: list = None,
        vio_qua_lw: list[float] = None,
        vio_s: list = None,
        vio_w: list[float] = None,
        vio_x: list[float] = None,
        vio_y: list[list[float]] = None,
        vio_z: list = None,
        **kwargs):
    """
    Plot violinplot

    May switch to seaborn
    https://seaborn.pydata.org/generated/seaborn.violinplot.html

    Input:
    ------
    :param ax: matplotlib Axes object
    :param vio_a: list, optional
        List of alpha for violinplot; e.g., vio_a = [0.5, 0.5]
    :param vio_b: list[{'scott', 'silverman'}], optional
        List of {'scott', 'silverman'} method used to calculate the estimator bandwidth for violinplot;
        e.g., vio_c = ['scott', 'scott']
    :param vio_bar_ec: list, optional
        List of edgecolor (color of the cbars of the violin) for violinplot; e.g., vio_bar_ec = ['k', 'b']
    :param vio_bar_ls: list, optional
        List of linestyle (linestyle of the cbars of the violin) for violinplot; e.g., vio_bar_ls = ['-', ':']
    :param vio_bar_lw: list[float], optional
        List of linewidth (linewidth of the cbars of the violin) for violinplot; e.g., vio_bar_lw = [1, 2]
    :param vio_bod_ec: list, optional
        List of edgecolor (color of the edge of the violin) for violinplot; e.g., vio_bod_ec = ['k', 'b']
        If not given, no edge
    :param vio_bod_fc: list, optional
        List of facecolor (color of the body of the violin) for violinplot; e.g., vio_bod_fc = ['k', 'b']
        If not given, vio_c is used
    :param vio_bod_ls: list, optional
        List of linestyle (linestyle of the edge of the violin) for violinplot; e.g., vio_bod_ls = ['-', ':']
        '-' if not given but vio_bod_ec is given, else no edge
    :param vio_bod_lw: list[float], optional
        List of linewidth (linewidth of the edge of the violin) for violinplot; e.g., vio_bod_lw = [1, 2]
        1 if not given but vio_bod_ec is given, else no edge
    :param vio_c: list, optional
        List of color for violinplot; e.g., vio_c = ['k', 'b']
    :param vio_flie: list[bool], optional
        List of boolean to show extrema with a line for violinplot, if True, show the extrema;
        e.g., vio_flie = [True, False]
    :param vio_ls: list, optional
        List of linestyle (used for cbars, cmaxes, cmeans, cmedians, cmins, cquantiles if not provided) for violinplot;
        e.g., vio_ls = ['-', ':']
    :param vio_lw: list[float], optional
        List of linewidth (used for cbars, cmaxes, cmeans, cmedians, cmins, cquantiles if not provided) for violinplot;
        e.g., vio_lw = [1, 2]
    :param vio_max_ec: list, optional
        List of edgecolor (color of the cmaxes of the violin) for violinplot; e.g., vio_max_ec = ['k', 'b']
    :param vio_max_ls: list, optional
        List of linestyle (linestyle of the cmaxes of the violin) for violinplot; e.g., vio_max_ls = ['-', ':']
    :param vio_max_lw: list[float], optional
        List of linewidth (linewidth of the cmaxes of the violin) for violinplot; e.g., vio_max_lw = [1, 2]
    :param vio_mea: list[bool], optional
        List of boolean to show mean with a line for violinplot, if True, show the mean; e.g., vio_mea = [True, False]
    :param vio_mea_ec: list, optional
        List of edgecolor (color of the cmeans of the violin) for violinplot; e.g., vio_mea_ec = ['k', 'b']
    :param vio_mea_ls: list, optional
        List of linestyle (linestyle of the cmeans of the violin) for violinplot; e.g., vio_mea_ls = ['-', ':']
    :param vio_mea_lw: list[float], optional
        List of linewidth (linewidth of the cmeans of the violin) for violinplot; e.g., vio_mea_lw = [1, 2]
    :param vio_med: list[bool], optional
        List of boolean to show meadian with a line for violinplot, if True, show the meadian;
        e.g., vio_mea = [True, False]
    :param vio_med_ec: list, optional
        List of edgecolor (color of the cmedians of the violin) for violinplot; e.g., vio_med_ec = ['k', 'b']
    :param vio_med_ls: list, optional
        List of linestyle (linestyle of the cmedians of the violin) for violinplot; e.g., vio_med_ls = ['-', ':']
    :param vio_med_lw: list[float], optional
        List of linewidth (linewidth of the cmedians of the violin) for violinplot; e.g., vio_med_lw = [1, 2]
    :param vio_min_ec: list, optional
        List of edgecolor (color of the cmins of the violin) for violinplot; e.g., vio_min_ec = ['k', 'b']
    :param vio_min_ls: list, optional
        List of linestyle (linestyle of the cmins of the violin) for violinplot; e.g., vio_min_ls = ['-', ':']
    :param vio_min_lw: list[float], optional
        List of linewidth (linewidth of the cmins of the violin) for violinplot; e.g., vio_min_lw = [1, 2]
    :param vio_o: list[{'vertical', 'horizontal'}], optional
        List of {'vertical', 'horizontal'} orientations for violinplot; e.g., vio_o = ['vertical', 'vertical']
    :param vio_p: list[int], optional
        List of number of points to evaluate each of the gaussian kernel density estimations at for violinplot;
        e.g., vio_p = [100, 100]
    :param vio_q: list[list[float, float]], optional
        List of quantiles [0, 1] for each violin, which stands for the quantiles that will be rendered for that violin.
        number of points to evaluate each of the gaussian kernel density estimations at for violinplot;
        e.g., vio_p = [100, 100]
    :param vio_qua_ec: list, optional
        List of edgecolor (color of the cquantiles of the violin) for violinplot; e.g., vio_qua_ec = ['k', 'b']
    :param vio_qua_ls: list, optional
        List of linestyle (linestyle of the cquantiles of the violin) for violinplot; e.g., vio_qua_ls = ['-', ':']
    :param vio_qua_lw: list[float], optional
        List of linewidth (linewidth of the cquantiles of the violin) for violinplot; e.g., vio_qua_lw = [1, 2]
    :param vio_s: list[{'both', 'low', 'high'}], optional
        List of {'both', 'low', 'high'} side for violinplot; e.g., vio_o = ['both', 'both']
    :param vio_w: list[float], optional
        List of width for violinplot; e.g., vio_w = [0.5, 0.5]
    :param vio_x: list[float], optional
        List of x values (positions) for violinplot; e.g., vio_x = [1, 2]
    :param vio_y: list[list[float]], optional
        List of y values (data for violin) for violinplot; e.g., vio_y = [ [1, 2, 3], [1, 2, 3] ]
    :param vio_z: list, optional
        List of zorder for violinplot (drawing order); e.g., vio_z = [2, 4]
    """
    if isinstance(vio_x, list) is True and isinstance(vio_y, list) is True and len(vio_x) == len(vio_y):
        # control violinplot parameters
        if not is_list_to_plot(vio_a, vio_x):
            vio_a = [default_plot["violinplot_alpha"]] * len(vio_x)
        if not is_list_to_plot(vio_b, vio_x):
            vio_b = [default_plot["violinplot_bw_method"]] * len(vio_x)
        if not is_list_to_plot(vio_c, vio_x):
            vio_c = [default_plot["color"]] * len(vio_x)
        if not is_list_to_plot(vio_bod_fc, vio_x):
            vio_bod_fc = [k for k in vio_c]
        if not is_list_to_plot(vio_flie, vio_x):
            vio_flie = [default_plot["violinplot_flier"]] * len(vio_x)
        if not is_list_to_plot(vio_ls, vio_x):
            vio_ls = [default_plot["linestyle"]] * len(vio_x)
        if not is_list_to_plot(vio_lw, vio_x):
            vio_lw = [default_plot["linewidth"]] * len(vio_x)
        if not is_list_to_plot(vio_mea, vio_x):
            vio_mea = [default_plot["violinplot_mean"]] * len(vio_x)
        if not is_list_to_plot(vio_med, vio_x):
            vio_med = [default_plot["violinplot_median"]] * len(vio_x)
        if not is_list_to_plot(vio_o, vio_x):
            vio_o = [default_plot["violinplot_orientation"]] * len(vio_x)
        if not is_list_to_plot(vio_p, vio_x):
            vio_p = [default_plot["violinplot_points"]] * len(vio_x)
        if not is_list_to_plot(vio_q, vio_x):
            vio_q = [default_plot["violinplot_quantiles"]] * len(vio_x)
        if not is_list_to_plot(vio_s, vio_x):
            vio_s = [default_plot["violinplot_side"]] * len(vio_x)
        if not is_list_to_plot(vio_w, vio_x):
            dx = 1e20
            for k1 in range(len(vio_x)):
                for k2 in range(k1 + 1, len(vio_x)):
                    if abs(vio_x[k1] - vio_x[k2]) < dx:
                        dx = abs(vio_x[k1] - vio_x[k2])
            if dx == 1e20:
                dx = 1
            vio_w = [dx * 0.7] * len(vio_x)
        if not is_list_to_plot(vio_z, vio_x):
            vio_z = [default_plot["violinplot_zorder"]] * len(vio_x)
        # -- combined param
        # cbars
        if not is_list_to_plot(vio_bar_ec, vio_x):
            vio_bar_ec = [None] * len(vio_x)
        if not is_list_to_plot(vio_bar_ls, vio_x):
            vio_bar_ls = [None] * len(vio_x)
        if not is_list_to_plot(vio_bar_lw, vio_x):
            vio_bar_lw = [k2 if k1 is True else 0 for k1, k2 in zip(vio_flie, vio_lw)]
        for lk, (ec, ls, lw, ) in enumerate(zip(vio_bar_ec, vio_bar_ls, vio_bar_lw)):
            if ec is None and (ls is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_bar_ec[lk] = vio_c[lk]
            if ls is None and (ec is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_bar_ls[lk] = vio_ls[lk]
            if (isinstance(lw, (float, int)) is False or lw == 0) and (ec is not None or ls is not None):
                vio_bar_lw[lk] = vio_lw[lk]
        # bodies
        if not is_list_to_plot(vio_bod_ec, vio_x):
            vio_bod_ec = [None] * len(vio_x)
        if not is_list_to_plot(vio_bod_ls, vio_x):
            vio_bod_ls = [None] * len(vio_x)
        if not is_list_to_plot(vio_bod_lw, vio_x):
            vio_bod_lw = [0] * len(vio_x)
        for lk, (ec, ls, lw) in enumerate(zip(vio_bod_ec, vio_bod_ls, vio_bod_lw)):
            if ec is None and (ls is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_bod_ec[lk] = vio_c[lk]
            if ls is None and (ec is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_bod_ls[lk] = vio_ls[lk]
            if (isinstance(lw, (float, int)) is False or lw == 0) and (ec is not None or ls is not None):
                vio_bod_lw[lk] = vio_lw[lk]
        # cmaxes
        if not is_list_to_plot(vio_max_ec, vio_x):
            vio_max_ec = [None] * len(vio_x)
        if not is_list_to_plot(vio_max_ls, vio_x):
            vio_max_ls = [None] * len(vio_x)
        if not is_list_to_plot(vio_max_lw, vio_x):
            vio_max_lw = [k2 if k1 is True else 0 for k1, k2 in zip(vio_flie, vio_lw)]
        for lk, (ec, ls, lw) in enumerate(zip(vio_max_ec, vio_max_ls, vio_max_lw)):
            if ec is None and (ls is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_max_ec[lk] = vio_c[lk]
            if ls is None and (ec is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_max_ls[lk] = vio_ls[lk]
            if (isinstance(lw, (float, int)) is False or lw == 0) and (ec is not None or ls is not None):
                vio_max_lw[lk] = vio_lw[lk]
        # cmeans
        if not is_list_to_plot(vio_mea_ec, vio_x):
            vio_mea_ec = [None] * len(vio_x)
        if not is_list_to_plot(vio_mea_ls, vio_x):
            vio_mea_ls = [None] * len(vio_x)
        if not is_list_to_plot(vio_mea_lw, vio_x):
            vio_mea_lw = [k2 if k1 is True else 0 for k1, k2 in zip(vio_mea, vio_lw)]
        for lk, (ec, ls, lw) in enumerate(zip(vio_mea_ec, vio_mea_ls, vio_mea_lw)):
            if ec is None and (ls is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_mea_ec[lk] = vio_c[lk]
            if ls is None and (ec is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_mea_ls[lk] = vio_ls[lk]
            if (isinstance(lw, (float, int)) is False or lw == 0) and (ec is not None or ls is not None):
                vio_mea_lw[lk] = vio_lw[lk]
        # cmedians
        if not is_list_to_plot(vio_med_ec, vio_x):
            vio_med_ec = [None] * len(vio_x)
        if not is_list_to_plot(vio_med_ls, vio_x):
            vio_med_ls = [None] * len(vio_x)
        if not is_list_to_plot(vio_med_lw, vio_x):
            vio_med_lw = [k2 if k1 is True else 0 for k1, k2 in zip(vio_med, vio_lw)]
        for lk, (ec, ls, lw) in enumerate(zip(vio_med_ec, vio_med_ls, vio_med_lw)):
            if ec is None and (ls is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_med_ec[lk] = vio_c[lk]
            if ls is None and (ec is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_med_ls[lk] = vio_ls[lk]
            if (isinstance(lw, (float, int)) is False or lw == 0) and (ec is not None or ls is not None):
                vio_med_lw[lk] = vio_lw[lk]
        # cmins
        if not is_list_to_plot(vio_min_ec, vio_x):
            vio_min_ec = [None] * len(vio_x)
        if not is_list_to_plot(vio_min_ls, vio_x):
            vio_min_ls = [None] * len(vio_x)
        if not is_list_to_plot(vio_min_lw, vio_x):
            vio_min_lw = [k2 if k1 is True else 0 for k1, k2 in zip(vio_flie, vio_lw)]
        for lk, (ec, ls, lw) in enumerate(zip(vio_min_ec, vio_min_ls, vio_min_lw)):
            if ec is None and (ls is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_min_ec[lk] = vio_c[lk]
            if ls is None and (ec is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_min_ls[lk] = vio_ls[lk]
            if (isinstance(lw, (float, int)) is False or lw == 0) and (ec is not None or ls is not None):
                vio_min_lw[lk] = vio_lw[lk]
        # cquantiles
        if not is_list_to_plot(vio_qua_ec, vio_x):
            vio_qua_ec = [None] * len(vio_x)
        if not is_list_to_plot(vio_qua_ls, vio_x):
            vio_qua_ls = [None] * len(vio_x)
        if not is_list_to_plot(vio_qua_lw, vio_x):
            vio_qua_lw = [k2 if k1 is True else 0 for k1, k2 in zip(vio_q, vio_lw)]
        for lk, (ec, ls, lw) in enumerate(zip(vio_qua_ec, vio_qua_ls, vio_qua_lw)):
            if ec is None and (ls is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_qua_ec[lk] = vio_c[lk]
            if ls is None and (ec is not None or (isinstance(lw, (float, int)) is True and lw > 0)):
                vio_qua_ls[lk] = vio_ls[lk]
            if (isinstance(lw, (float, int)) is False or lw == 0) and (ec is not None or ls is not None):
                vio_qua_lw[lk] = vio_lw[lk]
        # plot violinplot
        for lk, (aa, bw, ba_c, ba_s, ba_w, bo_ce, bo_cf, bo_s, bo_w, ff, ma_c, ma_s, ma_w, me, me_c, me_s, me_w, md,
                 md_c, md_s, md_w, mi_c, mi_s, mi_w, oo, pp, qu, qu_c, qu_s, qu_w, ss, ww, xx, yy, zz) in enumerate(
                zip(vio_a, vio_b, vio_bar_ec, vio_bar_ls, vio_bar_lw, vio_bod_ec, vio_bod_fc, vio_bod_ls, vio_bod_lw,
                    vio_flie, vio_max_ec, vio_max_ls, vio_max_lw, vio_mea, vio_mea_ec, vio_mea_ls, vio_mea_lw, vio_med,
                    vio_med_ec, vio_med_ls, vio_med_lw, vio_min_ec, vio_min_ls, vio_min_lw, vio_o, vio_p, vio_q,
                    vio_qua_ec, vio_qua_ls, vio_qua_lw, vio_s, vio_w, vio_x, vio_y, vio_z)):
            # cc
            parts = ax.violinplot(yy, bw_method=bw, orientation=oo, points=pp, positions=[xx], quantiles=qu,
                                  showextrema=ff, showmeans=me, showmedians=md, side=ss, widths=ww)
            parts["bodies"][0].set_alpha(aa)
            parts["bodies"][0].set_facecolor(bo_cf)
            if isinstance(bo_w, (float, int)) is True and bo_w > 0:
                parts["bodies"][0].set_edgecolor(bo_ce)
                parts["bodies"][0].set_linestyle(bo_s)
                parts["bodies"][0].set_linewidth(bo_w)
            parts["bodies"][0].set_zorder(zz)
            l1 = ["cbars", "cmaxes", "cmeans", "cmedians", "cmins", "cquantiles"]
            l2 = [ba_c, ma_c, me_c, md_c, mi_c, qu_c]
            l3 = [ba_s, ma_s, me_s, md_s, mi_s, qu_s]
            l4 = [ba_w, ma_w, me_w, md_w, mi_w, qu_w]
            for part_name, ec, ls, lw in zip(l1, l2, l3, l4):
                if part_name in list(parts.keys()) and (isinstance(lw, (float, int)) is True and lw > 0):
                    parts[part_name].set_color(ec)
                    parts[part_name].set_linestyle(ls)
                    parts[part_name].set_linewidth(lw)
                    parts[part_name].set_zorder(zz)
        # redo axes, the boxplot function causes troubles
        _ax_axis_x(ax, **kwargs)
        _ax_axis_y(ax, **kwargs)


def plot_main(ax, **kwargs):
    """
    Plot given boxplots, curves, histograms, markers and shadings

    Inputs:
    -------
    :param ax: matplotlib Axes object

    """
    # x-axis
    _ax_axis_x(ax, **kwargs)
    # y-axis
    _ax_axis_y(ax, **kwargs)
    # shading
    _plot_shading(ax, **kwargs)
    # boxplot
    _plot_boxplot(ax, **kwargs)
    # violinplot
    _plot_violin(ax, **kwargs)
    # curve
    _plot_curve(ax, **kwargs)
    # if "cur2_x" in list(kwargs.keys()) and "cur2_y" in list(kwargs.keys()):
    #     ax2 = ax.twinx()
    #     if "y2_lim" in list(kwargs.keys()) and "y2_tic" in list(kwargs.keys()):
    #         _ax_axis_y2(ax2, **kwargs)
    #     _plot_curve2(ax2, **kwargs)
    # density
    _plot_density(ax, **kwargs)
    # heatmap
    _plot_heatmap(ax, **kwargs)
    # histogram
    _plot_hist(ax, **kwargs)
    # marker
    # if "mar2_x" in list(kwargs.keys()) and "mar2_y" in list(kwargs.keys()):
    #     ax2 = ax.twinx()
    #     if "y2_lim" in list(kwargs.keys()) and "y2_tic" in list(kwargs.keys()):
    #         _ax_axis_y2(ax2, **kwargs)
    #     _plot_marker2(ax2, **kwargs)
    _plot_marker(ax, **kwargs)
    # text
    _ax_text(ax, **kwargs)
    # title
    _ax_titles(ax, **kwargs)
    # annotate
    _plot_annotate(ax, **kwargs)
    # legend
    _ax_legend(ax, **kwargs)


def plot_map(
        ax,
        fontsize: float = None,
        legend_position: str = "bottom",
        projection: object = None,
        region: list = None,
        map_c=None,
        map_c_cbar: str = None,
        map_c_lab: list = None,
        map_c_land: str = None,
        map_c_nam: str = "",
        map_c_tic: list = None,
        map_h=None,
        map_l1=None,
        map_l1_tic: list = None,
        map_l2=None,
        map_l2_tic: list = None,
        kwargs_h: dict = None,
        kwargs_l1: dict = None,
        kwargs_l2: dict = None,
        x_lab: list[str] = None,
        x_lim: list[Union[float, int]] = None,
        x_tic: list[Union[float, int]] = None,
        y_lab: list[str] = None,
        y_lim: list[Union[float, int]] = None,
        y_tic: list[Union[float, int]] = None,
        **kwargs):
    """
    Plot map

    :param ax: matplotlib Axes object
    :param fontsize: float, optional
        Size of the font; e.g., fontsize = 15.
    :param legend_position: str, optional
        Position of the legend; e.g., fig_legend_position = 'bottom'
        Two legend positions are accepted: 'bottom', 'right'
        Default is 'bottom'
    :param projection: object, optional
        Output from a cartopy projection (see https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html);
        e.g., projection = cartopy.crs.PlateCarree()
    :param region: list, optional
        Names of region to draw on the map, must be defined in default_region; e.g., region = ['n30e']
        Default is None (no region)
    :param map_c: array_like
        Array for map shading
    :param map_c_cbar: str, optional
        Name of a colorbar; e.g., map_c_bar = 'cmo.rain'
    :param map_c_lab: list, optional
        Tick labels for colorbar; e.g., map_c_lab = ['1', '2']
    :param map_c_land: str, optional
        Name of a color for the land; e.g., map_c_land = 'grey'
    :param map_c_nam: str, optional
        Label for the colorbar; e.g., map_c_nam = "sst"
        Default is ''
    :param map_c_tic: list, optional
        Ticks for colorbar; e.g., s_lab = [1, 2]
    :param map_h: array_like, optional
        Array for hatching (details must be passed in kwargs_h)
    :param map_l1: array_like, optional
        Array for the first set of contours (details must be passed in kwargs_l1)
    :param map_l1_tic: list, optional
        Tick level for map_l1 contours; e.g., map_l1_tic = [1, 2]
    :param map_l2: array_like, optional
        Array for the second set of contours (details must be passed in kwargs_l2)
    :param map_l2_tic: list, optional
        Tick level for map_l2 contours; e.g., map_l2_tic = [1, 2]
    :param kwargs_h: dict, optional
        Key arguments for matplotlib contourf to define map_h
        https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.contourf.html
    :param kwargs_l1: dict, optional
        Key arguments for matplotlib contour to define map_l1
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contour.html
    :param kwargs_l2: dict, optional
        Key arguments for matplotlib contour to define map_l2
    :param kwargs: dict, optional
        Extra keywords (for title), see _ax_titles
    """
    # set default values if not given
    fontsize = none_to_default(fontsize, default_plot["fontsize"])
    kwargs_h = none_to_default(kwargs_h, {})
    kwargs_l1 = none_to_default(kwargs_l1, {})
    kwargs_l2 = none_to_default(kwargs_l2, {})
    projection = none_to_default(projection, default_map["projection"])
    map_c_cbar = none_to_default(map_c_cbar, default_map["colorbar"])
    map_c_lab = none_to_default(map_c_lab, default_plot["lab"])
    map_c_tic = none_to_default(map_c_tic, default_plot["tic"])
    map_l1_tic = none_to_default(map_l1_tic, default_plot["tic"])
    map_l2_tic = none_to_default(map_l2_tic, default_plot["tic"])
    x_lab = none_to_default(x_lab, default_plot["lab"])
    x_tic = none_to_default(x_tic, default_plot["tic"])
    x_lim = none_to_default(x_lim, [min(x_tic), max(x_tic)])
    y_lab = none_to_default(y_lab, default_plot["lab"])
    y_tic = none_to_default(y_tic, default_plot["tic"])
    y_lim = none_to_default(y_lim, [min(y_tic), max(y_tic)])
    if map_c is not None:
        # select region
        # ax.coastlines()
        ax.set_extent(x_lim + y_lim, crs=projection)
        # draw coastlines
        # ax.coastlines()
        # fill continents
        if isinstance(map_c_land, str):
            ax.add_feature(cfeature.NaturalEarthFeature(
                "physical", "land", "110m", edgecolor="face", facecolor=map_c_land))
            # ax.add_feature(cfeature.NaturalEarthFeature(
            #     "physical", "land", "10m", edgecolor="face", facecolor=map_c_land))
            # ax.add_feature(cfeature.COASTLINE)
            # ax.coastlines()
        else:
            ax.add_feature(cfeature.COASTLINE)
        # x-y ticks and labels
        ax.set_xticks(x_tic, crs=projection)
        ax.set_yticks(y_tic, crs=projection)
        ax.set_xticklabels(x_lab)
        ax.set_yticklabels(y_lab)
        ax.tick_params(axis="both", direction="out", which="both", labelsize=fontsize, bottom=True, top=False,
                       left=True, right=True, labelbottom=True, labeltop=False, labelleft=True, labelright=False)
        # colorbar levels
        # dif, mul = round(float(map_c_tic[1] - map_c_tic[0]), 5), 20
        dif, mul = round(float(map_c_tic[1] - map_c_tic[0]), 5), 10
        delta = round(dif / mul, 5)
        lev = [round(k1 + k2 * delta, 5) for k1 in map_c_tic[:-1] for k2 in range(mul)] + [map_c_tic[-1]]
        # shading
        latitudes = map_c["latitude"].values if "latitude" in list(map_c.dims) else map_c["lat"].values
        longitudes = map_c["longitude"].values if "longitude" in list(map_c.dims) else map_c["lon"].values
        xx, yy = numpy__meshgrid(longitudes, latitudes)
        lc = ax.contourf(xx, yy, map_c, levels=lev, extend="both", transform=projection, cmap=map_c_cbar)
        # contours
        if map_l1 is not None:
            # contours
            latitudes = map_l1["latitude"].values if "latitude" in list(map_l1.dims) else map_l1["lat"].values
            longitudes = map_l1["longitude"].values if "longitude" in list(map_l1.dims) else map_l1["lon"].values
            xx, yy = numpy__meshgrid(longitudes, latitudes)
            # ll = ax.contour(xx, yy, map_l1, map_l1_tic, transform=ccrs.PlateCarree(), **kwargs_l1)
            # ax.clabel(ll, fontsize=fontsize, inline=True)
            z_smooth = gaussian_filter(map_l1, sigma=0.2)
            ll = ax.contour(xx, yy, z_smooth, map_l1_tic, transform=ccrs.PlateCarree(), **kwargs_l1)
            ax.clabel(ll, fontsize=fontsize, inline=True)
        if map_l2 is not None:
            # contours
            latitudes = map_l2["latitude"].values if "latitude" in list(map_l2.dims) else map_l2["lat"].values
            longitudes = map_l2["longitude"].values if "longitude" in list(map_l2.dims) else map_l2["lon"].values
            xx, yy = numpy__meshgrid(longitudes, latitudes)
            ll = ax.contour(xx, yy, map_l2, map_l2_tic, transform=ccrs.PlateCarree(), **kwargs_l2)
            ax.clabel(ll, fontsize=fontsize, inline=True)
        # hatches
        if map_h is not None:
            latitudes = map_h["latitude"].values if "latitude" in list(map_h.dims) else map_h["lat"].values
            longitudes = map_h["longitude"].values if "longitude" in list(map_h.dims) else map_h["lon"].values
            xx, yy = numpy__meshgrid(longitudes, latitudes)
            ax.contourf(xx, yy, map_h, transform=ccrs.PlateCarree(), colors="none", **kwargs_h)
        # title
        _ax_titles(ax, **kwargs)
        _ax_text(ax, **kwargs)
        _plot_annotate(ax, **kwargs)
        # region delimitation
        if isinstance(region, list) or isinstance(region, str):
            list_reg = copy__deepcopy(region) if isinstance(region, list) else [region]
            for k in list_reg:
                if len(default_region[k]["latitude"]) == 2 and len(default_region[k]["longitude"]) == 2:
                    lat_corners = numpy__array(default_region[k]["latitude"] +
                                               list(reversed(default_region[k]["latitude"])))
                    lon_corners = numpy__array([default_region[k]["longitude"][0]] * 2 +
                                               [default_region[k]["longitude"][1]] * 2)
                else:
                    lat_corners = numpy__array(default_region[k]["latitude"])
                    lon_corners = numpy__array(default_region[k]["longitude"])
                poly_corners = numpy__zeros((len(lat_corners), 2))
                poly_corners[:, 0] = lon_corners
                poly_corners[:, 1] = lat_corners
                pol = Polygon(poly_corners, alpha=default_region[k]["alpha"], linestyle=default_region[k]["linestyle"],
                              linewidth=default_region[k]["linewidth"], edgecolor=default_region[k]["color"],
                              facecolor=default_region[k]["facecolor"], zorder=default_region[k]["zorder"],
                              transform=projection)
                ax.add_patch(pol)
        # colorbar
        if legend_position == "bottom":
            x1, x2 = ax.get_position().x0, ax.get_position().x1
            y1, y2 = ax.get_position().y0, ax.get_position().y1
            dy = default_map["size_y"] / kwargs["y_size"]
            cax = pyplot.axes((x1, y1 - (y2 - y1) * dy / 6, x2 - x1, (y2 - y1) * dy / 20))
            cbar = pyplot.colorbar(lc, cax=cax, orientation="horizontal", ticks=map_c_tic, pad=0.3)
            cbar.ax.tick_params(labelsize=fontsize)
            cbar.set_label(map_c_nam, fontsize=fontsize, labelpad=2)
        elif legend_position == "right":
            x1, x2 = ax.get_position().x0, ax.get_position().x1
            y1, y2 = ax.get_position().y0, ax.get_position().y1
            dx = default_map["size_x"] / kwargs["x_size"]
            # cax = pyplot.axes((x2 + (x2 - x1) * dx / 20, y1, (x2 - x1) * dx / 30, y2 - y1))
            cax = pyplot.axes((x2 + (x2 - x1) * dx / 50, y1, (x2 - x1) * dx / 50, y2 - y1))
            cbar = pyplot.colorbar(lc, cax=cax, orientation="vertical", ticks=map_c_tic, pad=0.35)
            cbar.ax.tick_params(labelsize=fontsize)
            cbar.set_label(map_c_nam, fontsize=fontsize, rotation=90)
# ---------------------------------------------------------------------------------------------------------------------#
