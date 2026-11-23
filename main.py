# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Main program to generate figures
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# local functions
import enso_change.script_figures as fig
# ---------------------------------------------------#


# collect the figure scripts
figure_scripts = {
    "f1a": fig.f01a_remaining_iv_plot,
    "f1b": fig.f01b_remaining_iv_plot,
    "f2a": fig.f02a_change_plot,
    "f2b": fig.f02b_change_plot,
    "f3a": fig.f03a_bias_vs_change_plot,
    "f3b": fig.f03b_bias_vs_change_plot,
    "f3c": fig.f03c_bias_vs_change_plot,
    "s1": fig.s01_change_agreement_plot,
    "s2": fig.s02_signal_and_noise_plot,
    "s3": fig.s03_bias_vs_change_plot,
}
figure_names = sorted(list(figure_scripts.keys()), key=lambda v: v.lower())
figure_names_print = ", ".join(figure_names)


if __name__ == '__main__':
    # figure_number = input("Which figure do you want to plot?\n     Please enter one of: %s\n" % figure_names_print)
    figure_number = "f3c"
    while figure_number not in figure_names + ["all"]:
        figure_number = input("Given value %s does not correspond to a figure\n     Please enter one of: %s\n" % (
            figure_number, figure_names_print))
    if figure_number == "all":
        for k in figure_names:
            figure_scripts[k]()
    else:
        figure_scripts[figure_number]()
