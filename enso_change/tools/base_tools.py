# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Basic functions
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from glob import iglob as glob__iglob
from json import load as json__load
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join

# local functions
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def read_json() -> dict:
    """
    Read json files

    Output:
    -------
    :return dict_o: dict
        Dictionary with nine nested levels [metric, diagnostic or metadata, value or metadata_name, project, dataset,
        experiment, member, epoch_length, epoch], filled with a value
    """
    # data directory (relative to current file directory)
    data_directory = "/".join(os__path__dirname(__file__).split("/")[:-2])
    # path to input data file
    dict_o = {}
    # loop on statistic groups
    # for k1 in ["pr", "ts", "tx"]:
    #     json_file_path = os__path__join(data_directory,
    #                                     "data/ENSO_and_pacific_evolution_with_cc_300y_" + str(k1) + "_clean.json")
    #     if len(list(glob__iglob(json_file_path))) == 1:
    #         # load data
    #         with open(json_file_path) as ff:
    #             dict_t = json__load(ff)
    #         ff.close()
    #         # loop on diagnostics
    #         for k3 in list(dict_t["RESULTS"].keys()):
    #             dict_o[k3] = dict_t["RESULTS"][k3]
    # for k1 in ["ave", "ske", "std"]:
    for k1 in ["ave", "fbk", "ske", "std"]:
        for k2 in ["hf", "pr", "ts", "tx"]:
            json_file_path = os__path__join(
                data_directory, "data/enso_and_pacific_iv_" + str(k1) + "_" + str(k2) + ".json")
            if len(list(glob__iglob(json_file_path))) == 1:
                # load data
                with open(json_file_path) as ff:
                    dict_t = json__load(ff)
                ff.close()
                # loop on diagnostics
                for k3 in list(dict_t["RESULTS"].keys()):
                    dict_o[k3] = dict_t["RESULTS"][k3]
    return dict_o


def plural_s(list_i: list) -> str:
    """
    Return 's' if there are multiple values in the list

    Input:
    ------
    :param list_i: list

    Output:
    -------
    :return: str
        's' if there are multiple values in the list else ''
    """
    return "s" if isinstance(list_i, list) is True and len(list_i) > 1 else ""


def short_name(name_i: str) -> str:
    """
    Modify short name for plots

    Input:
    ------
    :param name_i: str
        Variable short name

    Output:
    -------
    :return: str
        Short name modified for plots
    """
    # replace statistic name
    name_o = name_i.replace("AVE", r"$\bar{x}$").replace("SKE", "g$_1$")
    name_o = name_o.replace("STD", r"$\sigma$").replace(" TRE", r" $\partial_{t}$")
    name_o = name_o.replace("VAR", r"$\sigma^2$")
    name_o = name_o.replace("COR", r"$\rho$").replace("REG", "b$_1$")
    # replace region name (nino3, nino3.4, nino4)
    # name_o = name_o.replace("N30", "Niño3").replace("N34", "Niño3.4")
    # name_o = name_o.replace("N40", "Niño4")
    # for k1 in ["Niño3", "Niño3.4", "Niño4"]:
    name_o = name_o.replace("N30", "N3").replace("N34", "N3.4")
    name_o = name_o.replace("N40", "N4")
    for k1 in ["N3", "N3.4", "N4"]:
        for k2 in ["E", "N", "S"]:
            if str(k1) + str(k2) in name_o:
                if k2 == "E":
                    name_o = name_o.replace(str(k1) + str(k2), k1)
                else:
                    name_o = name_o.replace(str(k1) + str(k2), str(k1) + "$_{" + str(k2) + "}$")
    # replace region name (other regions)
    name_o = name_o.replace("PENO", "EP$_N$").replace("PEQU", "EP")
    name_o = name_o.replace("PESO", "EP$_S$")
    name_o = name_o.replace("PEQE", "EEP").replace("PEQW", "WEP")
    name_o = name_o.replace("PTRO", "TP")
    # name_o = name_o.replace("GXIP", "$\Delta_{Ind-Pac}$")
    # name_o = name_o.replace("GXPB", "$\Delta_{N4-N3}$ [8°S-N]")
    # name_o = name_o.replace("GXPC", "$\Delta_{W-C}$")
    # name_o = name_o.replace("GXPD", "$\Delta_{W-C}$ [5°S-N]")
    # name_o = name_o.replace("GXPE", "$\Delta_{N4-N3}$")
    # name_o = name_o.replace("GXPL", "$\Delta_{W-E}$ (W21)")
    # name_o = name_o.replace("GXPN", "$\Delta_{N4-N3}$ [2.5°S-N]")
    # name_o = name_o.replace("GXPR", "$\Delta_{Pac-Tro}$")
    # name_o = name_o.replace("GXPS", "$\Delta_{W-E}$ (S22)")
    # name_o = name_o.replace("GXPW", "$\Delta_{W-E}$ (H21)")
    # name_o = name_o.replace("dzT", "$\Delta_{surf-sub}$T")
    name_o = name_o.replace("GXIP", "$\Delta_{Ind-Pac}$")
    name_o = name_o.replace("GXPB", "$\Delta_{x}$ [8°S-N]")
    name_o = name_o.replace("GXPC", "$\Delta_{W-C}$")
    name_o = name_o.replace("GXPD", "$\Delta_{W-C}$ [5°S-N]")
    name_o = name_o.replace("GXPE", "$\Delta_{x}$ [5°S-N]")
    name_o = name_o.replace("GXPL", "$\Delta_{x}$ (W21)")
    name_o = name_o.replace("GXPN", "$\Delta_{x}$ [2.5°S-N]")
    name_o = name_o.replace("GXPR", "$\Delta_{Pac-Tro}$")
    name_o = name_o.replace("GXPS", "$\Delta_{x}$ (S22)")
    name_o = name_o.replace("GXPW", "$\Delta_{x}$ (H21)")
    name_o = name_o.replace("dzT", "$\Delta_{z}$T")
    # variable
    if "TS" in name_o and "GLOB" not in name_o:
        name_o = name_o.replace("TS", "SST")
        name_o = name_o.replace("ts", "SST")
    if "sst" in name_o:
        name_o = name_o.replace("sst", "SST")
    if "ts" in name_o:
        name_o = name_o.replace("ts", "TS")
    return name_o


def sort_members(dataset: str, list_members: list) -> list:
    """
    Sort members

    Input:
    ------
    :param dataset: str
        Dataset name; e.g., dataset = 'ACCESS-CM2'
    :param list_members: list
        List of members; e.g., list_members = ['r10i1p1f1', 'r1i1p1f1', 'r20i1p1f1', 'r2i1p1f1']

    Output:
    -------
    :return list_o: list
        Sorted list of members; e.g., list_o = ['r1i1p1f1', 'r2i1p1f1', 'r10i1p1f1', 'r20i1p1f1']
    """
    # sort members
    list_o = copy__deepcopy(list_members)
    if dataset == "MPI-ESM" and len(list_members) > 1:
        list_o = sorted(list_members, key=lambda s: s.lower())
    elif len(list_members) > 1:
        # pad every variant-ID with 6 zeros
        list_t = []
        for k in list_members:
            # realisation (i.e. ensemble member)
            rr = str(k.split("i")[0].replace("r", "")).zfill(6)
            # initialisation method
            ii = str(k.split("i")[1].split("p")[0]).zfill(6)
            # physics
            pp = str(k.split("p")[1].split("f")[0]).zfill(6)
            # forcing
            if "f" in k:
                ff = "f" + str(k.split("f")[1]).zfill(6)
            else:
                ff = ""
            # save new member name
            list_t.append("r" + str(rr) + "i" + str(ii) + "p" + str(pp) + str(ff))
        # sort padded member names
        list_t = sorted(list_t, key=lambda s: s.lower())
        # remove pads
        list_o = []
        for k in list_t:
            tmp = copy__deepcopy(k)
            while "r0" in tmp:
                tmp = tmp.replace("r0", "r")
            while "i0" in tmp:
                tmp = tmp.replace("i0", "i")
            while "p0" in tmp:
                tmp = tmp.replace("p0", "p")
            while "f0" in tmp:
                tmp = tmp.replace("f0", "f")
            list_o.append(tmp)
    return list_o


def units(name_i: str) -> str:
    """
    Modify units for plots

    Input:
    ------
    :param name_i: str
        Variable units

    Output:
    -------
    :return: str
        Units modified for plots
    """
    name_o = name_i.replace("degC", "$^\circ$C")
    name_o = name_o.replace(" / cent", ".c**-1")
    name_o = name_o.replace("degree_north", "$^\circ$N")
    for k in range(-20, 21):
        name_o = name_o.replace("**" + str(k), "$^{" + str(k) + "}$")
    return name_o
# ---------------------------------------------------------------------------------------------------------------------#
