# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default_arguments = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical", "piControl"],
    # dictionary of desired observational datasets per diagnostic: dict[str, list[str]]
    "data_observational_datasets": {
        "ave_pr_val_n30e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ave_pr_val_n34e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ave_pr_val_n40e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ave_pr_val_glob": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ave_ts_val_n30e": ["COBE2", "ERSSTv5", "HadISST"],
        "ave_ts_val_n34e": ["COBE2", "ERSSTv5", "HadISST"],
        "ave_ts_val_n40e": ["COBE2", "ERSSTv5", "HadISST"],
        "ave_ts_val_glob": ["ERA5"],
        "cor_ts_n30e_to_tx_n40e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "cor_ts_n34e_to_tx_n40e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "fbk_ts_n30e_to_tx_n40e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "fbk_ts_n34e_to_tx_n40e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ske_pr_ano_n30e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ske_pr_ano_n34e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ske_pr_ano_n40e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ske_pr_ano_glob": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        # "ske_ts_ano_n30e": ["COBE2", "ERSSTv5", "HadISST"],
        # "ske_ts_ano_n34e": ["COBE2", "ERSSTv5", "HadISST"],
        # "ske_ts_ano_n40e": ["COBE2", "ERSSTv5", "HadISST"],
        "ske_ts_ano_n30e": ["COBE2", "HadISST"],
        "ske_ts_ano_n34e": ["COBE2", "HadISST"],
        "ske_ts_ano_n40e": ["COBE2", "HadISST"],
        "ske_ts_ano_glob": ["ERA5"],
        "std_pr_ano_n30e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "std_pr_ano_n34e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "std_pr_ano_n40e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "std_pr_ano_glob": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        # "std_ts_ano_n30e": ["COBE2", "ERSSTv5", "HadISST"],
        # "std_ts_ano_n34e": ["COBE2", "ERSSTv5", "HadISST"],
        # "std_ts_ano_n40e": ["COBE2", "ERSSTv5", "HadISST"],
        "std_ts_ano_n30e": ["COBE2", "HadISST"],
        "std_ts_ano_n34e": ["COBE2", "HadISST"],
        "std_ts_ano_n40e": ["COBE2", "HadISST"],
        "std_ts_ano_glob": ["ERA5"],
    },
    # list of projects: list[str]
    "data_projects": ["observations", "cmip6"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 5,
    # list of rejected SMILEs: list[str]
    #      - CAS-ESM2-0 has clear change in precipitation variability between piControl and historical runs
    #      - KACE-1-0-G's GMST from the piControl run has a clear trend of ~0.33°C/century throughout the simulation
    # "smile_rejected": ["EC-EARTH", "KACE-1-0-G"],
    "smile_rejected": ["CMCC-CM", "CMCC-CESM", "EC-EARTH", "GISS-E2-2-H", "KACE-1-0-G"],
    # "smile_rejected": ["CMCC-CM", "EC-EARTH", "KACE-1-0-G"],
    # "smile_rejected": ["CMCC-CM", "KACE-1-0-G"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # colors
    "fig_colors_dat": {
        # observations
        "COBE2": "turquoise",
        "ERSSTv5": "y",
        "GFDL-ECDAv3.1": "lightseagreen",
        "HadISST": "b",
        "ORAS5": "forestgreen",
    },
    "fig_colors_exp": {
        "historical": "forestgreen",
        # "piControl": "violet",
        "piControl": "b",
        "rcp26": "gold",
        "rcp45": "orange",
        "rcp60": "peru",
        "rcp85": "magenta",
        "scenario19": "darkgrey",
        "scenario26": "gold",
        "scenario45": "orange",
        "scenario60": "peru",
        "scenario70": "firebrick",
        "scenario85": "magenta",
        "ssp119": "darkgrey",
        "ssp126": "gold",
        "ssp245": "orange",
        "ssp370": "firebrick",
        "ssp585": "magenta",
    },
    # linestyles
    "fig_linestyles_dat": {
        # observations
        "COBE2": "-",
        "ERSSTv5": "--",
        "GFDL-ECDAv3.1": "--",
        "HadISST": "-.",
        "ORAS5": "-",
    },
    "fig_linestyles_exp": {
        "historical": "-",
        "piControl": "-",
        "rcp26": "-",
        "rcp45": "-",
        "rcp60": "-",
        "rcp85": "-",
        "scenario19": "-",
        "scenario26": "-",
        "scenario45": "-",
        "scenario60": "-",
        "scenario70": "-",
        "scenario85": "-",
        "ssp119": "-",
        "ssp126": "-",
        "ssp245": "-",
        "ssp370": "-",
        "ssp585": "-",
    },
    # markers
    "fig_markers_dat": {
        # observations
        "COBE2": "s",
        "ERSSTv5": "p",
        "HadISST": "D",
    },
    # names
    "fig_names_exp": {
        "historical": "Historical",  # "Historical",  #
        "piControl": "piControl",  # "piControl",  #
        "rcp26": "RCP2.6",
        "rcp45": "RCP4.5",
        "rcp60": "RCP6.0",
        "rcp85": "RCP8.5",
        "scenario19": "R/S1.9",
        "scenario26": "R/S2.6",
        "scenario45": "R/S4.5",
        "scenario60": "R/S6.0",
        "scenario70": "R/S7.0",
        "scenario85": "R/S8.5",
        "ssp119": "SSP1-1.9",
        "ssp126": "SSP1-2.6",
        "ssp245": "SSP2-4.5",
        "ssp370": "SSP3-7.0",
        "ssp585": "SSP5-8.5",
    },
    "fig_names_reg": {
        "glob": "global",
        "n30e": "Niño3",
        "n34e": "Niño3.4",
        "n40e": "Niño4",
    },
}
