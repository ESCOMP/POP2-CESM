""" Scripts to generate POP-specific MARBL-based diagnostic list
    (e.g. MARBL tracer state)
"""

def write_ecosys_diagnostics_file(active_tracers, autotroph_list, zooplankton_list, calcifier_list, ladjust_bury_coeff, ecosys_diag_filename):
    """ Subroutine to write a file in the same format as marbl_diagnostics containing
        a list of POP-generated diagnostics that should be included based on the
        MARBL configuration
    """

    fout = open(ecosys_diag_filename,"w")
    # Sort variables by subcategory
    fout.write("# This file contains a list of all ecosystem-related diagnostics POP output for a given MARBL configuration,\n")
    fout.write("# as well as the recommended frequency and operator for outputting each diagnostic.\n")
    fout.write("# Some diagnostics are computed in POP, while others are provided by MARBL.\n")
    fout.write("# The format of this file is:\n")
    fout.write("#\n")
    fout.write("# DIAGNOSTIC_NAME : frequency_operator\n")
    fout.write("#\n")
    fout.write("# And fields that should be output at multiple different frequencies will be comma-separated:\n")
    fout.write("#\n")
    fout.write("# DIAGNOSTIC_NAME : frequency1_operator1, frequency2_operator2, ..., frequencyN_operatorN\n")
    fout.write("#\n")
    fout.write("# Frequencies are never, low, medium, and high.\n")
    fout.write("# Operators are instantaneous, average, minimum, and maximum.\n")
    fout.write("#\n")
    fout.write("# To change BGC-related diagnostic output, copy this file to SourceMods/src.pop/\n")
    fout.write("# and edit as desired.\n")
    fout.write("#\n########################################\n")
    fout.write("#       POP-generated diagnostics      #\n")
    fout.write("########################################\n#\n")
    # Add forcing fields
    fout.write("# River Fluxes from the Coupler\n#\n")
    fout.write("FINE_DUST_FLUX_CPL : medium_average\n")
    fout.write("COARSE_DUST_FLUX_CPL : medium_average\n")
    fout.write("BLACK_CARBON_FLUX_CPL : medium_average\n")

    # Running means
    if ladjust_bury_coeff:
        fout.write("#\n# Running means computed for MARBL\n#\n")
        fout.write("MARBL_rmean_glo_scalar_POC_bury_coeff : medium_average\n")
        fout.write("MARBL_rmean_glo_scalar_POP_bury_coeff : medium_average\n")
        fout.write("MARBL_rmean_glo_scalar_bSi_bury_coeff : medium_average\n")

    # Per-tracer diagnostics
    for tracer_short_name in active_tracers:
        fout.write("#\n# Diagnostics for tracer %s\n#\n" % tracer_short_name)
        # Tracer state
        fout.write("%s : medium_average\n" % tracer_short_name)

        # River Flux
        if tracer_short_name in ['NO3', 'PO4', 'DON', 'DONr', 'DOP', 'DOPr', 'SiO3', 'Fe', 'DIC', 'ALK', 'DOC', 'DOCr', 'DIC_ALT_CO2', 'ALK_ALT_CO2', 'DI13C', 'DI14C', 'DO13C', 'DO14C']:
            fout.write("%s_RIV_FLUX : medium_average\n" % tracer_short_name)

        # STF
        if tracer_short_name in ['ALK']:
            fout.write("STF_%s : medium_average\n" % tracer_short_name)
        elif tracer_short_name in ['O2']:
            fout.write("STF_%s : medium_average, high_average\n" % tracer_short_name)
        else:
            # Now includes CISO
            fout.write("STF_%s : never_average\n" % tracer_short_name)

        # J
        if tracer_short_name in ['DIC', 'DI13C', 'DI14C']:
            fout.write("J_%s : medium_average\n" % tracer_short_name)
        elif tracer_short_name in ['NO3', 'NH4', 'PO4', 'Fe', 'SiO3', 'ALK']:
            fout.write("J_%s : low_average\n" % tracer_short_name)
        else:
            fout.write("J_%s : never_average\n" % tracer_short_name)

        # Jint
        fout.write("Jint_%s : never_average\n" % tracer_short_name)

        # Jint_100m
        if tracer_short_name in ['DIC', 'NO3', 'NH4', 'PO4', 'Fe', 'SiO3', 'ALK', 'O2', 'DOC','DI13C', 'DI14C', 'DO13C', 'DO14C']:
            fout.write("Jint_100m_%s : medium_average\n" % tracer_short_name)
        else:
            fout.write("Jint_100m_%s : never_average\n" % tracer_short_name)

        # zint_100m
        if tracer_short_name in [PFT+'C' for PFT in autotroph_list+zooplankton_list] + [calcifier+'CaCO3' for calcifier in calcifier_list]:
            fout.write("%s_zint_100m : high_average\n" % tracer_short_name)
        else:
            fout.write("%s_zint_100m : never_average\n" % tracer_short_name)

        # tend_zint_100m
        if tracer_short_name in ['DIC', 'DIC_ALT_CO2', 'NO3', 'NH4', 'PO4', 'Fe', 'SiO3', 'ALK', 'O2', 'DOC', 'DI13C', 'DI14C', 'DO13C', 'DO14C']:
            fout.write("tend_zint_100m_%s : medium_average\n" % tracer_short_name)
        else:
            fout.write("tend_zint_100m_%s : never_average\n" % tracer_short_name)

        # Surface Chlorophyll
        if tracer_short_name in [autotroph+'Chl' for autotroph in autotroph_list]:
            fout.write("%s_SURF : high_average\n" % tracer_short_name)

        # FvPER
        if tracer_short_name in ['DIC', 'ALK', 'DIC_ALT_CO2', 'ALK_ALT_CO2', 'DI13C', 'DI14C']:
            fout.write("FvPER_%s : medium_average\n" % tracer_short_name)

        # FvICE
        if tracer_short_name in ['DIC', 'ALK', 'DI13C', 'DI14C']:
            fout.write("FvICE_%s : medium_average\n" % tracer_short_name)

        # UE, VN, WT, DIA_IMPVF, HDIFE, HDIFN, HDIFB
        if tracer_short_name in ['O2', 'DOC', 'DIC', 'DIC_ALT_CO2', 'Fe']:
            fout.write("UE_%s : low_average\n" % tracer_short_name)
            fout.write("VN_%s : low_average\n" % tracer_short_name)
            fout.write("WT_%s : low_average\n" % tracer_short_name)
            fout.write("DIA_IMPVF_%s : low_average\n" % tracer_short_name)
            fout.write("HDIFE_%s : low_average\n" % tracer_short_name)
            fout.write("HDIFN_%s : low_average\n" % tracer_short_name)
            fout.write("HDIFB_%s : low_average\n" % tracer_short_name)
        else:
            fout.write("UE_%s : never_average\n" % tracer_short_name)
            fout.write("VN_%s : never_average\n" % tracer_short_name)
            fout.write("WT_%s : never_average\n" % tracer_short_name)
            fout.write("DIA_IMPVF_%s : never_average\n" % tracer_short_name)
            fout.write("HDIFE_%s : never_average\n" % tracer_short_name)
            fout.write("HDIFN_%s : never_average\n" % tracer_short_name)
            fout.write("HDIFB_%s : never_average\n" % tracer_short_name)

        # KPP_SRC
        if tracer_short_name in ['O2', 'DIC', 'DIC_ALT_CO2', 'Fe']:
            fout.write("KPP_SRC_%s : low_average\n" % tracer_short_name)
        else:
            fout.write("KPP_SRC_%s : never_average\n" % tracer_short_name)

    # Footer before MARBL diagnostics
    fout.write("#\n########################################\n")
    fout.write("#      MARBL-generated diagnostics     #\n")
    fout.write("########################################\n#\n")