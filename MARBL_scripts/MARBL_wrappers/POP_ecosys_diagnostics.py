""" Scripts to generate POP-specific diagnostic list in generic MARBL format
    (e.g. MARBL tracer state)
"""

def write_ecosys_diagnostics_file(active_tracers, autotroph_list, zooplankton_list, calcifier_list, ladjust_bury_coeff, ecosys_diag_filename):
    """ Subroutine to write a file in the same format as marbl_diagnostics containing
        a list of POP-generated diagnostics that should be included based on the
        MARBL configuration
    """

    fout = open(ecosys_diag_filename,"w")
    # File header with information on how to use generated file
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

    # File will contain POP and MARBL diagnostics, so we provide header to make
    # the provenance of each diagnostic clear
    fout.write("#\n########################################\n")
    fout.write("#       POP-generated diagnostics      #\n")
    fout.write("########################################\n#\n")

    # Add tracer-agnostic forcing fields to requested diagnostics
    fout.write("# Dust and Carbon Fluxes from the Coupler\n#\n")
    fout.write("ATM_FINE_DUST_FLUX_CPL : medium_average\n")
    fout.write("ATM_COARSE_DUST_FLUX_CPL : medium_average\n")
    fout.write("SEAICE_DUST_FLUX_CPL : medium_average\n")
    fout.write("ATM_BLACK_CARBON_FLUX_CPL : medium_average\n")
    fout.write("SEAICE_BLACK_CARBON_FLUX_CPL : medium_average\n")

    # If adjusting bury coefficients, add running means to requested diagnostics
    if ladjust_bury_coeff:
        fout.write("#\n# Running means computed for MARBL\n#\n")
        fout.write("MARBL_rmean_glo_scalar_POC_bury_coeff : medium_average\n")
        fout.write("MARBL_rmean_glo_scalar_POP_bury_coeff : medium_average\n")
        fout.write("MARBL_rmean_glo_scalar_bSi_bury_coeff : medium_average\n")

    # 1. Create dictionary with default tracer output for all tracers
    #    - This dictionary also stores some tracer properties (currently just for budget-specific diagnostics)
    #    NOTE: using OrderedDict to maintain alphabetical listing of tracers in keys()
    from collections import OrderedDict
    full_diag_dict = OrderedDict()
    for tracer_short_name in sorted(active_tracers):
        per_tracer_dict = dict()

        # Properties used to determine frequency of budget terms
        per_tracer_dict['properties'] = dict()
        per_tracer_dict['properties']['include budget terms'] = False
        per_tracer_dict['properties']['has surface flux'] = False

        # Default frequencies for per-tracer diagnostics
        # - tracer state should be output monthly
        # - everything else is off by default
        per_tracer_dict['diags'] = OrderedDict()
        per_tracer_dict['diags'][tracer_short_name] = 'medium_average'
        per_tracer_dict['diags']['STF_%s' % tracer_short_name] = 'never_average'
        per_tracer_dict['diags']['J_%s' % tracer_short_name] = 'never_average'
        per_tracer_dict['diags']['Jint_100m_%s' % tracer_short_name] = 'never_average'
        per_tracer_dict['diags']['Jint_%s' % tracer_short_name] = 'never_average'
        per_tracer_dict['diags']['%s_zint_100m' % tracer_short_name] = 'never_average'
        per_tracer_dict['diags']['tend_zint_100m_%s' % tracer_short_name] = 'never_average'

        # Some diagnostics are not defined for all tracers; diagnostics with 'none'
        # are not added to diagnostics file and will not show up in tavg_contents
        per_tracer_dict['diags']['%s_RIV_FLUX' % tracer_short_name] = 'none'
        per_tracer_dict['diags']['FvPER_%s' % tracer_short_name] = 'none'
        per_tracer_dict['diags']['FvICE_%s' % tracer_short_name] = 'none'
        full_diag_dict[tracer_short_name] = dict(per_tracer_dict)

    # 2. Update dictionary for tracers that don't just rely on default diagnostic output
    #    This is organized per-tracer, and specific blocks are ignored if MARBL is not
    #    configured to run with that particular tracer
    # PO4
    if 'PO4' in full_diag_dict.keys():
        full_diag_dict['PO4']['diags']['PO4_RIV_FLUX'] = 'medium_average'
        full_diag_dict['PO4']['diags']['J_PO4'] = 'low_average'
        full_diag_dict['PO4']['diags']['Jint_100m_PO4'] = 'medium_average'
        full_diag_dict['PO4']['diags']['tend_zint_100m_PO4'] = 'medium_average'
        full_diag_dict['PO4']['properties']['has surface flux'] = True
    # NO3
    if 'NO3' in full_diag_dict.keys():
        full_diag_dict['NO3']['diags']['NO3_RIV_FLUX'] = 'medium_average'
        full_diag_dict['NO3']['diags']['J_NO3'] = 'low_average'
        full_diag_dict['NO3']['diags']['Jint_100m_NO3'] = 'medium_average'
        full_diag_dict['NO3']['diags']['tend_zint_100m_NO3'] = 'medium_average'
        full_diag_dict['NO3']['properties']['has surface flux'] = True
    # SiO3
    if 'SiO3' in full_diag_dict.keys():
        full_diag_dict['SiO3']['diags']['SiO3_RIV_FLUX'] = 'medium_average'
        full_diag_dict['SiO3']['diags']['J_SiO3'] = 'low_average'
        full_diag_dict['SiO3']['diags']['Jint_100m_SiO3'] = 'medium_average'
        full_diag_dict['SiO3']['diags']['tend_zint_100m_SiO3'] = 'medium_average'
        full_diag_dict['SiO3']['properties']['has surface flux'] = True
    # NH4
    if 'NH4' in full_diag_dict.keys():
        full_diag_dict['NH4']['diags']['J_NH4'] = 'low_average'
        full_diag_dict['NH4']['diags']['Jint_100m_NH4'] = 'medium_average'
        full_diag_dict['NH4']['diags']['tend_zint_100m_NH4'] = 'medium_average'
        full_diag_dict['NH4']['properties']['has surface flux'] = True
    # Fe
    if 'Fe' in full_diag_dict.keys():
        full_diag_dict['Fe']['diags']['Fe_RIV_FLUX'] = 'medium_average'
        full_diag_dict['Fe']['diags']['J_Fe'] = 'low_average'
        full_diag_dict['Fe']['diags']['Jint_100m_Fe'] = 'medium_average'
        full_diag_dict['Fe']['diags']['tend_zint_100m_Fe'] = 'medium_average'
        full_diag_dict['Fe']['properties']['include budget terms'] = True
        full_diag_dict['Fe']['properties']['has surface flux'] = True
    # Lig
    if 'Lig' in full_diag_dict.keys():
        pass # Lig just uses default settings
    # O2
    if 'O2' in full_diag_dict.keys():
        full_diag_dict['O2']['diags']['STF_O2'] = 'medium_average, high_average'
        full_diag_dict['O2']['diags']['Jint_100m_O2'] = 'medium_average'
        full_diag_dict['O2']['diags']['tend_zint_100m_O2'] = 'medium_average'
        full_diag_dict['O2']['properties']['include budget terms'] = True
        full_diag_dict['O2']['properties']['has surface flux'] = True
    # DIC
    if 'DIC' in full_diag_dict.keys():
        full_diag_dict['DIC']['diags']['DIC_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DIC']['diags']['J_DIC'] = 'medium_average'
        full_diag_dict['DIC']['diags']['Jint_100m_DIC'] = 'medium_average'
        full_diag_dict['DIC']['diags']['tend_zint_100m_DIC'] = 'medium_average'
        full_diag_dict['DIC']['diags']['FvPER_DIC'] = 'medium_average'
        full_diag_dict['DIC']['diags']['FvICE_DIC'] = 'medium_average'
        full_diag_dict['DIC']['properties']['include budget terms'] = True
        full_diag_dict['DIC']['properties']['has surface flux'] = True
    # DIC_ALT_CO2
    if 'DIC_ALT_CO2' in full_diag_dict.keys():
        full_diag_dict['DIC_ALT_CO2']['diags']['DIC_ALT_CO2_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DIC_ALT_CO2']['diags']['J_DIC_ALT_CO2'] = 'medium_average'
        full_diag_dict['DIC_ALT_CO2']['diags']['Jint_100m_DIC_ALT_CO2'] = 'medium_average'
        full_diag_dict['DIC_ALT_CO2']['diags']['tend_zint_100m_DIC_ALT_CO2'] = 'medium_average'
        full_diag_dict['DIC_ALT_CO2']['diags']['FvPER_DIC_ALT_CO2'] = 'medium_average'
        full_diag_dict['DIC_ALT_CO2']['diags']['FvICE_DIC_ALT_CO2'] = 'medium_average'
        full_diag_dict['DIC_ALT_CO2']['properties']['include budget terms'] = True
        full_diag_dict['DIC_ALT_CO2']['properties']['has surface flux'] = True
    # ALK
    if 'ALK' in full_diag_dict.keys():
        full_diag_dict['ALK']['diags']['ALK_RIV_FLUX'] = 'medium_average'
        full_diag_dict['ALK']['diags']['STF_ALK'] = 'medium_average'
        full_diag_dict['ALK']['diags']['J_ALK'] = 'low_average'
        full_diag_dict['ALK']['diags']['Jint_100m_ALK'] = 'medium_average'
        full_diag_dict['ALK']['diags']['tend_zint_100m_ALK'] = 'medium_average'
        full_diag_dict['ALK']['diags']['FvPER_ALK'] = 'medium_average'
        full_diag_dict['ALK']['diags']['FvICE_ALK'] = 'medium_average'
        full_diag_dict['ALK']['properties']['has surface flux'] = True
    # ALK_ALT_CO2
    if 'ALK_ALT_CO2' in full_diag_dict.keys():
        full_diag_dict['ALK_ALT_CO2']['diags']['ALK_ALT_CO2_RIV_FLUX'] = 'medium_average'
        full_diag_dict['ALK_ALT_CO2']['diags']['STF_ALK_ALT_CO2'] = 'medium_average'
        full_diag_dict['ALK_ALT_CO2']['diags']['J_ALK_ALT_CO2'] = 'low_average'
        full_diag_dict['ALK_ALT_CO2']['diags']['Jint_100m_ALK_ALT_CO2'] = 'medium_average'
        full_diag_dict['ALK_ALT_CO2']['diags']['tend_zint_100m_ALK_ALT_CO2'] = 'medium_average'
        full_diag_dict['ALK_ALT_CO2']['diags']['FvPER_ALK_ALT_CO2'] = 'medium_average'
        full_diag_dict['ALK_ALT_CO2']['diags']['FvICE_ALK_ALT_CO2'] = 'medium_average'
        full_diag_dict['ALK_ALT_CO2']['properties']['has surface flux'] = True
    # DOC
    if 'DOC' in full_diag_dict.keys():
        full_diag_dict['DOC']['diags']['DOC_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DOC']['diags']['Jint_100m_DOC'] = 'medium_average'
        full_diag_dict['DOC']['diags']['tend_zint_100m_DOC'] = 'medium_average'
        full_diag_dict['DOC']['properties']['include budget terms'] = True
        full_diag_dict['DOC']['properties']['has surface flux'] = False # this should be True if EBM is off
    # DON
    if 'DON' in full_diag_dict.keys():
        full_diag_dict['DON']['diags']['DON_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DON']['properties']['has surface flux'] = True
    # DOP
    if 'DOP' in full_diag_dict.keys():
        full_diag_dict['DOP']['diags']['DOP_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DOP']['properties']['has surface flux'] = True
    # DOPr
    if 'DOPr' in full_diag_dict.keys():
        full_diag_dict['DOPr']['diags']['DOPr_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DOPr']['properties']['has surface flux'] = True
    # DONr
    if 'DONr' in full_diag_dict.keys():
        full_diag_dict['DONr']['diags']['DONr_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DONr']['properties']['has surface flux'] = True
    # DOCr
    if 'DOCr' in full_diag_dict.keys():
        full_diag_dict['DOCr']['diags']['DOCr_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DOCr']['diags']['Jint_100m_DOCr'] = 'medium_average'
        full_diag_dict['DOCr']['diags']['tend_zint_100m_DOCr'] = 'medium_average'
        full_diag_dict['DOCr']['properties']['include budget terms'] = True
        full_diag_dict['DOCr']['properties']['has surface flux'] = False # this should be True if EBM is off
    # DI13C
    if 'DI13C' in full_diag_dict.keys():
        full_diag_dict['DI13C']['diags']['DI13C_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DI13C']['diags']['J_DI13C'] = 'medium_average'
        full_diag_dict['DI13C']['diags']['Jint_100m_DI13C'] = 'medium_average'
        full_diag_dict['DI13C']['diags']['tend_zint_100m_DI13C'] = 'medium_average'
        full_diag_dict['DI13C']['diags']['FvPER_DI13C'] = 'medium_average'
        full_diag_dict['DI13C']['diags']['FvICE_DI13C'] = 'medium_average'
        full_diag_dict['DI13C']['properties']['has surface flux'] = True
    # DO13Ctot
    if 'DO13Ctot' in full_diag_dict.keys():
        full_diag_dict['DO13Ctot']['diags']['DO13Ctot_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DO13Ctot']['diags']['Jint_100m_DO13Ctot'] = 'medium_average'
        full_diag_dict['DO13Ctot']['diags']['tend_zint_100m_DO13Ctot'] = 'medium_average'
        full_diag_dict['DO13Ctot']['properties']['has surface flux'] = True
    # DI14C
    if 'DI14C' in full_diag_dict.keys():
        full_diag_dict['DI14C']['diags']['DI14C_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DI14C']['diags']['J_DI14C'] = 'medium_average'
        full_diag_dict['DI14C']['diags']['Jint_100m_DI14C'] = 'medium_average'
        full_diag_dict['DI14C']['diags']['tend_zint_100m_DI14C'] = 'medium_average'
        full_diag_dict['DI14C']['diags']['FvPER_DI14C'] = 'medium_average'
        full_diag_dict['DI14C']['diags']['FvICE_DI14C'] = 'medium_average'
        full_diag_dict['DI14C']['properties']['has surface flux'] = True
    # DO14Ctot
    if 'DO14Ctot' in full_diag_dict.keys():
        full_diag_dict['DO14Ctot']['diags']['DO14Ctot_RIV_FLUX'] = 'medium_average'
        full_diag_dict['DO14Ctot']['diags']['Jint_100m_DO14Ctot'] = 'medium_average'
        full_diag_dict['DO14Ctot']['diags']['tend_zint_100m_DO14Ctot'] = 'medium_average'
        full_diag_dict['DO14Ctot']['properties']['has surface flux'] = True

    # 3. Per-autotroph diagnostics
    for autotroph_name in autotroph_list:
        tracer_short_name = autotroph_name+'C'
        if tracer_short_name in full_diag_dict.keys():
            full_diag_dict[tracer_short_name]['diags']['%s_zint_100m' % tracer_short_name] = 'high_average'
        tracer_short_name = autotroph_name+'CaCO3'
        if tracer_short_name in full_diag_dict.keys():
            full_diag_dict[tracer_short_name]['diags']['%s_zint_100m' % tracer_short_name] = 'high_average'
        tracer_short_name = autotroph_name+'Chl'
        if tracer_short_name in full_diag_dict.keys():
            full_diag_dict[tracer_short_name]['diags']['%s_SURF' % tracer_short_name] = 'high_average'

    # 4. Per-zooplankton diagnostics
    for zooplankton_name in zooplankton_list:
        tracer_short_name = zooplankton_name+'C'
        if tracer_short_name in full_diag_dict.keys():
            full_diag_dict[tracer_short_name]['diags']['%s_zint_100m' % tracer_short_name] = 'high_average'

    # 5. Write tracer-specific diagnostics to file
    for tracer_short_name in full_diag_dict.keys():
        # (a) Process ['properties'] dictionary for budget terms
        #     Note that this step will add to the ['diags'] dictionary but will not change previously set values
        per_tracer_properties = full_diag_dict[tracer_short_name]['properties']
        if per_tracer_properties['include budget terms']:
            value = "low_average"
        else:
            value = "never_average"

        # These diagnostics should be included by default for tracers requested budget terms
        for key in ['UE', 'VN', 'WT', 'DIA_IMPVF', 'HDIFE', 'HDIFN', 'HDIFB', 'TEND', 'RF_TEND']:
            specific_key = '%s_%s' % (key, tracer_short_name)
            if specific_key not in full_diag_dict[tracer_short_name]['diags'].keys():
                full_diag_dict[tracer_short_name]['diags'][specific_key] = value

        if per_tracer_properties['include budget terms'] and per_tracer_properties['has surface flux']:
            value = "low_average"
        else:
            value = "never_average"

        # These diagnostics should be included by default for tracers requested budget terms
        # ONLY for tracers that have non-zero surface fluxes
        for key in ['KPP_SRC']:
            specific_key = '%s_%s' % (key, tracer_short_name)
            if specific_key not in full_diag_dict[tracer_short_name]['diags'].keys():
                full_diag_dict[tracer_short_name]['diags'][specific_key] = value

        # (b) Loop through ['diags'] dictionary and write all diagnostics to file
        fout.write("#\n# Diagnostics for tracer %s\n#\n" % tracer_short_name)
        for diag in full_diag_dict[tracer_short_name]['diags'].keys():
            per_tracer_dict = full_diag_dict[tracer_short_name]['diags']
            if per_tracer_dict[diag] != 'none':
                fout.write("%s : %s\n" % (diag, per_tracer_dict[diag]))

    # 6. Add section header for MARBL diagnostics
    #    (Another tool appends MARBL diagnostics to this file)
    fout.write("#\n########################################\n")
    fout.write("#      MARBL-generated diagnostics     #\n")
    fout.write("########################################\n#\n")
