# Moore configuration
from Moore import options, run_moore

from Hlt2Conf.lines.b_to_kpi import BToKpi0_line


# In a normal options file, we would import the line from Hlt2Conf where it is
# defined
# from Hlt2Conf.lines.LbToLcPi import lbtolcpi_lctopkpi_line

# Temporary workaround for TrackStateProvider
from RecoConf.global_tools import stateProvider_with_simplified_geom
from RecoConf.reconstruction_objects import upfront_reconstruction, reconstruction

public_tools = [stateProvider_with_simplified_geom()]

options.set_input_and_conds_from_testfiledb('Upgrade_MinBias_LDST')
options.input_raw_format = 4.3
#options.evt_max = 100
#options.control_flow_file = 'control_flow.gv'
#options.data_flow_file = 'data_flow.gv'

def all_lines():
    return [BToKpi0_line()]
    
run_moore(options, all_lines, public_tools)
    
