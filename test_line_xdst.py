# Moore configuration
from Moore import options, run_moore

from Hlt2Conf.lines.lambdab_to_lambdacpi import lbtolcpi_lctopkpi_line

from RecoConf.reconstruction_objects import upfront_reconstruction, reconstruction


# Temporary workaround for TrackStateProvider
from RecoConf.global_tools import stateProvider_with_simplified_geom
public_tools = [stateProvider_with_simplified_geom()]

input_files = ['root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00111132/0000/00111132_00000075_2.xdst']


options.input_files = input_files

def all_lines():
    return [lbtolcpi_lctopkpi_line()]


#options.set_input_and_conds_from_testfiledb('Upgrade_MinBias_LDST')
options.input_type = 'ROOT'
options.dddb_tag = 'dddb-20190223'
options.conddb_tag = 'sim-20180530-vc-md100'
options.input_raw_format = 4.3
#options.evt_max = 100
#options.control_flow_file = 'control_flow.gv'
#options.data_flow_file = 'data_flow.gv'

run_moore(options, all_lines, public_tools)

