from Moore import options, run_moore

from Hlt2Conf.lines.b_to_kpi.lines import BToKpi0_line
from RecoConf.global_tools import stateProvider_with_simplified_geom
from RecoConf.reconstruction_objects import upfront_reconstruction, reconstruction

def all_lines():
    return [BToKpi0_line()]

input_files = [
    #'b_to_kpi_events/00127044_00000011_1.xdigi'
    'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000011_1.xdigi'
]

#options.set_input_and_conds_from_testfiledb('Upgrade_MinBias_LDST')
options.input_files = input_files
options.input_type = "ROOT"
options.input_raw_format = 0.3

#options.set_input_and_conds_from_testfiledb('b_to_kpi_events/00127044_00000011_1.xdigi')                                                          
options.evt_max = 100                                                                       
options.dddb_tag = 'dddb-20171126' #?
options.conddb_tag = 'sim-20171127-vc-md100' #?

#options.control_flow_file = 'control_flow.gv'

public_tools = [stateProvider_with_simplified_geom()]

with reconstruction.bind(from_file=False):
    run_moore(options, all_lines, public_tools)

