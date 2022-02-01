# Moore configuration
from Moore import options

from HltEfficiencyChecker.config  import run_moore_with_tuples

# Temporary workaround for TrackStateProvider
from RecoConf.global_tools import stateProvider_with_simplified_geom
from RecoConf.reconstruction_objects import upfront_reconstruction, reconstruction

options.set_input_and_conds_from_testfiledb('upgrade_minbias_hlt1_filtered')
options.input_raw_format = 4.3
#options.dddb_tag = 'dddb-20171126'
#options.conddb_tag = 'sim-20171127-vc-md100'
#options.evt_max = 100

options.ntuple_file = "rate_ntuple.root"

#def all_lines():
#    return [BToKpi0_line()]

from RecoConf.hlt1_tracking import default_ft_decoding_version
default_ft_decoding_version.global_bind(value=2)

with reconstruction.bind(from_file=False):    
    run_moore_with_tuples(
        options, False, public_tools=[stateProvider_with_simplified_geom()])
    
