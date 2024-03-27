#Moore configurations
from Moore import options, run_moore
# from HltEfficiencyChecker.config import run_moore_with_tuples
from Hlt2Conf.lines.bnoc import all_lines

from RecoConf.global_tools import stateProvider_with_simplified_geom, trackMasterExtrapolator_with_simplified_geom
from RecoConf.reconstruction_objects import reconstruction
from RecoConf.hlt2_global_reco import (
    reconstruction as hlt2_reconstruction,
    make_light_reco_pr_kf_without_UT,
)
from RecoConf.ttrack_selections_reco import make_ttrack_reco
from RecoConf.calorimeter_reconstruction import make_digits
from RecoConf.hlt1_muonid import make_muon_hits

from Gaudi.Configuration import FileCatalog
from GaudiConf import IOHelper

options.output_file = 'test_b_to_kpi_MC_2403.ldst'
options.output_type = 'ROOT'
options.output_manifest_file = "test_b_to_kpi_MC_2403.tck.json"

options.set_input_and_conds_from_testfiledb("upgrade_minbias_hlt1_filtered")

options.persistreco_version = 0.0
options.input_raw_format = 4.3
make_digits.global_bind(calo_raw_bank=False)
options.evt_max = 2000

public_tools = [
    stateProvider_with_simplified_geom(),
]

def make_lines():
    lines = [all_lines['Hlt2BnoC_BdsToKSPi0_LL']()]
    # lines = [all_lines['Hlt2BnoC_BuToKpPi0']()]
    return lines
    
with reconstruction.bind(from_file=False):    
    config = run_moore(options, make_lines, public_tools)
