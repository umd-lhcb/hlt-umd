# Moore configuration
from Moore import options, run_moore
# from Moore.tcks import dump_hlt2_configuration
from Hlt2Conf.lines.bnoc.hlt2_bnoc import Bds_KSPi_LL_line, BuPiPi_line

# Temporary workaround for TrackStateProvider
from RecoConf.global_tools import stateProvider_with_simplified_geom
from RecoConf.reconstruction_objects import upfront_reconstruction, reconstruction

from Gaudi.Configuration import FileCatalog, ApplicationMgr
from GaudiConf import IOHelper

# catalog = FileCatalog().Catalogs = [ 'xmlcatalog_file:/afs/cern.ch/user/e/ejiang/my_work/public/stack/MC/catalogs/minbiasCatalog.xml' ]
# ApplicationMgr().ExtSvc.append(catalog)
# input_files = []

#options.evt_max = 1000
options.simulation = True

#options.input_files = input_files
options.input_type = 'ROOT'
#options.set_input_and_conds_from_testfiledb('Upgrade_MinBias_LDST')
options.input_raw_format = 4.3
options.data_type = 'Upgrade'
options.dddb_tag = 'dddb-20171126'
options.conddb_tag = 'sim-20171127-vc-md100'
options.geometry_version = 'trunk'
options.conditions_version = 'master'

options.output_file = 'test_b_to_kpi_minbias_tck.LDST'
options.output_type = 'ROOT'

def all_lines():
    return [Bds_KSPi_LL_line(), BuPiPi_line()]
    
public_tools = [stateProvider_with_simplified_geom()]
#run_moore(options, all_lines, public_tools)

config = run_moore(options, all_lines, public_tools)
#dump_hlt2_configuration(config, "my_test_b_to_kpi_minbias.tck.json")
options.output_manifest_file = "my_test_b_to_kpi_minbias.tck.json"    
