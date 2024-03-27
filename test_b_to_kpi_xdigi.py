#Moore configurations
from Moore import options, run_moore
from Hlt2Conf.lines.bnoc import all_lines

from RecoConf.global_tools import stateProvider_with_simplified_geom
from RecoConf.reconstruction_objects import upfront_reconstruction, reconstruction
# from RecoConf.hlt1_tracking import default_ft_decoding_version

from Gaudi.Configuration import FileCatalog
from GaudiConf import IOHelper

#set decoding version
# ft_decoding_version = 6
# default_ft_decoding_version.global_bind(value=ft_decoding_version)

public_tools = [stateProvider_with_simplified_geom()]

#Use this one for catalog
# input_files = ['LFN:/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000056_1.xdigi',
#               'LFN:/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000101_1.xdigi',]
#Use this one for just testing/no catalog
# input_files = ['root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000056_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000007_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000008_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000013_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000015_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000016_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000019_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000020_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000022_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000024_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000025_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000027_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000028_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000029_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000030_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000031_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000033_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000034_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000035_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000036_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000037_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000038_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000039_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000040_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000042_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000043_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000044_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000045_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000046_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000047_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000048_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000049_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000050_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000051_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000052_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000053_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000054_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000055_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000056_1.xdigi',
#                'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000057_1.xdigi',
# ]
options.input_files = input_files
options.input_type = "ROOT"
options.input_raw_format = 0.5
options.output_file = 'test_b_to_kpi_MC_2403.ldst'
options.output_type = 'ROOT'
options.output_manifest_file = "test_b_to_kpi_MC_2403.tck.json"
options.simulation = True
options.data_type = 'Upgrade'
# FileCatalog().Catalogs = ["xmlcatalog_file:/afs/cern.ch/work/e/ejiang/public/stack2/MCCatalog.xml"]
# options.evt_max = 2000
options.print_freq = 100
options.dddb_tag = 'dddb-20201211'
options.conddb_tag = 'sim-20201218-vc-md100'
# options.ntuple_file = "eff_ntuple_bnoc_hlt2.root"

from RecoConf.hlt1_muonid import make_muon_hits
make_muon_hits.global_bind(geometry_version=2)

from RecoConf.calorimeter_reconstruction import make_digits
make_digits.global_bind(calo_raw_bank=False)

def make_lines():
    lines = [all_lines['Hlt2BnoC_BdsToKSPi0_LL']()]
    # lines = [all_lines['Hlt2BnoC_BuToKpPi0']()]
    return lines
    
with reconstruction.bind(from_file=False):    
    config = run_moore(options, make_lines, public_tools)
