#Moore configurations
from Moore import options, run_moore
#from Hlt2Conf.lines.lambdab_to_lambdacpi import lbtolcpi_lctopkpi_line
from Hlt2Conf.lines.b_to_kpi import BToKpi0_line
from RecoConf.global_tools import stateProvider_with_simplified_geom
from RecoConf.reconstruction_objects import upfront_reconstruction, reconstruction
from RecoConf.hlt1_tracking import default_ft_decoding_version

from Gaudi.Configuration import FileCatalog
from GaudiConf import IOHelper

#set decoding version
ft_decoding_version = 6
default_ft_decoding_version.global_bind(value=ft_decoding_version)

public_tools = [stateProvider_with_simplified_geom()]

input_files = ['LFN:/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000056_1.xdigi',]
#input_files = ['root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00127044/0000/00127044_00000056_1.xdigi']
options.input_files = input_files
options.input_type = "ROOT"
options.input_raw_format = 0.3

FileCatalog().Catalogs = ["xmlcatalog_file:/afs/cern.ch/work/e/ejiang/public/stack/myCatalogDown.xml"]
    
#options.evt_max = 100                                                                       
options.dddb_tag = 'dddb-20201211' 
options.conddb_tag = 'sim-20201218-vc-md100' 

#options.control_flow_file = 'control_flow.gv'

def all_lines():
    #return [lbtolcpi_lctopkpi_line()]
    return [BToKpi0_line()]
    
with reconstruction.bind(from_file=False):
    run_moore(options, all_lines, public_tools)
        
    
    
