## --------------------------------------------------------------------------------
## Lines for modes without a reconstructible decay vertex
## Ie, B+ -> K+ pi0 and B0 -> K0 pi0
## Emily Jiang, ejiang@umd.edu
## --------------------------------------------------------------------------------

from Moore.config import register_line_builder
from Moore.lines import Hlt2Line
from PyConf import configurable
from GaudiKernel.SystemOfUnits import MeV, mm

#builders
from RecoConf.reconstruction_objects import make_pvs, upfront_reconstruction
from Hlt2Conf.standard_particles import make_has_rich_long_kaons, make_merged_pi0s, make_KsLL
from Hlt2Conf.lines.b_to_kpi.builders import make_hb2x

#filters
from Hlt2Conf.lines.b_to_kpi.filters import Pi0Filter, KaonFilter, KS0Filter

## Dict for each stage, plus 'Common'
cuts = {'Common'       : {'NTRACK_MAX':300,
                          'TRACK_TRCHI2DOF_MAX'  :3.0,
                          'TRACK_TRGHOSTPROB_MAX':0.5,
                          'L0FILTER'  :"L0_CHANNEL('Photon')|L0_CHANNEL('Electron')",
                          'PI0_TISTOS':'L0(Photon|Electron).*Decision%TOS'},

        'FilteredPi0s' : {'PI0_PT_MIN':3500 * MeV,
                          'PI0_P_MIN' :5000 * MeV },

        'FilteredKaons': {'TRACK_TISTOS'    :'Hlt1TrackMVA.*Decision%TOS',
                          'TRACK_PT_MIN'    :1200 * MeV,
                          'TRACK_P_MIN'     :12000 * MeV,
                          'TRACK_IPCHI2_MIN':50,
                          'TRACK_PIDK_MIN'  :-0.5 },

        'FilteredKSs'  : {'KS0_TISTOS':None,
                          'KS0_PT_MIN':500 * MeV,
                          'KS0_P_MIN' :8000 * MeV,
                          'KS0_ADMASS':15 * MeV,
                          'KS0_MASS':498 * MeV, #added because we don't have ADMASS in ThOr
                          'KS0_VCHI2PDOF_MAX':15,
                          'KS0_IPCHI2_MIN'   :10.0},

        'B2Kpi0'       : {'HLT1FILTER' :"HLT_PASS_RE('Hlt1TrackMVA.*Decision')",
                          'MASS_MIN'   :4000, ## units (MeV) are in cut string
                          'MASS_MAX'   :6200, ## units (MeV) are in cut string
                          'ASUM_PT_MIN':6500 * MeV,
                          'PT_MIN'     :5000 * MeV,
                          'MTDOCACHI2_MAX':8.0 },

        'B2K0pi0'      : {'HLT1FILTER' :"HLT_PASS_RE('Hlt1(Two)?TrackMVA.*Decision')",
                          'MASS_MIN'   :4000, ## units (MeV) are in cut string
                          'MASS_MAX'   :6200, ## units (MeV) are in cut string
                          'ASUM_PT_MIN':5000 * MeV,
                          'PT_MIN'     :4000 * MeV,
                          'MTDOCACHI2_MAX':10.0 },
        }

all_lines={}
@register_line_builder(all_lines)
@configurable

def BToKpi0_line(name='Hlt2BToKpi0_Line', prescale=1):
    functor_backend = 'ThOr'
    pvs= make_pvs()
    pi0s = make_merged_pi0s()
    kaons = make_has_rich_long_kaons()
    ksll = make_KsLL()

    filteredPi0s = Pi0Filter(pi0s, pvs, 
                             cuts['FilteredPi0s']['PI0_PT_MIN'],
                             cuts['FilteredPi0s']['PI0_P_MIN'])
    filteredKaons = KaonFilter(kaons, pvs,
                               cuts['Common']['TRACK_TRCHI2DOF_MAX'],
                               cuts['Common']['TRACK_TRGHOSTPROB_MAX'],
                               cuts['FilteredKaons']['TRACK_PT_MIN'],
                               cuts['FilteredKaons']['TRACK_P_MIN'],
                               cuts['FilteredKaons']['TRACK_IPCHI2_MIN'],
                               cuts['FilteredKaons']['TRACK_PIDK_MIN'])
    filteredKSs = KS0Filter(ksll, pvs,
                            cuts['FilteredKSs']['KS0_PT_MIN'],
                            cuts['FilteredKSs']['KS0_P_MIN'],
                            cuts['FilteredKSs']['KS0_MASS'],
                            cuts['FilteredKSs']['KS0_ADMASS'],
                            cuts['FilteredKSs']['KS0_VCHI2PDOF_MAX'],
                            cuts['FilteredKSs']['KS0_IPCHI2_MIN'],
                            cuts['Common']['TRACK_TRGHOSTPROB_MAX'],
                            cuts['Common']['TRACK_TRCHI2DOF_MAX'])
    
    B2Kpi0 = make_hb2x(kaons = filteredKaons, pions = filteredPi0s, pvs = pvs,
                       MTDOCACHI2_MAX = cuts['B2Kpi0']['MTDOCACHI2_MAX'],
                       MASS_MIN = cuts['B2Kpi0']['MASS_MIN'],
                       MASS_MAX = cuts['B2Kpi0']['MASS_MAX'],
                       ASUM_PT_MIN = cuts['B2Kpi0']['ASUM_PT_MIN'],
                       PT_MIN = cuts['B2Kpi0']['PT_MIN'])
    
    return Hlt2Line(name=name, algs = (upfront_reconstruction()+[B2Kpi0]))
    
    
