###############################################################################
# (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
import Functors as F
from Functors.math import in_range
from GaudiKernel.SystemOfUnits import MeV, mm

from Moore.config import register_line_builder
from Moore.lines import Hlt2Line

from RecoConf.reconstruction_objects import (
    make_pvs_v2 as make_pvs,
    upfront_reconstruction,
)
from Hlt2Conf.standard_particles import (
    make_has_rich_long_kaons,
    make_merged_pi0s,
    make_KsLL,
)

from Hlt2Conf.algorithms_thor import ParticleCombiner, ParticleFilter, require_all

all_lines = {}

def filter_kaons(particles, 
                 pvs, 
                 trchi2dof_max=3.0, 
                 trghostprob_max=0.5,
                 pt_min=1200 * MeV, 
                 p_min=12000 * MeV,
                 mipchi2_min=50, 
                 pid_k_min=-0.5):
    cut = require_all(
        F.CHI2DOF < trchi2dof_max,
        F.GHOSTPROB < trghostprob_max,
        F.PT > pt_min,
        F.P > p_min,
        F.MINIPCHI2(pvs) > mipchi2_min,
        F.PID_K > pid_k_min)

    return ParticleFilter(particles, F.FILTER(cut))


def filter_pions(particles, 
                 pvs, 
                 pt_min=3500 * MeV, 
                 p_min=5000 * MeV):
    cut = require_all(
        F.PT > pt_min,
        F.P > p_min,
    )
    return ParticleFilter(particles, F.FILTER(cut))

def make_bs(kaons,
            pions,
            pvs,
            two_body_comb_maxdocachi2=9.0,
            comb_m_min=4000 * MeV,
            comb_m_max=6200 * MeV,
            comb_pt_min=5000 * MeV,
            mtdocachi2_max=10.0,
            pt_min=4000 * MeV):
    two_body_combination_code = F.MAXDOCACHI2CUT(two_body_comb_maxdocachi2)
    combination_code = require_all(
        in_range(comb_m_min, F.MASS, comb_m_max),
        F.SUM(F.PT) > comb_pt_min,
    )
    composite_code = require_all(
        F.MTDOCACHI2(pvs)<mtdocachi2_max,
        #F.MTDOCACHI2(1)<mtdocachi2_max,
        #F.BPVVDZ(pvs, 1)<mtdocachi2_max,
        #F.BPVVDZ(pvs)<mtdocachi2_max,
        F.PT > pt_min,
    )
    return ParticleCombiner(
        [kaons, pions],
        ParticleCombiner="ParticleAdder",
        DecayDescriptor="[B+ -> K+ pi0]cc",
        CombinationCut=combination_code,
        CompositeCut=composite_code,
    )

@register_line_builder(all_lines)
def BToKpi0_line(name="Hlt2BToKpi0_Line", prescale=1):
    pvs = make_pvs()
    kaons = filter_kaons(make_has_rich_long_kaons(), pvs)
    pions = filter_pions(make_merged_pi0s(), pvs)
    bs = make_bs(kaons, pions, pvs)

    return Hlt2Line(
        name=name,
        algs=upfront_reconstruction() + [bs],
        prescale=prescale,
    )
