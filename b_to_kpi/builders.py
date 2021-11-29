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

from GaudiKernel.SystemOfUnits import GeV, MeV, mm, picosecond
from Hlt2Conf.algorithms import require_all, ParticleCombinerWithPVs, N3BodyCombinerWithPVs, N4BodyCombinerWithPVs
from PyConf import configurable
from Hlt2Conf.algorithms_thor import ParticleCombiner

import Functors as F
from Functors.math import in_range

 
@configurable
#Build the b hadron from hadrons and photons

def make_hb2x(kaons, pions, pvs,
        MTDOCACHI2_MAX,
        MASS_MIN,
        MASS_MAX,
        ASUM_PT_MIN,
        PT_MIN,
        name="Hb2XNeutralCombiner",
        functor_backend='ThOr',
        combiner='Th0rCombiner',
        alg_backend='Best',
        decay_descriptor ="B+ -> K+ pi0"):        
    combination_code = (in_range(MASS_MIN *MeV, F.MASS, MASS_MAX *MeV)) & (F.SUM(F.PT)> ASUM_PT_MIN) #not sure about ASUM(PT)-->F.SUM(F.PT)
    mother_code = (F.PT > PT_MIN)# & (MTDOCACHI2(1)<%(MTDOCACHI2_MAX)s)
    return ParticleCombiner(
        [kaons, pions],
        DecayDescriptor='[B+ -> K+ pi0]cc',
        CombinationCut=combination_code,
        CompositeCut=mother_code,
    )
