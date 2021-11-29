## --------------------------------------------------------------------------------
## Lines for modes without a reconstructible decay vertex
## Ie, B+ -> K+ pi0 and B0 -> K0 pi0
## Author: Jason Andrews, jea@umd.edu, Emily Jiang, ejiang@umd.edu (run 3 corrections)
## --------------------------------------------------------------------------------
import Functors as F
from Hlt2Conf.algorithms_thor import ParticleFilter, ParticleCombiner, require_all
from Functors.math import in_range

# The GEC, what is GEC?
#class TrackGEC(Hlt2VoidFilter):
#    """Apply the GEC in number of tracks."""
#    def __init__(self):
#        from HltTracking.Hlt2TrackingConfigurations import \
#            Hlt2BiKalmanFittedForwardTracking as Hlt2LongTracking #hmm
#        tracks = Hlt2LongTracking().hlt2PrepareTracks()
#        code = (CONTAINS('%s') % tracks.outputSelection()) + < %(NTRACK_MAX)s
#        super(TrackGEC, self).__init__('B2Kpi0TrackGEC',
#                                       code,
#                                       [tracks],
#                                       nickname='TrackGEC',
#                                       shared=True)

#inputs make_kaons --> make_long_cb_kaons() from Hlt2Conf.standard_particles        
def KaonFilter(particles, pvs,
               TRACK_TRCHI2DOF_MAX,
               TRACK_TRGHOSTPROB_MAX,
               TRACK_PT_MIN,
               TRACK_P_MIN,
               TRACK_IPCHI2_MIN,
               TRACK_PIDK_MIN):
    """Filter the Hlt1 TOS track"""
    cut = ((F.CHI2DOF < TRACK_TRCHI2DOF_MAX)
           & (F.GHOSTPROB < TRACK_TRGHOSTPROB_MAX)
           & (F.PT > TRACK_PT_MIN)
           & (F.P > TRACK_P_MIN)
           & (F.MINIPCHI2(pvs) > TRACK_IPCHI2_MIN)
           & (F.PID_K > TRACK_PIDK_MIN))
    return ParticleFilter(particles, F.FILTER(cut))
    
def KS0Filter(particles, pvs,
              KS0_PT_MIN, 
              KS0_P_MIN,
              KS0_MASS,
              KS0_ADMASS,
              KS0_VCHI2PDOF_MAX,
              KS0_IPCHI2_MIN,
              TRACK_TRGHOSTPROB_MAX,
              TRACK_TRCHI2DOF_MAX):
            
    """Filter the Hlt1 TOS K0"""
    cut = ((F.PT > KS0_PT_MIN) 
           & (F.P > KS0_P_MIN)
           & in_range((KS0_MASS-KS0_ADMASS), F.MASS, (KS0_MASS+KS0_ADMASS))
           #& (abs((KS0_MASS-F.MASS())/MeV) < KS0_ADMASS) #ADMASS (use DMASS?)
           #& (CHILDCUT((F.GHOSTPROB < %(TRACK_TRGHOSTPROB_MAX)s),1))#CHILDCUT
           #& (CHILDCUT((F.GHOSTPROB < %(TRACK_TRGHOSTPROB_MAX)s),2))
           #& (CHILDCUT((F.CHI2DOF < %(TRACK_TRCHI2DOF_MAX)s),1))
           #& (CHILDCUT((F.CHI2DOF < %(TRACK_TRCHI2DOF_MAX)s),2))
           & ((F.CHI2/F.NDOF) < KS0_VCHI2PDOF_MAX) #VFASPF removed, VCHI2DOF-->F.CHI2/F.NDOF
           & (F.MINIPCHI2CUT(IPChi2Cut=KS0_IPCHI2_MIN, Vertices=pvs)))
           #& ((F.MINIPCHI2(pvs) - F.MINIPCHI2CUT(pvs)) > KS0_IPCHI2_MIN))
    return ParticleFilter(particles, F.FILTER(cut))

#come back here later
#MergedPi0s from inputs --> make_merged_pi0s from Hlt2Conf.standard_particles
                
def Pi0Filter(particles, pvs,
              PI0_PT_MIN,
              PI0_P_MIN):
    """Filter pi0s"""
    cut = (F.PT > PI0_PT_MIN) & (F.P > PI0_P_MIN)
    return ParticleFilter(particles, F.FILTER(cut)) 
