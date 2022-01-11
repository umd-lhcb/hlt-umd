###############################################################################
# (c) Copyright 2019-2021 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
from Functors.grammar import Functor, BoundFunctor
from PyConf.dataflow import DataHandle
FILTER = Functor(
    'FILTER',
    'Filter',
    'Filter.h',
    'Adapt a predicate to filter a container.',
    Params=[('Functor', 'Predicate to filter the container with.',
             BoundFunctor)])
X = Functor('X', "Track::X", "TrackLike.h", 'X coordinate.')
Y = Functor('Y', "Track::Y", "TrackLike.h", 'Y coordinate.')
Z = Functor('Z', "Track::Z", "TrackLike.h", 'Z coordinate.')
TX = Functor('TX', "Track::TX", "TrackLike.h", 'X slope.')
TY = Functor('TY', "Track::TY", "TrackLike.h", 'Y slope.')
COV = Functor(
    'COV',
    'Track::Covariance',
    'TrackLike.h',
    'Get the specified (i, j)th element of the covariance matrix.',
    Params=[('Row', 'Row to access', int), ('Col', 'Column to access', int)])
P = Functor('P', "Track::Momentum", "TrackLike.h", "Momentum.")
PT = Functor('PT', "Track::TransverseMomentum", "TrackLike.h",
             "Transverse momentum.")
MYPT = Functor ('MYPT', "Track::MyTransverseMomentum", "TrackLike.h", "My transverse Momentum.")
PHI = Functor('PHI', "Track::Phi", "TrackLike.h", "Azimuthal angle, phi.")
ETA = Functor('ETA', "Track::PseudoRapidity", "TrackLike.h", "Pseudorapidity.")
ISMUON = Functor('ISMUON', "Track::IsMuon", "TrackLike.h", "IsMuon.")
NDOF = Functor('NDOF', 'Track::nDoF', 'TrackLike.h',
               'Number of degrees of freedom [for chi2]')
QOVERP = Functor('QOVERP', 'Track::QoverP', 'TrackLike.h', 'q/p')
CHI2DOF = Functor(
    'CHI2DOF', "Track::Chi2PerDoF", "TrackLike.h",
    "Vertex or track chi^2 per degree of freedom (works for both).")
CHI2 = Functor('CHI2', "Track::Chi2", "TrackLike.h",
               "Vertex or track chi^2 (works for both).")
GHOSTPROB = Functor('GHOSTPROB', "Track::GhostProbability", "TrackLike.h",
                    "Ghost probability.")
CLOSESTTOBEAM = Functor(
    'CLOSESTTOBEAM',
    'Track::ClosestToBeamState',
    'TrackLike.h',
    'Adapter to apply a functor to the closest-to-beam state of the argument.',
    Params=[('Functor', 'Functor to apply to the closest to beam state',
             BoundFunctor)])
NHITS = Functor('NHITS', "Track::nHits", "TrackLike.h",
                "Track number of hits.")
MINIP = Functor(
    'MINIP',
    "Track::MinimumImpactParameter",
    "TrackLike.h",
    """Calculate the minimum impact parameter w.r.t. any of the given vertices.
    MINIPCUT may be more efficient.""",
    Params=[('Vertices', 'TES location of input [primary] vertices',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
MINIPCHI2 = Functor(
    'MINIPCHI2',
    "Track::MinimumImpactParameterChi2",
    "TrackLike.h",
    """Calculate the minimum impact parameter chi2 w.r.t. any of the given
    vertices. MINIPCHI2CUT may be more efficient.""",
    Params=[('Vertices', 'TES location of input [primary] vertices',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
MINIPCUT = Functor(
    'MINIPCUT',
    "Track::MinimumImpactParameterCut",
    "TrackLike.h",
    """Require the minimum impact parameter w.r.t. any of the given vertices is
    greater than some threshold.""",
    Params=[('IPCut', 'The impact parameter cut value', float),
            ('Vertices', 'TES location of input [primary] vertices',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
MINIPCHI2CUT = Functor(
    'MINIPCHI2CUT',
    "Track::MinimumImpactParameterChi2Cut",
    "TrackLike.h",
    """Require the minimum impact parameter chi2 w.r.t. any of the given
    vertices is greater than some threshold.""",
    Params=[('IPChi2Cut', 'The impact parameter chi2 cut value', float),
            ('Vertices', 'TES location of input [primary] vertices',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
ALL = Functor('ALL', "AcceptAll", "Function.h",
              "Accept everything; always evaluates to 'true'.")
NONE = Functor('NONE', "AcceptNone", "Function.h",
               "Accept nothing; always evaluates to 'false'.")
SUM = Functor(
    'SUM',
    "Adapters::Accumulate",
    "Adapters.h",
    "Calculate the [scalar] sum of the given functor value.",
    Params=[('Functor', 'The functor to accumulate the return value of.',
             BoundFunctor)])
MIN = Functor(
    'MIN',
    "Adapters::Minimum",
    "Adapters.h",
    "Calculate the minimum of the given functor value.",
    Params=[('Functor', 'The functor to find the minimum return value of.',
             BoundFunctor)])
MAX = Functor(
    'MAX',
    "Adapters::Maximum",
    "Adapters.h",
    "Calculate the maximum of the given functor value.",
    Params=[('Functor', 'The functor to find the maximum return value of.',
             BoundFunctor)])
CHILD = Functor(
    'CHILD',
    "Adapters::Child",
    "Adapters.h",
    "Apply functor on a child.",
    Params=[
        ('Index',
         'The index of the child to apply the functor to (starting from 1).',
         int),
        ('Functor', 'The functor to apply on the child.', BoundFunctor),
    ],
    AllowMultiplePositionalArguments=True)


def idxs_formatter(idxs):
    return ", ".join([str(i) for i in idxs]), []


SUBCOMB = Functor(
    'SUBCOMB',
    "Adapters::SubCombination",
    "Adapters.h",
    "Apply functor on a SubCombination.",
    Params=[
        ('Functor', 'The functor to apply on the subcombination.',
         BoundFunctor),
    ],
    TemplateParams=[
        ('Indices',
         'List of indices to build the subcombination from (starting from 1).',
         idxs_formatter)
    ])

MASSWITHHYPOTHESES = Functor(
    'MASSWITHHYPOTHESES',
    'Composite::MassWithHypotheses',
    'Composite.h',
    '''Invariant mass of a combined particle given child mass hypotheses.''',
    Params=[('Masses', 'Masses of the children', list)])
MASS = Functor('MASS', 'Composite::Mass', 'Composite.h',
               'Get the particle (composite or basic) mass.')
VX = Functor('VX', 'Composite::VertexX', 'Composite.h',
             'Get the particle vertex X. NOT using PV as reference.')
VY = Functor('VY', 'Composite::VertexY', 'Composite.h',
             'Get the particle vertex Y. NOT using PV as reference.')
VZ = Functor('VZ', 'Composite::VertexZ', 'Composite.h',
             'Get the particle vertex Z. NOT using PV as reference.')
VRho = Functor(
    'VRho', 'Composite::VertexRho', 'Composite.h',
    'Get the particle vertex Rho=sqrt(X*X+Y*Y). NOT using PV as reference.')
DOCA = Functor(
    'DOCA',
    "Combination::DistanceOfClosestApproach",
    "Combination.h",
    """Compute the distance of closest approach between two 'states'.""",
    Params=[('Child1',
             'Index [starting from 1] of the first child to consider.', int),
            ('Child2',
             'Index [starting from 1] of the second child to consider.', int)],
    TemplateParams=[('DistanceCalculator',
                     'Distance calculator implementation to use.')],
    AllowMultiplePositionalArguments=True)
DOCACHI2 = Functor(
    'DOCACHI2',
    "Combination::DistanceOfClosestApproachChi2",
    "Combination.h",
    """Compute the significance of the distance of closest
    approach between two 'states'.""",
    Params=[('Child1',
             'Index [starting from 1] of the first child to consider.', int),
            ('Child2',
             'Index [starting from 1] of the second child to consider.', int)],
    TemplateParams=[('DistanceCalculator',
                     'Distance calculator implementation to use.')],
    AllowMultiplePositionalArguments=True)
MTDOCACHI2 = Functor(
    'MTDOCACHI2',
    "Combination::MotherTrajectoryDistanceOfClosestApproachChi2",
    "Combination.h",
    """Compute the significance of the distance of closest
    approach between mother and child.""",
    Params=[('Child1',
             'Index [starting from 1] of the first child to consider.', int),
            ('Child2',
             'Index [starting from 1] of the second child to consider.', int)],
    TemplateParams=[('DistanceCalculator',
                     'Distance calculator implementation to use.')],
    AllowMultiplePositionalArguments=True)
ALV = Functor(
    'ALV',
    "Combination::CosAngleBetweenDecayProducts",
    "Combination.h",
    """Compute the cosine value of angle between two decay products.""",
    Params=[('Child1',
             'Index [starting from 1] of the first child to consider.', int),
            ('Child2',
             'Index [starting from 1] of the second child to consider.', int)],
    AllowMultiplePositionalArguments=True)
MAXDOCA = Functor(
    'DOCA', "Combination::MaxDistanceOfClosestApproach", "Combination.h",
    "Compute the maximum pairwise distance of closest approach between members of a combination."
)
MAXDOCACHI2 = Functor(
    'MAXDOCACHI2', "Combination::MaxDistanceOfClosestApproachChi2",
    "Combination.h",
    "Compute the maximum pairwise significance of the distance of closest approach between members of a combination."
)
MAXDOCACUT = Functor(
    'MAXDOCACUT',
    "Combination::MaxDistanceOfClosestApproachCut",
    "Combination.h",
    """Cut on the the distance of closest approach between two 'states'.""",
    Params=[('thresh', 'Threshold for cut', float)])
MAXDOCACHI2CUT = Functor(
    'MAXDOCACHI2CUT',
    "Combination::MaxDistanceOfClosestApproachChi2Cut",
    "Combination.h",
    """Cut on the significance of the distance of closest
    approach between two 'states'.""",
    Params=[('thresh', 'Threshold for cut', float)])
CHARGE = Functor("CHARGE", "Combination::Charge", "Combination.h",
                 "Compute the charge")
PID_MU = Functor('PID_MU', "Track::PIDmu", "TrackLike.h", "CombDLLmu.")
PID_PI = Functor('PID_PI', "Track::PIDpi", "TrackLike.h", "CombDLLpi.")
PID_K = Functor('PID_K', "Track::PIDk", "TrackLike.h", "CombDLLk.")
PID_P = Functor('PID_P', "Track::PIDp", "TrackLike.h", "CombDLLp.")
PID_E = Functor('PID_E', "Track::PIDe", "TrackLike.h", "CombDLLe.")
PROBNN_D = Functor('PROBNN_D', "Track::PROBNN_D_t", "TrackLike.h", "PROBNN_D.")
PROBNN_E = Functor('PROBNN_E', "Track::PROBNN_E_t", "TrackLike.h", "PROBNN_E.")
PROBNN_GHOST = Functor('PROBNNGHOST', "Track::PROBNN_GHOST_t", "TrackLike.h",
                       "PROBNN_GHOST.")
PROBNN_K = Functor('PROBNN_K', "Track::PROBNN_K_t", "TrackLike.h", "PROBNN_K.")
PROBNN_MU = Functor('PROBNN_MU', "Track::PROBNN_MU_t", "TrackLike.h",
                    "PROBNN_MU.")
PROBNN_P = Functor('PROBNN_P', "Track::PROBNN_P_t", "TrackLike.h", "PROBNN_P.")
PROBNN_PI = Functor('PROBNN_PI', "Track::PROBNN_PI_t", "TrackLike.h",
                    "PROBNN_PI.")
SIZE = Functor(
    'SIZE',
    "TES::Size",
    "TES.h",
    "Get the size of the given container",
    Params=[('Container',
             'The TES location of the container whose size we return',
             DataHandle)],
    TemplateParams=[('ContainerType', 'Type of the container')])
BPVETA = Functor(
    'BPVETA',
    'Composite::PseudoRapidityFromVertex',
    'Composite.h',
    '''Compute the pseudorapidity of the vector connecting the associated
    [primary] vertex to the composite particle decay vertex.''',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
BPVCORRM = Functor(
    'BPVCORRM',
    'Composite::CorrectedMass',
    'Composite.h',
    'Compute the corrected mass of the composite using the associated [primary] vertex.',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
BPVDIRA = Functor(
    'BPVDIRA',
    'Composite::CosDirectionAngleToVertex',
    'Composite.h',
    '''Compute the cosine of the angle between the particle momentum and the
    vector from the associated [primary] vertex to the composite particle
    decay vertex.''',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
BPVIPCHI2 = Functor(
    'BPVIPCHI2',
    'Track::ImpactParameterChi2ToVertex',
    'TrackLike.h',
    '''Return the impact parameter chi2 w.r.t. the associated vertex, assuming
    that it is from the container of vertices that is passed. If no association
    is available, compute one.''',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
BPVFDCHI2 = Functor(
    'BPVFDCHI2',
    'Composite::FlightDistanceChi2ToVertex',
    'Composite.h',
    '''Return the flight distance chi2 w.r.t. the associated vertex, assuming
    that it is from the container of vertices that is passed. If no association
    is available, compute one.''',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
BPVVDZ = Functor(
    'BPVVDZ',
    'Composite::DeltaZToVertex',
    'Composite.h',
    '''Return the z component of flight distance w.r.t. the associated vertex, assuming
    that it is from the container of vertices that is passed.
    ''',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
BPVVDRHO = Functor(
    'BPVVDRHO',
    'Composite::DeltaRhoToVertex',
    'Composite.h',
    '''Return the rho component of flight distance w.r.t. the associated vertex, assuming
    that it is from the container of vertices that is passed.
    ''',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[('VerticesType', 'Input vertex container type')])
BPVLTIME = Functor(
    'BPVLTIME',
    'Composite::ComputeLifetime',
    'Composite.h',
    '''Return the particle lifetime w.r.t. the associated vertex, assuming
    that it is from the container of vertices that is passed. If no association
    is available, compute one.''',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[
        ('VerticesType', 'Input vertex container type'),
    ])
BPVDLS = Functor(
    'BPVDLS',
    'Composite::ComputeDecayLengthSignificance',
    'Composite.h',
    '''Return the decay length significance w.r.t. the associated vertex, assuming
    that it is from the container of vertices that is passed. If no association
    is available, compute one.''',
    Params=[('Vertices', 'TES location of input [primary] vertices.',
             DataHandle)],
    TemplateParams=[
        ('VerticesType', 'Input vertex container type'),
    ])
RUNNUMBER = Functor(
    'RUNNUMBER',
    'TES::RunNumber',
    'TES.h',
    '''Return the run number from ODIN.''',
    Params=[('ODIN', 'TES location of ODIN information.', DataHandle)],
    TemplateParams=[('ODINType', 'Type of the ODIN object')])
EVENTNUMBER = Functor(
    'EVENTNUMBER',
    'TES::EventNumber',
    'TES.h',
    '''Return the event number from ODIN.''',
    Params=[('ODIN', 'TES location of ODIN information.', DataHandle)],
    TemplateParams=[('ODINType', 'Type of the ODIN object')])
EVENTTYPE = Functor(
    'EVENTTYPE',
    'TES::EventType',
    'TES.h',
    '''Return the event type from ODIN.''',
    Params=[('ODIN', 'TES location of ODIN information.', DataHandle)],
    TemplateParams=[('ODINType', 'Type of the ODIN object')])

# Examples:
Ex_TimesTwo = Functor('TimesTwo', 'Examples::TimesTwo', 'Example.h',
                      ''' calculate the double ''')

Ex_GreaterThan = Functor(
    'GreaterThan',
    'Examples::GreaterThan',
    'Example.h',
    ''' value greater than v? ''',
    Params=[('v', 'reference value', float)])

Ex_TBL = Functor('TBL', 'Examples::ThorBeatsLoki', 'Example.h', ''' ... ''')


def mva_input_formatter(config):
    from Functors.grammar import python_to_cpp_str
    mva_inputs = []
    headers = []
    # Data dependencies
    inputs = []
    for key, val in list(config.items()):
        assert type(key) == str
        key_code, key_headers = python_to_cpp_str(key)
        val_code, val_headers = python_to_cpp_str(val)
        mva_inputs.append('MVAInput( ' + key_code + ', ' + val_code + ' )')
        headers += key_headers + val_headers
        inputs += val.data_dependencies()
    return ', '.join(mva_inputs), list(set(headers)), inputs


def mva_impl_formatter(config):
    return 'Sel::' + config, ['SelTools/' + config + '.h']


MVA = Functor(
    'MVA',
    'MVA',
    'MVA.h',
    '''Evaluate an MVA.''',
    Params=[('Config', 'MVA config', dict),
            ('Inputs', 'MVA inputs', dict, mva_input_formatter)],
    TemplateParams=[('MVAType', 'The MVA implementation to use',
                     mva_impl_formatter)])


def comb_locations_formatter(locations):
    from Functors.grammar import python_to_cpp_str, value_matches_type
    expressions, headers, inputs = [], [], []
    for loc in locations:
        # The functor framework is unable to introspect type requirements of
        # composites, e.g. 'list-of-DataHandle', only 'list'. We enforce the
        # list-of-DataHandle type requirement of the COMB functor here instead
        assert value_matches_type(loc, DataHandle)
        code, headers = python_to_cpp_str(loc)
        expressions.append(code)
        headers += headers
        inputs.append(loc)
    return ', '.join(expressions), list(set(headers)), inputs


def comb_types_formatter(type_list):
    return ', '.join([x[0] for x in type_list]), list(
        set([header for x in type_list for header in x[1:]]))


COMB = Functor(
    'COMB',
    'Adapters::CombinationFromComposite',
    'Adapters.h',
    '''Re-constructs a 'combination' object from a 'composite' object and
    returns the result of applying the given functor to that 'combination'
    object.''',
    Params=[('Functor', "The functor to apply to the 'combination' object.",
             BoundFunctor),
            ('ChildContainers',
             "List of data locations where the child objects can be found.",
             list, comb_locations_formatter)],
    TemplateParams=[('ChildContainerTypes',
                     'List of types of the given child locations.',
                     comb_types_formatter)])

POD = Functor(
    'POD',
    'Adapters::ConvertToPOD',
    'Adapters.h',
    'Try to convert an object representing a scalar number into a plain C++ type. For example, convert SIMDWrapper::scalar::float_v to float.',
    Params=[('Functor', 'The functor to convert the return value of.',
             BoundFunctor)])