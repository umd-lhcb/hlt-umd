###############################################################################
# (c) Copyright 2020-2021 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
import logging
from functools import partial
from Functors import (
    ALL, BPVCORRM, BPVDIRA, BPVETA, BPVFDCHI2, BPVIPCHI2, CHARGE, CHI2DOF,
    CLOSESTTOBEAM, COMB, COV, DOCA, DOCACHI2, ETA, EVENTNUMBER, EVENTTYPE,
    FILTER, GHOSTPROB, ISMUON, P, PHI, POD, PT, MYPT, MASS, MAX, MAXDOCA,
    MAXDOCACHI2, MAXDOCACUT, MAXDOCACHI2CUT, MIN, MINIP, MINIPCHI2,
    MINIPCHI2CUT, MINIPCUT, MTDOCACHI2, MVA, NDOF, NHITS, NONE, PID_E, PID_K, PID_MU,
    PID_P, PID_PI, PROBNN_D, PROBNN_E, PROBNN_GHOST, PROBNN_K, PROBNN_MU,
    PROBNN_P, PROBNN_PI, QOVERP, RUNNUMBER, SIZE, SUM, TX, TY, X, Y, Z)
from Functors.grammar import BoundFunctor
from PyConf.Algorithms import Gaudi__Examples__IntDataProducer
"""
This file defines various sets of functors based on what class of physics
object they should work on. This should be a useful resource for test cases.
"""

logger = logging.getLogger(__name__)


def do_not_execute(obj):
    obj._do_not_execute = True
    return obj


# Functors don't check the C++ types of their data dependencies, so we can use
# the output of any algorithm as a dummy placeholder
DUMMY_DATA_DEP = Gaudi__Examples__IntDataProducer().OutputLocation
DUMMY_DATA_DEP_WARNING_TRACKER = set()


def simple_wrap(FUN, input_type, path_name, type_name=None,
                expects_list=False):
    if type_name is None:
        type_name = path_name + 'Type'

    def helper(inputs={}):
        df = inputs.get(input_type, None)
        if df is None:
            if not logger.disabled:
                log_message = "Using a dummy value for property '{}' of category '{}'".format(
                    path_name, input_type)
                if log_message not in DUMMY_DATA_DEP_WARNING_TRACKER:
                    logger.info(log_message)
                    DUMMY_DATA_DEP_WARNING_TRACKER.add(log_message)
            return do_not_execute(
                FUN(
                    **{
                        path_name:
                        [DUMMY_DATA_DEP] if expects_list else DUMMY_DATA_DEP
                    }))
        else:

            def handle_single_dh(dh):
                # This is clearly horrible and must be improved!
                header = {
                    'LHCb::Pr::Fitted::Forward::Tracks':
                    'Event/PrFittedForwardTracks.h'
                }[dh.type]
                return dh, (dh.type, header)

            # if 'df' is iterable, try handling it as a list of datahandles
            try:
                path_value, type_value = [], []
                for dh in df:
                    path, type_str = handle_single_dh(dh)
                    path_value.append(path)
                    type_value.append(type_str)
            except TypeError:
                # try handling it as a single datahandle
                path_value, type_value = handle_single_dh(df)

            return FUN(**{path_name: path_value, type_name: type_value})

    return helper


def test_mva_with_four_inputs(in1, in2, in3, in4):
    """
    This is a bit of a cheat, but we can abuse the Hlt1TwoTrackMVA as a general
    test of an MVA that takes four inputs.
    """
    return MVA(
        MVAType='MatrixNet',
        Config={'MatrixnetFile': "${PARAMFILESROOT}/data/Hlt1TwoTrackMVA.mx"},
        Inputs={
            'chi2': in1,
            'fdchi2': in2,
            'sumpt': in3,
            'nlt16': in4,
        })


DICT = {
    # These ones should be valid for all input types, including void
    'Generic': {
        'Functors': [
            ALL,
            NONE,
        ],
    },
    # These are void functors (i.e. they take no argument)
    'Void': {
        'Functors': [
            simple_wrap(RUNNUMBER, 'ODIN', 'ODIN'),
            simple_wrap(EVENTTYPE, 'ODIN', 'ODIN'),
            simple_wrap(EVENTNUMBER, 'ODIN', 'ODIN'),
            simple_wrap(SIZE, 'Container', 'Container'),
        ],
        'Includes': ['Generic'],
    },
    'Particle': {
        'Functors': [
            CHI2DOF,
            NDOF,
            PT,
            MYPT,
            ETA,
            P,
            PHI,
            CHARGE,
            simple_wrap(BPVIPCHI2, 'PVs', 'Vertices'),
            simple_wrap(MINIP, 'PVs', 'Vertices'),
            simple_wrap(MINIPCHI2, 'PVs', 'Vertices'),
            simple_wrap(partial(MINIPCUT, IPCut=7.0), 'PVs', 'Vertices'),
            simple_wrap(
                partial(MINIPCHI2CUT, IPChi2Cut=7.0), 'PVs', 'Vertices'),
            test_mva_with_four_inputs(PT, CHI2DOF, PT, CHI2DOF),
        ],
        'Includes': ['Generic'],
    },
    'Composite': {
        'Functors': [
            simple_wrap(
                partial(COMB, Functor=SUM(ETA)),
                'Children',
                'ChildContainers',
                type_name='ChildContainerTypes',
                expects_list=True),
            simple_wrap(BPVETA, 'PVs', 'Vertices'),
            simple_wrap(BPVDIRA, 'PVs', 'Vertices'),
            simple_wrap(BPVFDCHI2, 'PVs', 'Vertices'),
            simple_wrap(BPVCORRM, 'PVs', 'Vertices'),
        ],
        'Includes': ['Particle'],
    },
    'Composite2Body': {
        'Functors': [
            MASS,
        ],
        'Includes': ['Composite'],
    },
    # Functors that make sense for track objects
    'Track': {
        'Functors': [
            CLOSESTTOBEAM(X),
            CLOSESTTOBEAM(Y),
            CLOSESTTOBEAM(Z),
            CLOSESTTOBEAM(TX),
            CLOSESTTOBEAM(TY),
            CLOSESTTOBEAM(QOVERP),
            QOVERP,
            NHITS,
        ] + [
            CLOSESTTOBEAM(COV(Row=row, Col=col)) for row in range(5)
            for col in range(row, 5)
        ],
        'Includes': ['Particle'],
    },
    # Functors that make sense for ChargedBasics, i.e. tracks + PIDs
    'TrackWithMuonID': {
        'Functors': [
            ISMUON,
            # FIXME the following cannot be here
            # PID_MU, PID_PI, PID_K, PID_P, PID_E, PROBNN_D, PROBNN_E,
            # PROBNN_GHOST, PROBNN_K, PROBNN_MU, PROBNN_P,
        ],
        'Includes': ['Particle', 'Track'],
    },
    # Functors that make sense in a CombinationCut
    'Combination': {
        'Functors': [
            MAXDOCA,
            MAXDOCACHI2,
            DOCA(1, 2),
            DOCACHI2(1, 2),
            MTDOCACHI2(1, 2),
            MAXDOCACUT(10.),
            MAXDOCACHI2CUT(10.),
            SUM(PT),
            MIN(PT),
            MAX(PT),
        ],
        'Includes': ['Generic'],
    },

    # Eventually this list should be empty. These are functors that are not valid
    # for any of the current input types
    'Untestable': {
        'Functors': [
            FILTER(ALL),
            GHOSTPROB,
            POD(ALL),
        ],
    },
}


def functors_for_class(object_class,
                       inputs={},
                       exclusions=[],
                       cannot_execute=[],
                       skip_includes=False):
    """Return a set of functors that should work on objects of a particular
    class, defined in a rather loose/physics sense. e.g. Particle,
    ChargedBasic, Composite, Track, ...

    The optional dict of inputs contains DataHandles to various things the
    functors might depend on, e.g. PVs, ODIN, children...
    If a DataHandle is available for e.g. PVs, it will be passed to the
    relevant functor, otherwise a default, dummy value will be passed.

    The optional list of exclusions indicates a set of functors that should be
    removed from the list. This should be a list of strings

    The optional list `cannot_execute` indicates a set of functors that should
    not be included in the list of "executable" functors.

    The return value is divided into two subsets:
     - functors with no data dependencies, or where the data dependency was met
       in `inputs` -- these should be executable in tests
     - functors where a dummy value was used -- here only compilation, types
       and initialisation can be safely tested
    """
    assert object_class in DICT
    class_info = DICT[object_class]
    functors = class_info.get('Functors', [])
    # apply 'inputs' to 'functors'
    bound_can_execute, bound_cannot_execute = [], []
    for f in functors:
        if isinstance(f, BoundFunctor):
            bound_f = f
        else:
            bound_f = f(inputs=inputs)
            assert isinstance(bound_f, BoundFunctor)
        # split by whether it's executable
        if hasattr(bound_f, '_do_not_execute'):
            bound_cannot_execute.append(bound_f)
        else:
            bound_can_execute.append(bound_f)

    includes = [] if skip_includes else class_info.get('Includes', [])
    for included_class in includes:
        can_exec, cannot_exec = functors_for_class(
            included_class, inputs=inputs)
        bound_can_execute += can_exec
        bound_cannot_execute += cannot_exec

    # apply the list of exclusions
    exclude_names = set(exclusions)
    removed_names = set()

    def check(f):
        retain = f.name() not in exclude_names
        if not retain: removed_names.add(f.name())
        return retain

    for functor_list in [bound_can_execute, bound_cannot_execute]:
        functor_list[:] = [f for f in functor_list if check(f)]

    # shuffle functors from `bound_can_execute` to `bound_cannot_execute if needed`
    cannot_execute_names = set(cannot_execute)
    moved_names = set()

    def can_execute(f):
        move = f.name() in cannot_execute_names
        if move:
            moved_names.add(f.name())
            bound_cannot_execute.append(f)
        return not move

    bound_can_execute = [f for f in bound_can_execute if can_execute(f)]

    # check for exclusions that didn't match anything
    unused_exclusions = exclude_names - removed_names
    if len(unused_exclusions):
        logger.info(
            'There were redundant exclusions: {}'.format(unused_exclusions))

    # check for `cannot_execute` entries that didn't match anything
    unused_execute_overrides = cannot_execute_names - moved_names
    if len(unused_execute_overrides):
        logger.info('There were redundant "can execute" overrides: {}'.format(
            unused_execute_overrides))
    return bound_can_execute, bound_cannot_execute


# Do some token all-Python testing -- first collect all the functors that are
# listed in this file
ALL_CLASSES = DICT.keys()
logger.disabled = True
ALL_FUNCTORS = [
    f for class_name in ALL_CLASSES
    for functor_list in functors_for_class(class_name, skip_includes=True)
    for f in functor_list
]
logger.disabled = False

# Do some basic checks
# TODO: include a round-trip test of the pretty representations? perhaps in a separate unit test?
for func in ALL_FUNCTORS:
    assert len(func.code()) > 0 and len(func.headers()) > 0 and len(
        func.code_repr()) > 0 and len(func.name()) > 0
