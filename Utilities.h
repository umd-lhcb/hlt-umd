/*****************************************************************************\
 * (c) Copyright 2020 CERN for the benefit of the LHCb Collaboration           *
 *                                                                             *
 * This software is distributed under the terms of the GNU General Public      *
 * Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
 *                                                                             *
 * In applying this licence, CERN does not waive the privileges and immunities *
 * granted to it by virtue of its status as an Intergovernmental Organization  *
 * or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once
#include "SelKernel/State4.h"
#include "SelKernel/Utilities.h"

#include "LHCbMath/MatVec.h"

namespace Sel {
  /** @fn    stateVectorFromComposite
   *  @brief Convert a particle (4-momentum + position) to 4-state (at some z).
   *
   *  Adapted from ParticleVertexFitter.cpp / Wouter Huslbergen, templated for
   *  new event model types and explicit vectorisation.
   */
  template <typename particle_t>
    auto stateVectorFromComposite( particle_t const& p ) {
    using Sel::Utils::covMatrix;
    using Sel::Utils::deref_if_ptr;
    using Sel::Utils::posCovMatrix;
    using Sel::Utils::threeMomCovMatrix;
    using Sel::Utils::threeMomPosCovMatrix;
    auto const invpz = 1.f / deref_if_ptr( p ).pz();
    using float_v    = std::decay_t<decltype( invpz )>;

    State4<float_v> s{};
    s.x()  = deref_if_ptr( p ).x();
    s.y()  = deref_if_ptr( p ).y();
    s.z()  = deref_if_ptr( p ).z();
    s.tx() = deref_if_ptr( p ).px() * invpz;
    s.ty() = deref_if_ptr( p ).py() * invpz;

    auto Jxpos    = LHCb::LinAlg::initialize_with_zeros<LHCb::LinAlg::Mat<float_v, 2, 3>>();
    Jxpos( 0, 0 ) = 1.f;
    Jxpos( 1, 1 ) = 1.f;
    Jxpos( 0, 2 ) = -s.tx();
    Jxpos( 1, 2 ) = -s.ty();

    auto Jtxmom    = LHCb::LinAlg::initialize_with_zeros<LHCb::LinAlg::Mat<float_v, 2, 3>>();
    Jtxmom( 0, 0 ) = invpz;
    Jtxmom( 1, 1 ) = invpz;
    Jtxmom( 0, 2 ) = -s.tx() * invpz;
    Jtxmom( 1, 2 ) = -s.ty() * invpz;

    s.covXX() = similarity( Jxpos, posCovMatrix( deref_if_ptr( p ) ) );
    s.covTT() = similarity( Jtxmom, threeMomCovMatrix( deref_if_ptr( p ) ) );
    s.covXT() = Jxpos * threeMomPosCovMatrix( deref_if_ptr( p ) ).transpose() * Jtxmom.transpose();

    

    std::cout<<"posSlopeCov from utils = "<<s.posSlopeCovariance()<<std::endl;
    /*
    //Commenting for now. This ought to moved into a test somewhere, and not calculated on the fly every time one wants to get a state-like vector from a composite.
    // For the computation of the Jacobian it is important to understand the following.
    //  x' = x + (z' - z) * tx
    // --> dx/dz = - tx ;
    //
    // Notation:
    //      J_{4,6} = ( Jxpos      Jxmom )
    //                ( Jtxpos     Jtxmom )
    // Jtxpos = 0 and Jxmom=0, so we rather do not introduce them. Instead, we'll compute Jxpos and Jxtxmom
    //
    // to test this, the full calculation is included below
    auto Jxpos    = LHCb::LinAlg::initialize_with_zeros<LHCb::LinAlg::Mat<float_v, 2, 3>>();
    Jxpos( 0, 0 ) = 1.f;
    Jxpos( 1, 1 ) = 1.f;
    Jxpos( 0, 2 ) = -s.tx();
    Jxpos( 1, 2 ) = -s.ty();

    auto Jtxmom    = LHCb::LinAlg::initialize_with_zeros<LHCb::LinAlg::Mat<float_v, 2, 3>>();
    Jtxmom( 0, 0 ) = invpz;
    Jtxmom( 1, 1 ) = invpz;
    Jtxmom( 0, 2 ) = -s.tx() * invpz;
    Jtxmom( 1, 2 ) = -s.ty() * invpz;

    if constexpr ( false ) {
    // Do the full calculation
    auto J    = LHCb::LinAlg::initialize_with_zeros<LHCb::LinAlg::Mat<float_v, 4, 6>>();
    J( 0, 0 ) = 1.f;
    J( 1, 1 ) = 1.f;
    J( 0, 2 ) = -s.tx();
    J( 1, 2 ) = -s.ty();
    J( 2, 3 ) = invpz;
    J( 3, 4 ) = invpz;
    J( 2, 5 ) = -s.tx() * invpz;
    J( 3, 5 ) = -s.ty() * invpz;
    // Get the {x, y, z, px, py, pz, pe} covariance matrix
    auto const cov77 = covMatrix( deref_if_ptr( p ) );
    // Strip off the pe part
    auto const cov66 = cov77.template sub<LHCb::LinAlg::MatSym<float_v, 6>, 0, 0>();
    // Transform to {x, y, tx, ty}
    auto const cov44 = similarity( J, cov66 );
    // Assemble XX, XT and TT into a 4x4 matrix
    auto const cov44_alt = s.posSlopeCovariance();
    // Check if that matches the result of the full calculation
    auto const cov44_del = cov44 - cov44_alt;
    // std::cout << "cov44_del = " << cov44_del << std::endl;
      }
    */
    return s;
  }
} // namespace Sel
