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
#include "SelKernel/VertexRelation.h"
#include "SelTools/Utilities.h"

#include "Event/Track_v3.h"
#include "LHCbMath/MatVec.h"
#include "LHCbMath/MatrixTransforms.h"
#include "LHCbMath/MatrixUtils.h"
#include "LHCbMath/SIMDWrapper.h"

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "GaudiKernel/Vector4DTypes.h"

namespace Sel {
  /** @class DistanceCalculator
   *  @brief Collection of functions for calculating quantities like IP,
   *         IPCHI2, DOCA, DOCACHI2, ...
   */
  using StateLocation = LHCb::Event::v3::Tracks::StateLocation;
  struct DistanceCalculator {
    DistanceCalculator( Gaudi::Algorithm* ) {}

  private:
    template <typename Particle>
    static auto getClosestToBeamState( Particle const& p ) {
      if constexpr ( Sel::Utils::has_tracklike_API<Particle> ) {
        return p.state( StateLocation::ClosestToBeam );
      } else if constexpr ( Sel::type_traits::has_closestToBeamState_v<Particle> ) {
        return p.closestToBeamState();
      } else {
        // Make a state vector from a particle (position + 4-momentum)
        return stateVectorFromComposite( p );
      }
    }

  public:
    /** @fn    stateDOCA
     *  @brief Calculate the distance of closest approach between two states.
     *
     *  This is copied from LHCb::TrackVertexUtils::doca so it can be templated
     *  on the state type, allowing the two states to be of different types. In
     *  the TrackVertexUtils version this would lead to a collision with the doca
     *  function taking a state and a point.
     */
    template <typename StateA, typename StateB>
    auto stateDOCA( StateA const& stateA, StateB const& stateB ) const {
      using std::sqrt; // allows sqrt() below to work with basic types + ADL
      // first compute the cross product of the directions.
      auto const txA = stateA.tx();
      auto const tyA = stateA.ty();
      auto const txB = stateB.tx();
      auto const tyB = stateB.ty();
      auto const nx  = tyA - tyB;             //   y1 * z2 - y2 * z1
      auto const ny  = txB - txA;             // - x1 * z2 + x2 * z1
      auto const nz  = txA * tyB - tyA * txB; //   x1 * y2 - x2 * y1
      auto const n   = sqrt( nx * nx + ny * ny + nz * nz );
      // compute the doca
      auto const dx    = stateA.x() - stateB.x();
      auto const dy    = stateA.y() - stateB.y();
      auto const dz    = stateA.z() - stateB.z();
      auto const ndoca = dx * nx + dy * ny + dz * nz;
      return ndoca / n;
    }

    /** @fn    stateDOCAChi2
     *  @brief Significance of DOCA between two states.
     */
    template <typename StateA, typename StateB>
    auto stateDOCAChi2( StateA const& sA, StateB const& sB ) const {
      // first compute the cross product of the directions. we'll need this in any case
      using float_v     = std::decay_t<decltype( sA.tx() )>;
      float_v const txA = sA.tx();
      float_v const tyA = sA.ty();
      float_v const txB = sB.tx();
      float_v const tyB = sB.ty();
      float_v const nx  = tyA - tyB;             //   y1 * z2 - y2 * z1
      float_v const ny  = txB - txA;             // - x1 * z2 + x2 * z1
      float_v const nz  = txA * tyB - tyA * txB; //   x1 * y2 - x2 * y1

      // compute doca. we don't divide by the normalization to save time. we call it
      // 'ndoca'
      float_v const dx    = sA.x() - sB.x();
      float_v const dy    = sA.y() - sB.y();
      float_v const dz    = sA.z() - sB.z();
      float_v const ndoca = dx * nx + dy * ny + dz * nz;

      // the hard part: compute the jacobians :-)
      using Vector4 = LHCb::LinAlg::Vec<float_v, 4>;
      Vector4 jacA, jacB;
      jacA( 0 ) = nx;
      jacA( 1 ) = ny;
      jacA( 2 ) = -dy + dz * tyB;
      jacA( 3 ) = dx - dz * txB;
      jacB( 0 ) = -nx;
      jacB( 1 ) = -ny;
      jacB( 2 ) = dy - dz * tyA;
      jacB( 3 ) = -dx + dz * txA;

      // compute the variance on ndoca
      float_v const varndoca =
          similarity( jacA, get::posSlopeCovariance( sA ) ) + similarity( jacB, get::posSlopeCovariance( sB ) );

      /*std::cout<<"varndoca = "<<varndoca<<std::endl;
      std::cout<<"simA = "<<similarity( jacA, get::posSlopeCovariance( sA ) )<<std::endl;
      std::cout<<"simB = "<<similarity( jacB, get::posSlopeCovariance( sB ) )<<std::endl;
      */
      // return the chi2
      return ndoca * ndoca / varndoca;
    }

    /** @fn    particleDOCA
     *  @brief Calculate the distance of closest approach between two particles.
     *  @todo  This needs to do something sensible when the particle is a
     *         downstream track.
     */
    template <typename ParticleA, typename ParticleB>
    auto particleDOCA( ParticleA const& pA, ParticleB const& pB ) const {
      return stateDOCA( getClosestToBeamState( pA ), getClosestToBeamState( pB ) );
    }

    /** @fn    particleDOCAChi2
     *  @brief Significance of DOCA between two particles.
     */
    template <typename ParticleA, typename ParticleB>
    auto particleDOCAChi2( ParticleA const& pA, ParticleB const& pB ) const {
      return stateDOCAChi2( getClosestToBeamState( pA ), getClosestToBeamState( pB ) );
    }
  };

  namespace helper {

    // FIXME: we should _not_ need this version...
    template <typename OUT, typename float_v, auto N, auto M>
    auto toSMatrix( LHCb::LinAlg::Mat<float_v, N, M> const& inMat ) {
      ROOT::Math::SMatrix<OUT, N, M> outMat;
      LHCb::Utils::unwind<0, N>(
          [&]( auto i ) { LHCb::Utils::unwind<0, M>( [&]( auto j ) { outMat( i, j ) = inMat( i, j ).cast(); } ); } );
      return outMat;
    }
    template <typename float_v, auto N, auto M>
    auto toSMatrix( LHCb::LinAlg::Mat<float_v, N, M> const& inMat ) {
      return toSMatrix<float_v, float_v, N, M>( inMat );
    }

  } // end of namespace helper

  struct LifetimeFitter {
    // FIXME: remove default argument...
    LifetimeFitter( Gaudi::Algorithm* owning_algorithm = nullptr ) : m_owning_algorithm{owning_algorithm} {}

  private:
    Gaudi::Algorithm*      m_owning_algorithm = nullptr;
    static constexpr float NonPhysicalValue   = std::numeric_limits<float>::quiet_NaN();

  public:
    /** @fn   DecayLengthSignificance
     *  @brief Calculate the significance of a non-zero decay length.
     *
     */

    template <typename Particle, typename VContainer>
    auto DecayLengthSignificance( VContainer const& primary, Particle const& particle ) const {
      using Sel::Utils::covMatrix;
      using Sel::Utils::endVertexPos;
      using Sel::Utils::posCovMatrix;
      using Sel::Utils::threeMomentum;
      using std::sqrt;

      // Calculate the distance between the particle and the vertex we hold.
      const auto flight = endVertexPos( particle ) - endVertexPos( primary );

      // Get the 3 momentum vectors for following calculations.
      const auto p3  = threeMomentum( particle );
      const auto dir = p3 / p3.mag();

      // Get the covariance matrix of the vertex.
      const auto vertexCov = posCovMatrix( primary );

      // Get the fulll particle covariance matrix.
      const auto pCov = covMatrix( particle );

      // Project the momentum of the particle onto its distance to the vertex
      const auto a = dot( dir, flight ) / p3.mag();

      // Update the covariance matrix
      std::decay_t<decltype( vertexCov )> W;
      for ( size_t row = 0; row < 3; ++row ) {
        for ( size_t col = 0; col <= row; ++col ) {
          W( row, col ) = vertexCov( row, col ) + pCov( row, col ) + pCov( row + 3, col + 3 ) * a * a -
                          ( pCov( row, col + 3 ) + pCov( col + 3, row ) ) * a;
        }
      }

      auto W_Inv          = W.invChol(); // FIXME: what about failures???
      auto halfdChi2dLam2 = similarity( dir, W_Inv );
      auto decayLength    = dot( dir, W_Inv * flight ) / halfdChi2dLam2;
      auto decayLengthErr = sqrt( 1 / halfdChi2dLam2 );

      return decayLength / decayLengthErr;
    }

    /** @fn    Lifetime
     *  @brief Calculate the lifetime between the particle and the best PV.
     *
     * The implementation here is a very close reproduction of the one in
     * LoKi::DirectionFitBase::iterate, with some condensing of various helper
     * functions into three main parts:
     *
     * 1. The `ctau0` call, which computes a first-order approximation of the
     * lifetime.
     * 2. The `iterate` call, which refines the approximation based on an
     * iterative fit.
     * 3. The `ctau_step` call, which is used by `iterate` to compute the
     * momentum and position updates made during each fit step.
     */
    template <typename Particle, typename VContainer>
    auto Lifetime( VContainer const& primary, Particle const& particle ) const {

      // LoKi fit runs on a 'transported' particle, transporting p1 to position
      // z, saving result to p2; the lifetime fit then acts on p2
      // Transporter is an instance of ParticleTransporter
      // LoKi calls the transported particle 'good'
      // auto [status, transported] = m_transporter->transport( p1, z, p2 );
      // LoKi fitter defines a 'decay' variable as particle.endVertex
      auto ctau = ctau0( primary, particle );

      using float_v = decltype( ctau );
      float_v error = -1.e+10 * Gaudi::Units::mm;
      float_v chi2  = -1.e+10;
      iterate( primary, particle, ctau, error, chi2 );

      auto lifetime = ctau / Gaudi::Units::c_light;
      error /= Gaudi::Units::c_light;
      return std::tuple{lifetime, error, chi2};
    }

  private:
    /** @fn   iterate
     *  @brief Calculate the lifetime between the particle and the best PV.
     *
     */

    template <typename Particle, typename VContainer, typename float_v>
    auto iterate( VContainer const& primary, Particle const& particle, float_v& ctau, float_v& error,
                  float_v& chi2 ) const {

      using Sel::Utils::all;
      using Sel::Utils::endVertexPos;
      using Sel::Utils::fourMomentum;
      using std::abs;
      using std::sqrt;

      // convergence parameters
      float_v delta_chi2 = 0.001;
      float_v delta_ctau = 0.01 * Gaudi::Units::micrometer;
      int     m_max_iter = 20;

      // invariants which are not changed during iteration
      const auto momCov        = momCovMatrix( particle );
      const auto posCov        = posCovMatrix( particle );
      const auto momPosCov     = momPosCovMatrix( particle );
      const auto initMom       = fourMomentum( particle );
      const auto initPos       = endVertexPos( particle );
      const auto primaryPosCov = posCovMatrix( primary );
      const auto primaryPos    = endVertexPos( primary );

      // Copies which will be modified during the iteration
      auto momentum   = initMom;
      auto decvertex  = initPos;
      auto primvertex = primaryPos;

      auto converged = false;
      for ( auto iter = 0; iter < m_max_iter; iter++ ) {
        const auto& [new_ctau, new_chi2, new_error] =
            ctau_step( primaryPos, primaryPosCov, initMom, initPos, momCov, posCov, momPosCov, momentum, decvertex,
                       primvertex, ctau );
        converged = all( ( abs( chi2 - new_chi2 ) < delta_chi2 ) || ( abs( ctau - new_ctau ) < delta_ctau ) );
        if ( converged ) { break; }
        ctau  = new_ctau;
        chi2  = new_chi2;
        error = new_error;
      }
      if ( !converged ) m_owning_algorithm->warning() << "Lifetime fit did not converge. Aborting." << endmsg;
    }

    /** @fn  ctau_step
     *  @brief Calculate one step of the var-fit.
     *
     * The `momentum`, `decvertex`, and `primvertex` inputs are mutated by this
     * method based on updates computed in a single iteration of the fit. The
     * `primary` and `particle` inputs hold the 'reference' points against which
     * the updates will be applied, e.g. `momentum = particle.momentum() +
     * momentum_update;`.
     */
    template <typename MomCov, typename PosCov, typename MomPosCov, typename Vec4D, typename Vec3D, typename float_v>
    auto ctau_step( Vec3D const& primaryPos, PosCov const& primaryPosCov, Vec4D const& initMom, Vec3D const& initPos,
                    MomCov const& momCov, PosCov const& posCov, MomPosCov const& momPosCov, Vec4D& momentum,
                    Vec3D& decvertex, Vec3D& primvertex, float_v const& ctau ) const {

      using std::sqrt;

      auto const px = X( momentum );
      auto const py = Y( momentum );
      auto const pz = Z( momentum );
      auto const e  = E( momentum );
      auto const m2 = e * e - ( px * px + py * py + pz * pz );
      auto const m  = sqrt( m2 );

      // LoKi::Fitters::e_ctau
      auto const vec_E = LHCb::LinAlg::Vec3{px, py, pz} / m;

      LHCb::LinAlg::Mat<float_v, 3, 4> mat_W;
      mat_W( 0, 0 ) = ( 1.0 + px * px / m2 );
      mat_W( 0, 1 ) = ( px * py / m2 );
      mat_W( 0, 2 ) = ( px * pz / m2 );
      mat_W( 1, 0 ) = mat_W( 0, 1 );
      mat_W( 1, 1 ) = ( 1.0 + py * py / m2 );
      mat_W( 1, 2 ) = ( py * pz / m2 );
      mat_W( 2, 0 ) = mat_W( 0, 2 );
      mat_W( 2, 1 ) = mat_W( 1, 2 );
      mat_W( 2, 2 ) = ( 1.0 + pz * pz / m2 );

      mat_W( 0, 3 ) = ( -px * e / m2 );
      mat_W( 1, 3 ) = ( -py * e / m2 );
      mat_W( 2, 3 ) = ( -pz * e / m2 );
      mat_W         = mat_W * ctau / m;

      auto const m_d = LHCb::LinAlg::convert( vec_E * ctau + ( primvertex - decvertex ) );

      auto const mat_VD_tmp = similarity( mat_W, momCov ) + posCov + primaryPosCov - mat_W * momPosCov;
      // FIXME: why is mat_VD_tmp not a symmetrical matrix type anymore?????? Because if it was, we could just use
      // mat_VD_tmp.invChol()...
      auto mat_VD_tmp_gaudi =
          helper::toSMatrix<double>( mat_VD_tmp ); // FIXME: this can only work in the scalar case...
      auto const [inversion_success, mat_VD_gaudi] =
          Sel::vector::invert{}( mat_VD_tmp_gaudi ); // No inverse for LHCb::LinAlg::Mat (yet), so taking
                                                     // converting to a Gaudi3X3 for the inversion for now
      if ( !inversion_success ) {
        m_owning_algorithm->warning() << "Error in Matrix Inversion Type 2." << endmsg;
        return std::tuple<float_v, float_v, float_v>{NonPhysicalValue, NonPhysicalValue, NonPhysicalValue};
      }

      LHCb::LinAlg::Mat<float_v, 3, 3> mat_VD = LHCb::LinAlg::convert<float_v>( mat_VD_gaudi );

      const auto m_Da0 = mat_W * ( initMom - momentum ) - LHCb::LinAlg::convert( initPos - decvertex ) +
                         LHCb::LinAlg::convert( primaryPos - primvertex );

      const auto m_l0 = mat_VD * ( m_Da0 + m_d );

      const auto ctau_variance = 1. / similarity( vec_E, mat_VD );
      if ( ctau_variance < 0. ) {
        m_owning_algorithm->warning() << "Negative variance produced in lifetime fit iteration." << endmsg;
        return std::tuple<float_v, float_v, float_v>{NonPhysicalValue, NonPhysicalValue, NonPhysicalValue};
      }

      const auto delta_ctau = -ctau_variance * dot( vec_E, m_l0 );

      const auto  m_D1 = momCov * mat_W.transpose() - momPosCov;
      const auto  m_D2 = ( mat_W * momPosCov ).transpose() - posCov;
      const auto& m_D3 = primaryPosCov;

      const auto m_l = m_l0 + LHCb::LinAlg::convert( mat_VD * vec_E * delta_ctau );

      const auto delta_momentum    = m_D1 * m_l * -1;
      const auto delta_decay_pos   = toVec3( m_D2 * m_l * -1 );
      const auto delta_primary_pos = toVec3( m_D3 * m_l * -1 );

      // Fill in 'output' values
      auto const updated_ctau = ctau + delta_ctau;
      auto const chi2         = dot( m_l, m_Da0 + m_d );
      auto const error        = sqrt( ctau_variance );

      // Update for the next iteration
      momentum   = initMom + delta_momentum;
      decvertex  = initPos + delta_decay_pos;
      primvertex = primaryPos + delta_primary_pos;

      return std::tuple{updated_ctau, chi2, error};
    }

    /**  @brief Fast, approximate evaluation of c * tau.
     *
     * This neglects the particle momentum covariance, taking into account only
     * the covariances of the primary and secondary vertex positions.
     *
     * Implementation from LoKi::DirectionFitBase::ctau0.
     */
    template <typename Particle, typename VContainer>
    auto ctau0( VContainer const& primary, Particle const& particle ) const {
      // Retrieve position and position covariance of the decay and primary vertices
      using Sel::Utils::endVertexPos;
      using Sel::Utils::mass2;
      using Sel::Utils::posCovMatrix;
      using Sel::Utils::threeMomentum;
      using std::sqrt;
      const auto decay_pos       = endVertexPos( particle );
      const auto decay_pos_cov   = posCovMatrix( particle );
      const auto primary_pos     = endVertexPos( primary );
      const auto primary_pos_cov = posCovMatrix( primary );
      const auto mat_VD_tmp      = primary_pos_cov + decay_pos_cov;
      const auto mat_VD          = mat_VD_tmp.invChol(); // what to do if inversion fails?
      // if ( !inversion_success ) {
      //  m_owning_algorithm->warning() << "Error in Matrix Inversion Type 1." <<
      //  endmsg; return NonPhysicalValue;
      //}

      auto const vec_E = threeMomentum( particle ) / sqrt( mass2( particle ) );
      auto const lam0  = mat_VD * ( primary_pos - decay_pos );
      return -1.0 * dot( vec_E, lam0 ) / similarity( vec_E, mat_VD );
    }
  };
} // namespace Sel
