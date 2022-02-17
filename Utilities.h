/*****************************************************************************\
* (c) Copyright 2019-20 CERN for the benefit of the LHCb Collaboration        *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once
#include "Event/Particle.h"
#include "Event/Particle_v2.h"
#include "Event/Proxy.h"
#include "Event/SOAZip.h"
#include "Event/State.h"
#include "Event/Track_v3.h"
#include "Kernel/Variant.h"
#include "LHCbMath/MatVec.h"
#include "LHCbMath/MatrixTransforms.h"
#include "LHCbMath/MatrixUtils.h"
#include "LHCbMath/SIMDWrapper.h"
#include "SelKernel/ParticleTraits.h"
#include "TrackKernel/TrackVertexUtils.h"

#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/detected.h"

#include <cassert>

/** @file  Utilities.h
 *  @brief Helper types and functions for selections.
 */

// Forward Declaration of type-checks
template <typename... InputTypes>
struct ParticleCombination;

namespace Sel {
  template <typename... child_ts>
    struct ParticleCombination;
} // namespace Sel

namespace LHCb {
  class Particle;
  class VertexBase;
  namespace TrackKernel {
    template <std::size_t NDAUGHTERS, class FTYPE, class TRACKPTR>
    struct TrackCompactVertex;
  } // namespace TrackKernel
} // namespace LHCb

namespace Sel::Utils {
  /** Backwards compatibility -- these moved to Kernel/Variant.h in LHCb
   */
  template <typename... Ts>
  using variant = LHCb::variant<Ts...>;

  template <typename T>
  using is_variant = LHCb::is_variant<T>;

  template <typename T>
  inline constexpr bool is_variant_v = LHCb::is_variant_v<T>;

  template <typename... Variants>
  inline constexpr bool are_all_variants_v = LHCb::are_all_variants_v<Variants...>;

  template <typename... Args>
  [[nodiscard, gnu::always_inline]] inline decltype( auto ) invoke_or_visit( Args&&... args ) {
    return LHCb::invoke_or_visit( std::forward<Args>( args )... );
  }

  /** Backwards compatibility -- this moved to LHCbMath/Utils.h in LHCb
   */
  template <typename Data>
  inline constexpr auto as_arithmetic( Data&& x ) {
    return LHCb::Utils::as_arithmetic( std::forward<Data>( x ) );
  }

  /** Helper to determine if the given type has a static mask_true() method. */
  template <typename T>
  using has_static_mask_true_ = decltype( T::mask_true() );

  template <typename T>
  inline constexpr bool has_static_mask_true_v = Gaudi::cpp17::is_detected_v<has_static_mask_true_, T>;

  template <typename T>
  using has_m_children_ = decltype( std::declval<T>().m_children );

  template <typename T>
  inline constexpr bool has_m_children_v = Gaudi::cpp17::is_detected_v<has_m_children_, T>;

  template <typename T>
  using has_dType = typename T::dType;

  template <typename T>
  inline constexpr bool has_dType_v = Gaudi::cpp17::is_detected_v<has_dType, T>;

  /** Define plain bool versions that mirror the functions defined for
   *  SIMDWrapper's mask_v types. These are useful when writing generic
   *  functor code that works both for scalar and vector types.
   */
  constexpr bool all( bool x ) { return x; }
  constexpr bool any( bool x ) { return x; }
  constexpr bool none( bool x ) { return !x; }
  constexpr int  popcount( bool x ) { return x; }
  template <typename T>
  constexpr T select( bool x, T a, T b ) {
    return x ? a : b;
  }
  template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  constexpr T hmin( T x, bool ) {
    return x;
  }
  template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  constexpr T hmax( T x, bool ) {
    return x;
  }

  template <typename T, typename M>
  constexpr auto hmin( T const& x, M const& m ) {
    return x.hmin( m );
  }
  template <typename T, typename M>
  constexpr auto hmax( T const& x, M const& m ) {
    return x.hmax( m );
  }

  template <typename T>
  auto&& deref_if_ptr( T&& x ) {
    if constexpr ( std::is_pointer_v<std::remove_reference_t<T>> ) {
      assert( x ); // turn on DEBUG flags if you want to check this
      return *x;
    } else {
      return std::forward<T>( x );
    }
  }

  /** Helper to determine if the given type has a bestPV() method. */
  template <typename T>
  using has_bestPV_ = decltype( std::declval<T>().bestPV() );

  template <typename T>
  inline constexpr bool has_bestPV_v = Gaudi::cpp17::is_detected_v<has_bestPV_, T>;

  /** Helper to determine if the given type has a size() method, to determine if one deals with a int_v or int. */
  template <typename T>
  using has_size_ = decltype( std::declval<T>().size() );

  template <typename T>
  inline constexpr bool has_size_v = Gaudi::cpp17::is_detected_v<has_size_, T>;

  template <typename T>
  constexpr bool is_legacy_particle =
      std::is_same_v<LHCb::Particle, std::decay_t<std::remove_pointer_t<std::decay_t<T>>>>;

  namespace details {

    template <std::size_t N, SIMDWrapper::InstructionSet... is>
    constexpr SIMDWrapper::InstructionSet match_size() {
      for ( auto [i, n] : std::array{std::pair{is, SIMDWrapper::type_map<is>::type::size}...} ) {
        if ( n == N ) return i;
      }
      throw std::invalid_argument( "unable to map vectorwidth to instruction set" );
      return SIMDWrapper::InstructionSet( -1 );
    }

    template <std::size_t N>
    constexpr SIMDWrapper::InstructionSet
        instructionSet_for_ = match_size<N, SIMDWrapper::InstructionSet::AVX512, SIMDWrapper::InstructionSet::AVX256,
                                         SIMDWrapper::InstructionSet::AVX2, SIMDWrapper::InstructionSet::SSE>();

    template <typename Target, typename... Ts>
    struct contains : std::disjunction<std::is_same<Target, Ts>...> {};

  } // namespace details

  template <typename T>
  constexpr bool is_lhcb_vertexbase = std::is_base_of_v<LHCb::VertexBase, T>;

  template <typename T>
  struct is_trackcompactvertex : std::false_type {};

  template <std::size_t NDAUGHTERS, class FTYPE, class TRACKPTR>
  struct is_trackcompactvertex<LHCb::TrackKernel::TrackCompactVertex<NDAUGHTERS, FTYPE, TRACKPTR>> : std::true_type {};

  template <typename>
    struct is_selparticlecombination : std::false_type {};
  template <typename... child_ts>
    struct is_selparticlecombination<const Sel::ParticleCombination<child_ts...>> : std::true_type {};

  template <typename>
    struct is_combparticlecombination : std::false_type {};

  template <typename... InputTypes>
    struct is_combparticlecombination<const ParticleCombination<InputTypes...>> : std::true_type {};

  inline constexpr auto get_track_property_from_particle = []( auto const& p, auto&& accessor, auto&& invalid ) {
    auto const* pp  = deref_if_ptr( p ).proto();
    auto const* trk = pp ? pp->track() : nullptr;
    return trk ? accessor( trk ) : invalid;
  };

  /* when we have different accessor names, we use ADL, with as default this indirection layer */
  template <typename T>
  auto trackState( T const& item ) -> decltype( item.trackState() ) {
    return item.trackState();
  }
  template <typename T>
  auto slopes( T const& item ) -> decltype( item.slopes() ) {
    return item.slopes();
  }
  template <typename T>
  auto threeMomentum( T const& item ) -> decltype( item.threeMomentum() ) {
    return item.threeMomentum();
  }
  template <typename T>
  auto fourMomentum( T const& item ) -> decltype( item.momentum() ) {
    return item.momentum();
  }
  template <typename T>
  auto mass2( T const& item ) -> decltype( item.mass2() ) {
    return item.mass2();
  }
  template <typename T>
  auto endVertexPos( T const& item ) -> decltype( item.endVertex() ) {
    return item.endVertex();
  }

  template <typename T, std::size_t N>
  auto endVertexPos( std::array<T const*, N> const& items ) {
    using float_v = typename SIMDWrapper::type_map<details::instructionSet_for_<N>>::type::float_v;
    return std::apply(
        []( const auto&... i ) {
          constexpr auto nan     = std::numeric_limits<float>::quiet_NaN();
          constexpr auto invalid = LHCb::LinAlg::Vec3<SIMDWrapper::scalar::float_v>{nan, nan, nan};
          using LHCb::LinAlg::gather;
          return gather<float_v>( std::array{( i ? endVertexPos( *i ) : invalid )...} );
        },
        items );
  }

  template <typename T>
  auto referencePoint( T const& item ) -> decltype( item.referencePoint() ) {
    return item.referencePoint();
  }
  template <typename T>
  auto threeMomCovMatrix( T const& item ) -> decltype( item.threeMomCovMatrix() ) {
    return item.threeMomCovMatrix();
  }

  // TODO: make sure this predicate becomes obsolete by moving the 'dispatch' code
  //       into the relevant event model classes -- and the code to which it dispatches
  //       (which is written in terms of more basic quantities like positions, directions,
  //       and covariance matrices) into LHCbMath
  template <typename T>
  constexpr bool hasVertex_v = std::remove_pointer_t<std::decay_t<T>>::hasVertex;

  template <typename T>
  auto threeMomPosCovMatrix( T const& item ) -> decltype( item.threeMomPosCovMatrix() ) {
    return item.threeMomPosCovMatrix();
  }
  template <typename T>
  auto momPosCovMatrix( T const& item ) -> decltype( item.momPosCovMatrix() ) {
    return item.momPosCovMatrix();
  }
  template <typename T>
  auto momCovMatrix( T const& item ) -> decltype( item.momCovMatrix() ) {
    return item.momCovMatrix();
  }
  template <typename T>
  auto posCovMatrix( T const& item ) -> decltype( item.posCovMatrix() ) {
    return item.posCovMatrix();
  }
  template <typename T, std::size_t N>
  auto posCovMatrix( std::array<T const*, N> const& items ) {
    using float_v = typename SIMDWrapper::type_map<details::instructionSet_for_<N>>::type::float_v;
    return std::apply(
        []( const auto&... i ) {
          constexpr auto nan     = std::numeric_limits<float>::quiet_NaN();
          constexpr auto invalid = LHCb::LinAlg::MatSym<SIMDWrapper::scalar::float_v, 3>{nan, nan, nan, nan, nan, nan};
          using LHCb::LinAlg::gather;
          return gather<float_v>( std::array{( i ? posCovMatrix( *i ) : invalid )...} );
        },
        items );
  }

  template <typename T>
  __attribute__( ( flatten ) ) auto covMatrix( T const& item ) {
    // LHCb::LinAlg::resize_t<decltype( posCovMatrix( obj ) ), 6> cov6{};
    // cov6 = cov6.template place_at<0, 0>( posCovMatrix( obj ) );
    // cov6 = cov6.template place_at<3, 0>( threeMomPosCovMatrix( obj ) );
    // cov6 = cov6.template place_at<3, 3>( threeMomCovMatrix( obj ) );
    // return cov6;
    auto                                              tempPosMat    = posCovMatrix( item );
    auto                                              tempMomMat    = threeMomCovMatrix( item );
    auto                                              tempMomPosMat = threeMomPosCovMatrix( item );
    LHCb::LinAlg::resize_t<decltype( tempPosMat ), 6> covMat;
    LHCb::Utils::unwind<0, 3>( [&]( auto i ) {
      LHCb::Utils::unwind<0, 3>( [&]( auto j ) {
        covMat( i, j )         = tempPosMat( i, j );
        covMat( i + 3, j + 3 ) = tempMomMat( i, j );
        covMat( i, j + 3 )     = tempMomPosMat( i, j );
        covMat( j, i + 3 )     = tempMomPosMat( j, i );
      } );
    } );
    return covMat;
  }

  /** Helpers for dispatching to the right fdchi2 calculation. */
  template <typename Vertex1, typename Vertex2>
  auto flightDistanceChi2( Vertex1 const& v1, Vertex2 const& v2 ) {
    auto cov = posCovMatrix( v1 ) + posCovMatrix( v2 );
    cov      = cov.invChol(); // what if it fails?
    return similarity( endVertexPos( v1 ) - endVertexPos( v2 ), cov );
  }

  /** @fn    impactParameterChi2
   *  @brief Helper for dispatching to the correct ipchi2 calculation.
   *
   * The input be particle-like, with an (x, y, z) position and 3/4-momentum, or
   * it could be state-like, with a (x, y, tx, ty[, q/p]) vector and covariance
   * matrix.
   *
   * LHCb::TrackVertexUtils::vertexChi2() has the ipchi2 calculation for a
   * state and [primary] vertex position/covariance.
   *
   * @todo Add a [template] version of LHCb::TrackVertexUtils::vertexChi2()
   *       that takes a particle-like? It only uses the 4x4 upper corner of
   *       the state-like covariance matrix, so we might be able to save
   *       something?
   */

  template <typename TrackOrParticle, typename Vertex>
  __attribute__( ( flatten ) ) auto impactParameterChi2( TrackOrParticle const& obj, Vertex const& vertex ) {
    if constexpr ( std::is_pointer_v<TrackOrParticle> ) {
      assert( obj );
      return impactParameterChi2( *obj, vertex );
    } else if constexpr ( !hasVertex_v<TrackOrParticle> ) {
      auto chi2 =
          LHCb::TrackVertexUtils::vertexChi2( trackState( obj ), endVertexPos( vertex ), posCovMatrix( vertex ) );
      return chi2;
    } else {
      // composite with a vertex
      auto [chi2, decaylength, decaylength_err] =
          LHCb::TrackVertexUtils::computeChiSquare( referencePoint( obj ), threeMomentum( obj ), covMatrix( obj ),
                                                    endVertexPos( vertex ), posCovMatrix( vertex ) );
      return chi2;
    }
  }

  template <typename Vertex>
  auto impactParameterChi2( LHCb::Particle const& particle, Vertex const& vertex ) {
    auto const* pp    = particle.proto();
    auto const* track = ( pp ? pp->track() : nullptr );
    if ( track ) return impactParameterChi2( *track, vertex );
    auto [chi2, decaylength, decaylength_err] = LHCb::TrackVertexUtils::computeChiSquare(
        referencePoint( particle ), threeMomentum( particle ), covMatrix( particle ), endVertexPos( vertex ),
        posCovMatrix( vertex ) );
    return chi2;
  }
} // namespace Sel::Utils
