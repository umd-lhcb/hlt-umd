/*****************************************************************************\
* (c) Copyright 2019-2021 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once
#include "Event/ProtoParticle.h"
#include "Event/Track_v3.h"
#include "Functors/Function.h"
#include "Functors/Utilities.h"
#include "GaudiKernel/detected.h"
#include "SelKernel/ParticleAccessors.h"
#include "SelKernel/Utilities.h"
#include "SelKernel/VertexRelation.h"

/** @file  TrackLike.h
 *  @brief Definitions of functors for track-like objects.
 */

/** @namespace Functors::Track
 *
 *  This defines some basic functors that operate on track-like, or maybe more
 *  accurately charged-particle-like, objects. Both predicates and functions
 *  are defined in the same file for the moment.
 */

namespace Functors {
  static constexpr int invalid_value = -1000;
  using Sel::Utils::get_kinematic_state;
  using Sel::Utils::has_tracklike_API;
  using Sel::Utils::is_legacy_particle;
} // namespace Functors

namespace Functors::detail {
  /*
   * Helper for ProbNN
   */
  enum struct Pid { electron, muon, pion, kaon, proton, deuteron, ghost };

  constexpr LHCb::ProtoParticle::additionalInfo to_ppai( Pid pid ) {
    switch ( pid ) {
    case Pid::electron:
      return LHCb::ProtoParticle::additionalInfo::ProbNNe;
    case Pid::muon:
      return LHCb::ProtoParticle::additionalInfo::ProbNNmu;
    case Pid::pion:
      return LHCb::ProtoParticle::additionalInfo::ProbNNpi;
    case Pid::kaon:
      return LHCb::ProtoParticle::additionalInfo::ProbNNk;
    case Pid::proton:
      return LHCb::ProtoParticle::additionalInfo::ProbNNp;
    case Pid::deuteron:
      return LHCb::ProtoParticle::additionalInfo::ProbNNd;
    case Pid::ghost:
      return LHCb::ProtoParticle::additionalInfo::ProbNNghost;
    }
    __builtin_unreachable(); // suppress gcc warning -Wreturn-type
  }

  /** @brief helpers for probNN
   */
  template <Pid id>
  struct has_probNN {
    template <typename T>
    using check_for_probNN_id = decltype( std::declval<T const&>().template probNN<id>() );
    template <typename T>
    static constexpr bool value = Gaudi::cpp17::is_detected_v<check_for_probNN_id, T>;
  };

  template <Pid id, typename T>
  constexpr bool has_probNN_v = has_probNN<id>::template value<T>;

  template <typename StatePos_t, typename StateDir_t, typename VertexPos_t>
  auto impactParameterSquared( StatePos_t const& state_pos, StateDir_t const& state_dir, VertexPos_t const& vertex ) {
    return ( state_pos - vertex ).Cross( state_dir ).mag2() / state_dir.mag2();
  }
  template <typename T>
  using check_for_proto = decltype( std::declval<T const&>().proto() );

  template <typename T>
  constexpr bool has_proto_v = Gaudi::cpp17::is_detected_v<check_for_proto, T>;

  template <typename T>
  constexpr bool probnn_always_false = false;

} // namespace Functors::detail

namespace Functors::Track {
  /** @brief x coordinate, as defined by the x() accessor.
   */
  struct X : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.x( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).x();
    }
  };

  /** @brief y coordinate, as defined by the y() accessor.
   */
  struct Y : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.y( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).y();
    }
  };

  /** @brief z coordinate, as defined by the z() accessor.
   */
  struct Z : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.z( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).z();
    }
  };

  /** @brief x slope, as defined by the tx() accessor.
   */
  struct TX : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.tx( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).tx();
    }
  };

  /** @brief y slope, as defined by the ty() accessor.
   */
  struct TY : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.ty( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).ty();
    }
  };

  /** @brief Given element of the covariance matrix.
   */
  struct Covariance : public Function {
    Covariance( int i, int j ) : m_i{i}, m_j{j} {}

    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.covariance( Sel::Utils::get_kinematic_state( d ) )( m_i, m_j );
      else
        return Sel::Utils::deref_if_ptr( d ).covariance()( m_i, m_j );
    }

  private:
    int m_i{0}, m_j{0};
  };

  /** @brief Charge, as defined by the charge() accessor.
   */
  struct Charge : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      return Sel::Utils::deref_if_ptr( d ).charge();
    }
  };

  /** @brief Momentum, as defined by the p() accessor.
   */
  struct Momentum : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.p( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).p();
    }
  };

  /** @brief Azimuthal angle, as defined by the phi() accessor.
   */
  struct Phi : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> )
        return Sel::Utils::deref_if_ptr( d ).momentum().Phi();
      else if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.phi( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).phi();
    }
  };

  /** @brief Transverse momentum, as defined by the pt() accessor.
   */
  struct TransverseMomentum : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.pt( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).pt();
    }
  };

  struct MyTransverseMomentum : public Function {
    template <typename Data>
      auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
		     return d.pt( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).pt();
    }
  };

  /** @brief IsMuon, as defined by the accessor of the same name.
   */
  struct IsMuon : public Predicate {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> ) {
        auto const* pp = Sel::Utils::deref_if_ptr( d ).proto();
        if ( not pp ) { return false; }
        auto pid = pp->muonPID();
        if ( pid ) return pid->IsMuon();
        return bool( pp->info( std::decay_t<decltype( *pp )>::additionalInfo::MuonPIDStatus, 0 ) );
      } else
        return Sel::Utils::deref_if_ptr( d ).IsMuon();
    }
  };

  /** @brief eta/pseudorapidity, as defined by the pseudoRapidity accessor.
   */
  struct PseudoRapidity : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> ) {
        return Sel::Utils::deref_if_ptr( d ).momentum().eta();
      } else if constexpr ( Sel::Utils::has_tracklike_API<Data> ) {
        return d.pseudoRapidity( Sel::Utils::get_kinematic_state( d ) );
      } else {
        return d.pseudoRapidity();
      }
    }
  };

  /** @brief Number of degrees of freedom, as defined by the nDoF accessor.
   */
  struct nDoF : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      return Sel::Utils::deref_if_ptr( d ).nDoF();
    }
  };

  /** @brief chi^2/d.o.f., as defined by the chi2PerDoF accessor.
   *
   * If the input is a legacy LHCb::Particle with a track, the track object's
   * accessor is used. If a track is not present but a vertex is the vertex's
   * accessor is used.
   */
  struct Chi2PerDoF : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      auto const& particle = Sel::Utils::deref_if_ptr( d );
      using Particle       = std::decay_t<decltype( particle )>;
      if constexpr ( Sel::Utils::is_legacy_particle<Particle> ) {
        // Try accessing a track first
        const auto value = Sel::Utils::get_track_property_from_particle(
            d, []( auto&& t ) { return t->chi2PerDoF(); }, invalid_value );
        // If that failed, try accessing a vertex
        if ( value == invalid_value && particle.endVertex() ) { return particle.endVertex()->chi2PerDoF(); }
        // Input is either a charged or a neutral basic
        return value;
      } else {
        return Sel::get::chi2PerDoF( particle );
      }
    }
  };

  /** @brief chi^2, as defined by the chi2 accessor.
   *
   * If the input is a legacy LHCb::Particle with a track, the track object's
   * accessor is used. If a track is not present but a vertex is the vertex's
   * accessor is used.
   */
  struct Chi2 : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      auto const& particle = Sel::Utils::deref_if_ptr( d );
      using Particle       = std::decay_t<decltype( particle )>;
      if constexpr ( is_legacy_particle<Particle> ) {
        // Try accessing a track first
        const auto value = Sel::Utils::get_track_property_from_particle(
            particle, []( auto&& t ) { return t->chi2(); }, invalid_value );
        // If that failed, try accessing a vertex
        if ( value == invalid_value && particle.endVertex() ) { return particle.endVertex()->chi2(); }
        // Input is either a charged or a neutral basic
        return value;
      } else {
        return particle.chi2();
      }
    }
  };

  /** @brief q/p, as defined by the qOverP accessor.
   */
  struct QoverP : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( Sel::Utils::has_tracklike_API<Data> )
        return d.qOverP( Sel::Utils::get_kinematic_state( d ) );
      else
        return Sel::Utils::deref_if_ptr( d ).qOverP();
    }
  };

  /** @brief Ghost probability, as defined by the ghostProbability accessor.
   */
  struct GhostProbability : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> ) {
        return Sel::Utils::get_track_property_from_particle(
            d, []( auto&& t ) { return t->ghostProbability(); }, invalid_value );
      } else {
        return Sel::Utils::deref_if_ptr( d ).ghostProbability();
      }
    }
  };

  template <typename F>
  struct ClosestToBeamState : public std::conditional_t<detail::is_functor_predicate_v<F>, Predicate, Function> {
    static_assert( detail::is_functor_function_v<F>, "ClosestToBeamState adapter must wrap a functor!" );
    ClosestToBeamState( F f ) : m_f{std::move( f )} {}

    /* Improve error messages. */
    constexpr auto name() const { return "ClosestToBeamState( " + detail::get_name( m_f ) + " )"; }

    void bind( TopLevelInfo& top_level ) { detail::bind( m_f, top_level ); }

    auto prepare( EventContext const& evtCtx, TopLevelInfo const& top_level ) const {
      return [f = detail::prepare( m_f, evtCtx, top_level )]( const auto& track ) {
        using TrackType = std::remove_cv_t<std::remove_reference_t<decltype( track )>>;
        if constexpr ( Sel::Utils::has_tracklike_API<TrackType> )
          return f( track.state( TrackType::StateLocation::ClosestToBeam ) );
        else
          return f( track.closestToBeamState() );
      };
    }

  private:
    F m_f;
  };

  // template <typename F>
  // struct ClosestToBeamState : public std::conditional_t<detail::is_functor_predicate_v<F>, Predicate, Function> {
  //  static_assert( detail::is_functor_function_v<F>, "ClosestToBeamState adapter must wrap a functor!" );
  //  ClosestToBeamState( F f ) : m_f{std::move( f )} {}

  //  /* Improve error messages. */
  //  constexpr auto name() const { return "ClosestToBeamState( " + detail::get_name( m_f ) + " )"; }

  //  void bind( TopLevelInfo& top_level ) { detail::bind( m_f, top_level ); }

  //  auto prepare( EventContext const& evtCtx, TopLevelInfo const& top_level ) const {
  //    return [f = detail::prepare( m_f, evtCtx, top_level )]( auto const& track ) {
  //        return f( track.closestToBeamState() );
  //    };
  //  }

  // private:
  //  F m_f;
  //};

  /** @brief Ghost probability, as defined by the ghostProbability accessor.of the track
   */
  struct TrGhostProbability : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      return Sel::Utils::get_track_property_from_particle(
          d, []( auto&& t ) { return t->ghostProbability(); }, invalid_value );
    }
  };

  /** @brief PIDmu, as defined by the CombDLLmu variable
   */
  struct PIDmu : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> ) {
        auto const* pp = Sel::Utils::deref_if_ptr( d ).proto();
        if ( not pp ) { return double( invalid_value ); }
        return pp->info( std::decay_t<decltype( *pp )>::additionalInfo::CombDLLmu, invalid_value );
      } else {
        return d.CombDLLmu();
      }
    }
  };

  /** @brief PIDp, as defined by the CombDLLp variable
   */
  struct PIDp : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> ) {
        auto const* pp = Sel::Utils::deref_if_ptr( d ).proto();
        if ( not pp ) { return double( invalid_value ); }
        return pp->info( std::decay_t<decltype( *pp )>::additionalInfo::CombDLLp, invalid_value );
      } else {
        return d.CombDLLp();
      }
    }
  };

  /** @brief PIDe, as defined by the CombDLLe variable
   */
  struct PIDe : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> ) {
        auto const* pp = Sel::Utils::deref_if_ptr( d ).proto();
        if ( not pp ) { return double( invalid_value ); }
        return pp->info( std::decay_t<decltype( *pp )>::additionalInfo::CombDLLe, invalid_value );
      } else {
        return d.CombDLLe();
      }
    }
  };

  /** @brief PIDk, as defined by the CombDLLk variable
   */
  struct PIDk : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> ) {
        auto const* pp = Sel::Utils::deref_if_ptr( d ).proto();
        if ( not pp ) { return double( invalid_value ); }
        return pp->info( std::decay_t<decltype( *pp )>::additionalInfo::CombDLLk, invalid_value );
      } else {
        return d.CombDLLk();
      }
    }
  };

  /** @brief PIDpi, as defined by the CombDLLpi variable
   */
  struct PIDpi : public Function {
    template <typename Data>
    auto operator()( Data const& d ) const {
      if constexpr ( is_legacy_particle<Data> ) {
        auto const* pp = Sel::Utils::deref_if_ptr( d ).proto();
        if ( not pp ) { return double( invalid_value ); }
        return pp->info( std::decay_t<decltype( *pp )>::additionalInfo::CombDLLpi, invalid_value );
      } else {
        return d.CombDLLpi();
      }
    }
  };

  /** @brief ProbNN templated definition
   */
  template <detail::Pid id>
  struct ProbNN : public Function {
    template <typename T>
    auto operator()( const T& d ) const {
      if constexpr ( std::is_pointer_v<T> ) {
        return operator()( *d ); // return itself with class
      } else if constexpr ( detail::has_probNN_v<id, T> ) {
        return std::invoke( &T::template probNN<id>, d );
      } else if constexpr ( detail::has_proto_v<T> ) {
        LHCb::ProtoParticle const* pp = d.proto();
        return pp ? pp->info( to_ppai( id ), invalid_value ) : invalid_value;
      } else {
        static_assert( detail::probnn_always_false<T>, "The type T neither has a `proto()` member function nor a "
                                                       "`probNN<id>`() member function -- sorry, not supported" );
        return -1.;
      }
    }
  };

  /** @brief The explicit definition for the probNN quantities
   * the interpreter doesn't like that. See
   * https://gitlab.cern.ch/lhcb/Rec/-/merge_requests/2471#note_4863872
   *
   *    constexpr auto PROBNN_D     = ProbNN<detail::Pid::deuteron>{};
   *    constexpr auto PROBNN_E     = ProbNN<detail::Pid::electron>{};
   *    constexpr auto PROBNN_GHOST = ProbNN<detail::Pid::ghost>{};
   *    constexpr auto PROBNN_K     = ProbNN<detail::Pid::kaon>{};
   *    constexpr auto PROBNN_MU    = ProbNN<detail::Pid::muon>{};
   *    constexpr auto PROBNN_PI    = ProbNN<detail::Pid::pion>{};
   *    constexpr auto PROBNN_P     = ProbNN<detail::Pid::proton>{};
   */
  /** @brief The explicit definition for the probNN quantities
   * Warning: these are types, so you need a {} to get an instance.
   */
  using PROBNN_D_t     = ProbNN<detail::Pid::deuteron>;
  using PROBNN_E_t     = ProbNN<detail::Pid::electron>;
  using PROBNN_GHOST_t = ProbNN<detail::Pid::ghost>;
  using PROBNN_K_t     = ProbNN<detail::Pid::kaon>;
  using PROBNN_MU_t    = ProbNN<detail::Pid::muon>;
  using PROBNN_PI_t    = ProbNN<detail::Pid::pion>;
  using PROBNN_P_t     = ProbNN<detail::Pid::proton>;

  /** @brief nHits, as defined by the nHits accessor.
   */
  struct nHits : public Function {
    static constexpr auto name() { return "nHits"; }
    template <typename Data>
    auto operator()( Data const& d ) const {
      return Sel::Utils::deref_if_ptr( d ).nHits();
    }
  };
} // namespace Functors::Track

namespace Functors::detail {

  struct MinimumImpactParameter : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "MinimumImpactParameter"; }

    template <typename VContainer, typename TrackChunk>
    auto operator()( VContainer const& vertices, TrackChunk const& track_chunk ) const {
      auto const& track     = Sel::Utils::get_track_from_particle_or_track<TrackChunk>( track_chunk );
      auto&&      state     = track.closestToBeamState();
      auto&&      state_pos = state.position();
      auto&&      state_dir = state.slopes();
      using std::min;  // so min() below works when float_v is plain float or double
      using std::sqrt; // ...
      using vector_t = std::decay_t<decltype( state_pos )>;
      using float_v  = decltype( state_pos.X() );
      auto min_ip2 =
          std::accumulate( std::begin( vertices ), std::end( vertices ), float_v{std::numeric_limits<float>::max()},
                           [&]( float_v ip2, auto const& vertex ) {
                             auto const& PV = vertex.position();
                             auto        this_ip2 =
                                 impactParameterSquared( state_pos, state_dir, vector_t( PV.x(), PV.y(), PV.z() ) );
                             static_assert( std::is_same_v<decltype( this_ip2 ), float_v> );
                             return min( ip2, this_ip2 );
                           } );
      return sqrt( min_ip2 );
    }
  };

  struct MinimumImpactParameterChi2 : public Function {
    template <typename VContainer, typename TrackChunk>
    auto operator()( VContainer const& vertices, TrackChunk const& track_chunk ) const {
      // Make sure min() works for plain floating point types
      using std::min;
      // Accept generic containers
      using std::begin;
      using std::end;
      auto const& track = Sel::Utils::get_track_from_particle_or_track<TrackChunk>( track_chunk );
      // Figure out the return type (probably float or SIMDWrapper's float_v for now)
      using float_v = decltype( Sel::Utils::impactParameterChi2( track, vertices.front() ) );
      static_assert( std::numeric_limits<float_v>::is_specialized,
                     "Dangerous use of std::numeric_limits<T> for unspecialised T." );
      return std::accumulate( begin( vertices ), end( vertices ), std::numeric_limits<float_v>::max(),
                              [&track]( float_v ipchi2, auto const& vertex ) {
                                return min( ipchi2, Sel::Utils::impactParameterChi2( track, vertex ) );
                              } );
    }
  };

  struct ImpactParameterChi2ToVertex : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "ImpactParameterChi2ToVertex"; }

    template <typename VContainer, typename Particle>
    auto operator()( VContainer const& vertices, Particle const& particle ) const {
      // Get the associated PV -- this uses a link if it's available and
      // computes the association if it's not.
      return Sel::getBestPVRel( particle, vertices ).ipchi2();
    }
  };

  struct MinimumImpactParameterCut : public Predicate {
    /** Flag that this adapter expects to be invoked with an explicit mask as
     *  its first argument.
     */
    static constexpr bool requires_explicit_mask = true;

    MinimumImpactParameterCut( float cut_value ) : m_cut_value_squared{cut_value * cut_value} {}

    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "MinimumImpactParameterCut"; }

    template <typename mask_t, typename VContainer, typename Particle>
    auto operator()( mask_t const& particle_mask, VContainer const& vertices, Particle const& particle ) const {
      auto const& state     = particle.closestToBeamState();
      auto const& state_pos = state.position();
      auto const& state_dir = state.slopes();
      using Sel::Utils::none; // so none() below works when mask_v is plain bool
      using std::min;         // so min() below works when float_v is plain float or double
      using vector_t = std::decay_t<decltype( state_pos )>;
      using float_v  = decltype( state_pos.X() );
      // This is very important, as for unspecialised types then ::max() returns
      // T(), which might be an uninitialised value
      static_assert( std::numeric_limits<float_v>::is_specialized,
                     "Dangerous use of std::numeric_limits<T> for unspecialised T." );
      // Return the mask we were given if there are no vertices
      auto    local_mask = particle_mask;
      float_v min_d{std::numeric_limits<float_v>::max()};
      for ( auto const& vertex : vertices ) {
        auto const& PV = vertex.position();
        min_d =
            min( min_d, detail::impactParameterSquared( state_pos, state_dir, vector_t( PV.x(), PV.y(), PV.z() ) ) );
        local_mask = local_mask && ( min_d > m_cut_value_squared );
        // Not completely obvious the overhead of short-circuiting is a net
        // gain for the vectorised version, but let's do it anyway for now.
        if ( none( local_mask ) ) { break; }
      }
      return local_mask;
    }

  private:
    float m_cut_value_squared = 0.f;
  };

  struct MinimumImpactParameterChi2Cut : public Predicate {
    /** Flag that this adapter expects to be invoked with an explicit mask as
     *  its first argument.
     */
    static constexpr bool requires_explicit_mask = true;

    MinimumImpactParameterChi2Cut( float cut_value ) : m_cut_value{cut_value} {}

    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "MinimumImpactParameterChi2Cut"; }

    template <typename mask_t, typename VContainer, typename Particle>
    auto operator()( mask_t const& mask, VContainer const& vertices, Particle const& particle ) const {
      // mask tells us how much of the 'width' of 'track_chunk' we should
      // operate on.
      // Figure out what the types are. It would be nice if we could retire this part...
      using float_v = decltype( Sel::Utils::impactParameterChi2( particle, vertices.front() ) );
      // This is very important, as for unspecialised types then ::max() returns
      // T(), which might be an uninitialised value
      static_assert( std::numeric_limits<float_v>::is_specialized,
                     "Dangerous use of std::numeric_limits<T> for unspecialised T." );
      // Keep track of the smallest ipchi2 as we loop, abort if/when all values
      // are below the cut value
      float_v min_ipchi2{std::numeric_limits<float_v>::max()};
      // Make sure min() works for basic arithmetic types
      using std::min;
      // Make sure none() works for plain bool
      using Sel::Utils::none;
      // pass through the mask we were given if there are no vertices
      mask_t local_mask = mask;
      // Unfortunately it seems difficult to avoid the explicit loop...
      for ( auto const& vertex : vertices ) {
        min_ipchi2 = min( min_ipchi2, Sel::Utils::impactParameterChi2( particle, vertex ) );
        local_mask = local_mask && ( min_ipchi2 > m_cut_value );
        // Check if all members of 'particle' have had an ipchi2 below
        // threshold w.r.t. one of the vertices we've already looped over.
        // If so, we can safely abort the loop over vertices.
        if ( none( local_mask ) ) { break; }
      }
      return local_mask;
    }

  private:
    float m_cut_value = 0.f;
  };

} // namespace Functors::detail

namespace Functors::Track {
  /** @brief Minimum impact parameter w.r.t. any of the given vertices.
   *
   *  Note that for the typical MINIP > x cut it is more efficient to use the
   *  dedicated functor that can abort the loop over vertices early.
   */
  template <typename VContainer = detail::DefaultPVContainer_t>
  auto MinimumImpactParameter( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::MinimumImpactParameter, VContainer>{std::move( vertex_location )};
  }

  /** @brief Minimum impact parameter chi2 w.r.t. any of the given vertices.
   *
   *  Note that for the typical MINIPCHI2 > x cut it is more efficiency to use
   *  the dedicated functor that can abort the loop over vertices early.
   */
  template <typename VContainer = detail::DefaultPVContainer_t>
  auto MinimumImpactParameterChi2( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::MinimumImpactParameterChi2, VContainer>{
        std::move( vertex_location )};
  }

  /** @brief Impact parameter chi2 to the "best" one of the given vertices.
   *
   *  Note that if the given data object contains a vertex link then that will
   *  be checked for compatibility with the given vertex container and, if it
   *  matches, be used. The link caches the ipchi2 value, so this should be the
   *  most efficient way of getting the ipchi2 value.
   */
  template <typename VContainer = detail::DefaultPVContainer_t>
  auto ImpactParameterChi2ToVertex( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::ImpactParameterChi2ToVertex, VContainer>{
        std::move( vertex_location )};
  }

  /** @brief Require the minimum impact parameter w.r.t. any of the given
   *         vertices is greater than some threshold.
   *
   *  @param cut_value    Minimum impact parameter value.
   *  @param tes_location TES location of vertex container.
   *
   *  Note this is more efficient than computing the minimum and then cutting
   *  on it, as this functor can short circuit.
   */
  template <typename VContainer = detail::DefaultPVContainer_t>
  auto MinimumImpactParameterCut( float cut_value, std::string vertex_location ) {
    return detail::add_data_deps<Predicate, VContainer>( detail::MinimumImpactParameterCut{cut_value},
                                                         std::move( vertex_location ) );
  }

  /** @brief Require the minimum impact parameter chi2 w.r.t. any of the given
   *         vertices is greater than some threshold.
   *
   *  @param cut_value    Minimum impact parameter chi2 value.
   *  @param tes_location TES location of vertex container.
   *
   *  Note this is more efficient than computing the minimum and then cutting
   *  on it, as this functor can short circuit.
   */
  template <typename VContainer = detail::DefaultPVContainer_t>
  auto MinimumImpactParameterChi2Cut( float cut_value, std::string vertex_location ) {
    return detail::add_data_deps<Predicate, VContainer>( detail::MinimumImpactParameterChi2Cut{cut_value},
                                                         std::move( vertex_location ) );
  }
} // namespace Functors::Track
