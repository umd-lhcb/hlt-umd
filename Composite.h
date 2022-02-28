/***************************************************************************** \
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
#include "Event/State.h"
#include "Functors/Core.h"
#include "Functors/Function.h"
#include "Functors/Utilities.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ParticleProperty.h"
#include "LHCbMath/MatVec.h"
#include "SelKernel/Utilities.h"
#include "SelKernel/VertexRelation.h"

//from old Composite.h
#include "Event/Particle.h"
#include "Event/Particle_v2.h"
//#include "SelKernel/ParticleTraits.h"
//#include "SelKernel/VectorOps.h"
//#include "TrackKernel/TrackCompactVertex.h"
#include "Functors/TrackLike.h"

/** @file  Composite.h
 *  @brief Definitions of functors for composite-particle-like objects.
 */

/** @namespace Functors::Composite
 *
 *  Functors that make sense for composite particles (i.e. ones that have a
 *  vertex as well as a trajectory)
 */
namespace Functors::detail {
  /**MTDOCACHI2**/
  template<int N>
  struct MotherTrajectoryDistanceOfClosestApproachChi2 : public Function {
    MotherTrajectoryDistanceOfClosestApproachChi2( std::integral_constant<int, N> ) {}
    void bind( TopLevelInfo& top_level ){
      m_dist_calc.emplace(top_level.algorithm() );
    }
    template <typename VContainer, typename Particle>
      auto operator()( VContainer const& vertices, Particle const& composite ) const {
      auto const bestPV = Sel::getBestPV(composite, vertices);
      if constexpr (Functors::is_legacy_particle<Particle>) {
	  const auto children( composite.daughtersVector() );
          const auto pN = *children[0];
          //move mother to PV                                                                
	  std::unique_ptr<LHCb::Particle> tempMother( composite.clone() );
          tempMother->setReferencePoint( bestPV.position() );
          tempMother->setPosCovMatrix( bestPV.covMatrix() );
          const auto& dist_calc = *m_dist_calc;
          return dist_calc.particleDOCAChi2(pN, *tempMother);
        }
      else {
        throw GaudiException{"v2 particles not yet supported.", "Functors::detail:MotherTrajectoryDistanceOfClosestApproachChi2", StatusCode::FAILURE};
      }
    }
  private:
    std::optional<Functors::detail::DefaultDistanceCalculator_t> m_dist_calc;
  };


  struct FlightDistanceChi2ToVertex : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "FlightDistanceChi2ToVertex"; }

    template <typename VContainer, typename Particle>
      auto operator()( VContainer const& vertices, Particle const& composite ) const {
      // Get the associated PV and calculate the flight distance chi2 between 'composite' and 'bestPV'
      return Sel::Utils::flightDistanceChi2( composite, Sel::getBestPV( composite, vertices ) );
    }
  };

  /** BPVVDZ */
  /** origin version in Loki defined in PHYS/PHYS_v25r1/Phys/LoKiPhys/src/Particles20.cpp,
      VertexZDistanceWithTheBestPV function */
  struct DeltaZToVertex : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "DeltaZToVertex"; }

    template <typename VContainer, typename Particle>
      auto operator()( VContainer const& vertices, Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return ( endVertexPos( composite ) - endVertexPos( Sel::getBestPV( composite, vertices ) ) ).Z();
    }
  };

  /** BPVVDRHO */
  /** origin version in Loki defined in PHYS/PHYS_v25r1/Phys/LoKiPhys/src/Particles20.cpp,
      VertexRhoDistanceWithTheBestPV function */
  struct DeltaRhoToVertex : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "DeltaRhoToVertex"; }

    template <typename VContainer, typename Particle>
      auto operator()( VContainer const& vertices, Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return ( endVertexPos( composite ) - endVertexPos( Sel::getBestPV( composite, vertices ) ) ).Rho();
    }
  };

  /** BPVDIRA */
  struct CosDirectionAngleToVertex : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "CosDirectionAngleToVertex"; }

    template <typename VContainer, typename Composite>
      auto operator()( VContainer const& vertices, Composite const& composite ) const {
      // Calculate the angle between the momentum vector of the composite and
      // the vector connecting its decay vertex to the origin
      using Sel::Utils::endVertexPos;
      using Sel::Utils::threeMomentum;
      using std::sqrt;
      auto const flight = endVertexPos( composite ) - endVertexPos( Sel::getBestPV( composite, vertices ) );
      auto const mom    = threeMomentum( composite );
      return dot( mom, flight ) / sqrt( flight.mag2() * mom.mag2() );
    }
  };

  struct PseudoRapidityFromVertex : public Function {
    static constexpr auto name() { return "PseudoRapidityFromVertex"; }

    template <typename VContainer, typename Composite>
      auto operator()( VContainer const& vertices, Composite const& composite ) const {
      using Sel::Utils::endVertexPos;
      return ( endVertexPos( composite ) - endVertexPos( Sel::getBestPV( composite, vertices ) ) ).eta();
    }
  };

  struct CorrectedMass : public Function {
    static constexpr auto name() { return "CorrectedMass"; }

    template <typename VContainer, typename Composite>
      auto operator()( VContainer const& vertices, Composite const& composite ) const {

      using Sel::Utils::endVertexPos;
      using Sel::Utils::mass2;
      using Sel::Utils::threeMomentum;
      using std::sqrt;

      // TODO this is lifted from LoKi but some checks were dropped.

      // Calculate the corrected mass using the 4-momentum (p4) and the
      // flight vector:
      auto const d = endVertexPos( composite ) - endVertexPos( Sel::getBestPV( composite, vertices ) );

      // Get the pT variable that we need
      auto const mom  = threeMomentum( composite );
      auto const perp = mom - d * ( dot( mom, d ) / d.mag2() );
      auto const pt   = perp.mag();

      // Calculate the corrected mass
      return pt + sqrt( mass2( composite ) + pt * pt );
    }
  };

  struct Lifetime : public Function {
    static constexpr auto name() { return "Lifetime"; }

    void bind( TopLevelInfo& top_level ) {
      m_lifetime_calc = Functors::detail::DefaultLifetimeFitter_t( top_level.algorithm() );
    }

    template <typename VContainer, typename Composite>
      auto operator()( VContainer const& vertices, Composite const& composite ) const {
      auto const& bestPV = Sel::getBestPV( composite, vertices );
      return std::get<0>( m_lifetime_calc.Lifetime( bestPV, composite ) );
    }

  private:
    Functors::detail::DefaultLifetimeFitter_t m_lifetime_calc;
  };

  struct LifetimeChi2 : public Function {
    static constexpr auto name() { return "LifetimeChi2"; }

    template <typename VContainer, typename Composite>
      auto operator()( VContainer const& vertices, Composite const& composite ) const {
      auto const& bestPV = Sel::getBestPV( composite, vertices );
      return std::get<2>( m_lifetime_calc.Lifetime( bestPV, composite ) );
    }

  private:
    Functors::detail::DefaultLifetimeFitter_t m_lifetime_calc;
  };

  struct ComputeDecayLengthSignificance : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "ComputeDecayLengthSignificance"; }

    void bind( TopLevelInfo& top_level ) {
      m_lifetime_calc = Functors::detail::DefaultLifetimeFitter_t( top_level.algorithm() );
    }

    template <typename VContainer, typename Composite>
      auto operator()( VContainer const& vertices, Composite const& composite ) const {
      auto const& bestPV = Sel::getBestPV( composite, vertices );
      return m_lifetime_calc.DecayLengthSignificance( bestPV, composite );
    }

  private:
    Functors::detail::DefaultLifetimeFitter_t m_lifetime_calc;
  };

} // namespace Functors::detail

namespace Functors::Composite {
  /** @brief Flight distance chi2 to the "best" one of the given vertices.
   *
   *  Note that if the given data object contains a vertex link then that will
   *  be checked for compatibility with the given vertex container and, if it
				     *  matches, be used.
   */

  template <int N, typename VContainer = detail::DefaultPVContainer_t>
    auto MotherTrajectoryDistanceOfClosestApproachChi2( std::integral_constant<int, N>, std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::MotherTrajectoryDistanceOfClosestApproachChi2<N>, VContainer>{
      std::move( vertex_location )};
  }
  
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto FlightDistanceChi2ToVertex( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::FlightDistanceChi2ToVertex, VContainer>{
      std::move( vertex_location )};
  }

  /** @brief Z component of flight distance of the given composite. */
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto DeltaZToVertex( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::DeltaZToVertex, VContainer>{std::move( vertex_location )};
  }

  /** @brief Rho component of flight distance of the given composite. */
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto DeltaRhoToVertex( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::DeltaRhoToVertex, VContainer>{std::move( vertex_location )};
  }
  /** @brief Calculate the cosine of the angle between the momentum vector of
   *         the composite and the vector connecting its decay vertex to the
   *         "best" one of the given vertices.
   *
   *  Note that if the given data object contains a vertex link then that will
   *  be checked for compatibility with the given vertex container and, if it
   *  matches, be used.
   */
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto CosDirectionAngleToVertex( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::CosDirectionAngleToVertex, VContainer>{
      std::move( vertex_location )};
  }

  /** @brief Pseudorapidity of the flight vector of the given composite. */
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto PseudoRapidityFromVertex( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::PseudoRapidityFromVertex, VContainer>{std::move( vertex_location )};
  }

  /** @brief Calculate the corrected mass. */
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto CorrectedMass( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::CorrectedMass, VContainer>{std::move( vertex_location )};
  }

  /** @brief Calculate the composite mass using the given child mass hypotheses. */
  template <typename... MassInputs>
    struct MassWithHypotheses : public Function {
    /** Create a mass functor with a list of mass hypotheses that can be a mix
     *  of floating point values (in MeV) and names of particles. Particle
     *  names are translated into mass values using the ParticlePropertySvc.
     */
  MassWithHypotheses( std::tuple<MassInputs...> mass_inputs ) : m_mass_inputs{std::move( mass_inputs )} {}

    void bind( TopLevelInfo& top_level ) { bind( top_level, std::index_sequence_for<MassInputs...>{} ); }

    template <typename CombinationType>
      auto operator()( CombinationType const& combination ) const {
      // Calculate the mass from the child 3-momenta and the given
      // mass hypotheses. Start by checking we have the correct number.
      using Sel::Utils::decayProducts;
      auto children    = decayProducts( combination );
      auto NumChildren = children.size();
      if ( sizeof...( MassInputs ) != NumChildren || m_mass_values.size() != NumChildren ) {
        throw GaudiException{
	  "Mismatch between number of mass values given (" + std::to_string( sizeof...( MassInputs ) ) +
	    ") and the number of children in the given object (" + std::to_string( NumChildren ) + ")",
            "Functors::Composite::Mass", StatusCode::FAILURE};
      }
      using Sel::Utils::threeMomentum;
      using std::sqrt;
      using float_t = decltype( threeMomentum( Sel::Utils::deref_if_ptr( children[0] ) ).mag2() );
      float_t E{0};
      for ( const auto& [i, child] : LHCb::range::enumerate( children ) ) {
        E += sqrt( threeMomentum( Sel::Utils::deref_if_ptr( child ) ).mag2() + m_mass_values[i] * m_mass_values[i] );
      }
      return sqrt( E * E - threeMomentum( combination ).mag2() );
    }

  private:
    template <std::size_t... Ns>
      void bind( TopLevelInfo& top_level, std::index_sequence<Ns...> ) {
      // Avoid setting up the service if we don't need it (e.g. all hypotheses)
      // were specified numerically)
      std::optional<ServiceHandle<LHCb::IParticlePropertySvc>> pp_svc{std::nullopt};
      // Helper to convert a member of m_mass_inputs to a numeric value.
      auto const converter = [&]( auto const& mass_or_name ) {
        if constexpr ( std::is_same_v<float, std::decay_t<decltype( mass_or_name )>> ) {
	    return mass_or_name;
	  } else {
          if ( !pp_svc ) {
            pp_svc.emplace( top_level.algorithm(), top_level.generate_property_name(), "LHCb::ParticlePropertySvc" );
          }
          auto const* pp = pp_svc.value()->find( mass_or_name );
          if ( !pp ) {
            throw GaudiException{"Couldn't get ParticleProperty for particle '" + mass_or_name + "'",
		"Functors::Composite::Mass::bind()", StatusCode::FAILURE};
          }
          return pp->mass();
        }
      };
      ( ( m_mass_values[Ns] = converter( std::get<Ns>( m_mass_inputs ) ) ), ... );
    }

    std::tuple<MassInputs...>                  m_mass_inputs;
    std::array<float, sizeof...( MassInputs )> m_mass_values{};
  };

  /** @brief Return the input object's mass as defined by an accessor. */
  struct Mass : public Function {
    static constexpr auto name() { return "Mass"; }

    template <typename Particle>
      auto operator()( Particle const& particle ) const {
      using Sel::Utils::mass2;
      using std::sqrt;
      return sqrt( mass2( Sel::Utils::deref_if_ptr( particle ) ) );
    }
  };

  /** VX */
  /** @brief Calculate Vertex X position in the detector. NOT using PV as reference. */
  struct EndVertexX : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "EndVertexX"; }

    template <typename Particle>
      auto operator()( Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return endVertexPos( Sel::Utils::deref_if_ptr( composite ) ).X();
    }
  };

  /** VY */
  /** @brief Calculate Vertex Y position in the detector. NOT using PV as reference. */
  struct EndVertexY : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "EndVertexY"; }

    template <typename Particle>
      auto operator()( Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return endVertexPos( Sel::Utils::deref_if_ptr( composite ) ).Y();
    }
  };

  /** VZ */
  /** @brief Calculate Vertex Z position in the detector. NOT using PV as reference. */
  struct EndVertexZ : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "EndVertexZ"; }

    template <typename Particle>
      auto operator()( Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return endVertexPos( Sel::Utils::deref_if_ptr( composite ) ).Z();
    }
  };

  /** VRho */
  /** @brief Calculate Vertex Rho = sqrt(x*x+y*y) position in the detector. NOT using PV as reference. */
  struct EndVertexRho : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "EndVertexRho"; }

    template <typename Particle>
      auto operator()( Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return endVertexPos( Sel::Utils::deref_if_ptr( composite ) ).rho();
    }
  };

  /** @brief Calculate the lifetime of the particle with respect to the PV. */
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto Lifetime( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::Lifetime, VContainer>{std::move( vertex_location )};
  }

  /** @brief Calculate the lifetime chi^2 of the particle with respect to the PV. */
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto LifetimeChi2( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::LifetimeChi2, VContainer>{std::move( vertex_location )};
  }

  /** @brief Calculate the lifetime of the particle with respect to the PV. */
  template <typename VContainer = detail::DefaultPVContainer_t>
    auto ComputeDecayLengthSignificance( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::ComputeDecayLengthSignificance, VContainer>{
      std::move( vertex_location )};
  }
} // namespace Functors::Composite
