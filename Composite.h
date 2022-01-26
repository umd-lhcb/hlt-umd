

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
#include "Event/State.h"
#include "Functors/Function.h"
#include "Functors/Utilities.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ParticleProperty.h"
#include "LHCbMath/MatVec.h"
#include "SelKernel/ParticleTraits.h"
#include "SelKernel/Utilities.h"
#include "SelKernel/VectorOps.h"
#include "SelKernel/VertexRelation.h"
#include "TrackKernel/TrackCompactVertex.h"

#include "GaudiKernel/detected.h"

#include <cassert>

//let's try
#include "SelTools/DistanceCalculator.h"
#include "Gaudi/Algorithm.h"
#include "SelTools/Utilities.h"                                                              

//added to hardcode DistanceCalculator, will delete
/*#include "SelKernel/State4.h"
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
*/

/** @file  Composite.h
 *  @brief Definitions of functors for composite-particle-like objects.
 */

/** @namespace Functors::Composite
 *
 *  Functors that make sense for composite particles (i.e. ones that have a
 *  vertex as well as a trajectory)
 */
namespace Functors::detail {
  template <typename T>
  using has_threeMomentum_ = decltype( std::declval<T>().threeMomentum() );
  template <typename T>
  inline constexpr bool has_threeMomentum_v = Gaudi::cpp17::is_detected_v<has_threeMomentum_, T>;

  template <typename>
  struct is_array_of_v2_RecVertex_pointers : std::false_type {};

  template <std::size_t N>
  struct is_array_of_v2_RecVertex_pointers<std::array<LHCb::Event::v2::RecVertex const*, N>> : std::true_type {};

  template <typename T>
  inline constexpr bool is_array_of_v2_RecVertex_pointers_v = is_array_of_v2_RecVertex_pointers<T>::value;

  template <typename T>
  using has_static_NumChildren_ = decltype( T::NumChildren );
  template <typename T>
  inline constexpr bool has_static_NumChildren_v = Gaudi::cpp17::is_detected_v<has_static_NumChildren_, T>;

  /** Get the position of the best associated PV
   */
  template <typename Composite, typename Vertices>
  decltype( auto ) getBestPVPosition( Composite const& composite, Vertices const& vertices ) {
    // Get the best PV, which would ideally be a nice proxy thing with some
    // gather operations under the hood, but which may actually be a vector
    // of RecVertex pointers...
    auto const& bestPV = Sel::getBestPV( composite, vertices );
    using bestPV_t     = std::decay_t<decltype( bestPV )>;
    if constexpr ( detail::is_array_of_v2_RecVertex_pointers_v<bestPV_t> ) {
      // Assume that 'composite' is a new-style proxy object
      using float_v = typename Composite::dType::float_v;
      std::array<float, float_v::size()> x, y, z;
      std::size_t const                  num_valid = popcount( composite.loop_mask() );
      for ( auto i = 0ul; i < num_valid; ++i ) {
        x[i] = bestPV[i]->position().x();
        y[i] = bestPV[i]->position().y();
        z[i] = bestPV[i]->position().z();
      }
      for ( auto i = num_valid; i < x.size(); ++i ) { x[i] = y[i] = z[i] = std::numeric_limits<float>::lowest(); }
      return LHCb::LinAlg::Vec3{float_v{x.data()}, float_v{y.data()}, float_v{z.data()}};
    } else {
      // If only...
      using Sel::Utils::endVertexPos;
      return endVertexPos( bestPV );
    }
  }


  struct FlightDistanceChi2ToVertex : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "FlightDistanceChi2ToVertex"; }

    template <typename VContainer, typename Particle>
      auto operator()( VContainer const& vertices, Particle const& composite ) const {
      // Get the associated PV -- this uses a link if it's available and
      // computes the association if it's not.
      auto const& bestPV = Sel::getBestPV( composite, vertices );

      // Now calculate the flight distance chi2 between 'composite' and 'bestPV'
      return Sel::Utils::flightDistanceChi2( composite, bestPV );
    }
  };  

  /**MTDOCACHI2**/
  struct MotherTrajectoryDistanceOfClosestApproachChi2 : public Function {
    void bind( TopLevelInfo& top_level ){ 
      dist_calc = Sel::DistanceCalculator(top_level.algorithm() );
    }
    template <typename VContainer, typename Particle> 
      auto operator()( VContainer const& vertices, Particle const& composite ) const {
      auto const& bestPV = Sel::getBestPV(composite, vertices);
      if constexpr (Sel::Utils::is_legacy_particle<Particle>) {
	  const auto& children = composite.daughtersVector();
	  const auto& pN = children[0];
	  const auto& tempMother = composite.clone();
	  tempMother->setReferencePoint( bestPV.position() );
	  tempMother->setPosCovMatrix( bestPV.covMatrix() );
	  return dist_calc.particleDOCAChi2(*pN, *tempMother)
	}
      else {
	return 1.0;
      }
    }
  private:
    Sel::DistanceCalculator dist_calc;
  };
  
  /** BPVVDZ */
  /** origin version in Loki defined in PHYS/PHYS_v25r1/Phys/LoKiPhys/src/Particles20.cpp,
  VertexZDistanceWithTheBestPV function */
  //template <typename T, T N>
  struct DeltaZToVertex : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "DeltaZToVertex"; }

    template <typename VContainer, typename Particle>
    auto operator()( VContainer const& vertices, Particle const& composite ) const {
      // Get the position of associated PV
      // Get the position of end vertex
      using Sel::Utils::endVertexPos;
      return ( endVertexPos( composite ) - getBestPVPosition( composite, vertices ) ).Z();
      //return 1.0;
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
      // Get the position of associated PV
      // Get the position of end vertex
      using Sel::Utils::endVertexPos;
      return ( endVertexPos( composite ) - getBestPVPosition( composite, vertices ) ).Rho();
    }
  };

  /** BPVDIRA */
  struct CosDirectionAngleToVertex : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "CosDirectionAngleToVertex"; }

    template <typename VContainer, typename Composite>
    auto operator()( VContainer const& vertices, Composite const& composite ) const {
      // Get the associated PV -- this uses a link if it's available and
      // computes the association if it's not.

      // Calculate the angle between the momentum vector of the composite and
      // the vector connecting its decay vertex to the origin
      using Sel::Utils::endVertexPos;
      using Sel::Utils::threeMomentum;
      auto const& mom             = threeMomentum( composite );
      auto const& bestPV_position = getBestPVPosition( composite, vertices );
      auto const& decay_vertex_v2 = endVertexPos( composite );
      // The flight direction vector
      auto const flight = decay_vertex_v2 - bestPV_position;
      using std::sqrt;
      return dot( mom, flight ) / sqrt( flight.mag2() * mom.mag2() );
    }
  };

  struct PseudoRapidityFromVertex : public Function {
    static constexpr auto name() { return "PseudoRapidityFromVertex"; }

    template <typename VContainer, typename Composite>
    auto operator()( VContainer const& vertices, Composite const& composite ) const {
      using Sel::Utils::endVertexPos;
      auto const& bestPV_position = getBestPVPosition( composite, vertices );
      auto        flight          = endVertexPos( composite ) - bestPV_position;
      return flight.eta();
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

      auto const decay_vertex_position = endVertexPos( composite );

      // Get the associated [primary] vertex position
      auto const bestPV_position = getBestPVPosition( composite, vertices );

      // Three-momentum
      auto mom = threeMomentum( composite );

      // TODO this is lifted from LoKi but some checks were dropped.

      // Calculate the corrected mass using the 4-momentum (p4) and the
      // flight vector:
      auto const d = decay_vertex_position - bestPV_position;

      // Get the pT variable that we need
      auto const dmag2 = d.mag2();
      auto const perp  = mom - d * ( dot( mom, d ) / dmag2 );
      auto const pt    = perp.mag();

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
      auto const& bestPV  = Sel::getBestPV( composite, vertices );
      auto        results = m_lifetime_calc.Lifetime( bestPV, composite );
      return std::get<0>( results );
    }

  private:
    Functors::detail::DefaultLifetimeFitter_t m_lifetime_calc;
  };

  struct LifetimeChi2 : public Function {
    static constexpr auto name() { return "LifetimeChi2"; }

    template <typename VContainer, typename Composite>
    auto operator()( VContainer const& vertices, Composite const& composite ) const {
      auto const& bestPV  = Sel::getBestPV( composite, vertices );
      auto        results = m_lifetime_calc.Lifetime( bestPV, composite );
      return std::get<2>( results );
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
  template <typename VContainer = detail::DefaultPVContainer_t>
  auto FlightDistanceChi2ToVertex( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::FlightDistanceChi2ToVertex, VContainer>{
      std::move( vertex_location )};
  }

  template <typename VContainer = detail::DefaultPVContainer_t>
    auto MotherTrajectoryDistanceOfClosestApproachChi2( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::MotherTrajectoryDistanceOfClosestApproachChi2, VContainer>{std::move( vertex_location )};
  }
  

  /** @brief Z component of flight distance of the given composite. */
  template < typename VContainer = detail::DefaultPVContainer_t >
  //template <
  auto DeltaZToVertex( std::string vertex_location ) {
    return detail::DataDepWrapper<Function, detail::DeltaZToVertex, VContainer>
    {std::move( vertex_location )};
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
      int NumChildren;
      // If its a v1 Particle
      if constexpr ( Sel::Utils::is_legacy_particle<CombinationType> ) {
        combination.daughters().size();
      }
      // If it's anything else
      else {
        NumChildren = combination.numChildren();
      }
      if ( sizeof...( MassInputs ) != NumChildren ) {
        throw GaudiException{
            "Mismatch between number of mass values given (" + std::to_string( sizeof...( MassInputs ) ) +
                ") and the number of children in the given object (" + std::to_string( NumChildren ) + ")",
            "Functors::Composite::Mass", StatusCode::FAILURE};
      }
      float E{0.f}; // TODO make SIMD-friendly?
      using Sel::Utils::threeMomentum;
      using std::sqrt;
      std::size_t idau = 0;
      // If it's a v1 Particle
      if constexpr ( Sel::Utils::is_legacy_particle<CombinationType> ) {
        for ( std::size_t idau = 0; idau < m_mass_values.size(); ++idau ) {
          E += sqrt( threeMomentum( *combination.daughters()[idau] ).mag2() +
                     m_mass_values[idau] * m_mass_values[idau] );
        }
        return sqrt( E * E - threeMomentum( combination ).mag2() );
      }
      // Is TrackCompactVertex
      else if constexpr ( Sel::Utils::is_trackcompactvertex<CombinationType>::value ) {
        for ( idau = 0; idau < m_mass_values.size(); ++idau ) {
          E += sqrt( std::pow( combination.daughters.P( idau ), 2 ) + m_mass_values[idau] * m_mass_values[idau] );
        }
        return sqrt( E * E - threeMomentum( combination ).mag2() );
      }
      // If it's a Combiner::ParticleCombination of v1 particles
      else if constexpr ( Sel::Utils::is_particlecombination<CombinationType>::value ) {
        throw std::invalid_argument(
            "Particle combinations not supported in MassWithHypotheses. Use a composite or vertex instead." );
      } else {
        throw std::invalid_argument( "v2::Particles cannot currently support accessing daughter information." );
        /*TODO This will break for v2::Particles, and there isn't functionality yet to support it!
        Update once we have a way to retrieve v2::Particle descendent info.
        Once we can access the children the code will look something like:
        */
      }
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

      if constexpr ( Sel::Utils::is_legacy_particle<Particle> ) {
        return Sel::Utils::deref_if_ptr( particle ).momentum().M();
      } else {
        return particle.mass();
      }
    }
  };

  /** VX */
  /** @brief Calculate Vertex X position in the detector. NOT using PV as reference. */
  struct VertexX : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "VertexX"; }

    template <typename Particle>
    auto operator()( Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return endVertexPos( composite ).X();
    }
  };

  /** VY */
  /** @brief Calculate Vertex Y position in the detector. NOT using PV as reference. */
  struct VertexY : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "VertexY"; }

    template <typename Particle>
    auto operator()( Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return endVertexPos( composite ).Y();
    }
  };

  /** VZ */
  /** @brief Calculate Vertex Z position in the detector. NOT using PV as reference. */
  struct VertexZ : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "VertexZ"; }

    template <typename Particle>
    auto operator()( Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return endVertexPos( composite ).Z();
    }
  };

  /** VRho */
  /** @brief Calculate Vertex Rho = sqrt(x*x+y*y) position in the detector. NOT using PV as reference. */
  struct VertexRho : public Function {
    /** This allows some error messages to be customised. It is not critical. */
    static constexpr auto name() { return "VertexRho"; }

    template <typename Particle>
    auto operator()( Particle const& composite ) const {
      using Sel::Utils::endVertexPos;
      return endVertexPos( composite ).rho();
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
