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
#include "Functors/Function.h"
#include "Functors/Utilities.h"

/** @file  Combination.h
 *  @brief Definitions of functors for tuple-likes of track-like objects.
 *
 *  Typical argument types for these functors are Sel::ParticleCombination and
 *  Sel::ParticleCombinationSpan.
 */

namespace Functors::detail {
  template <typename DistanceCalculator, typename T, T N, T M>
  struct DistanceOfClosestApproach : public Function {
    static_assert( N >= 1 && M >= 1, "Indices start from 1 for LoKi compatibility." );
    DistanceOfClosestApproach( std::integral_constant<T, N>, std::integral_constant<T, M> ) {}

    void bind( TopLevelInfo& top_level ) { m_dist_calc.emplace( top_level.algorithm() ); }

    template <typename CombinationType>
    auto operator()( CombinationType const& combination ) const {
      return combination.template doca<N - 1, M - 1>( *m_dist_calc );
    }

  private:
    std::optional<DistanceCalculator> m_dist_calc;
  };

  template <typename DistanceCalculator, typename T, T N, T M>
    struct DistanceOfClosestApproachChi2 : public Function {
    static_assert( N >= 1 && M >= 1, "Indices start from 1 for LoKi compatibility." );
    DistanceOfClosestApproachChi2( std::integral_constant<T, N>, std::integral_constant<T, M> ) {}
    
    void bind( TopLevelInfo& top_level ) { m_dist_calc.emplace( top_level.algorithm() ); }
    
    template <typename CombinationType>
      auto operator()( CombinationType const& combination ) const {
      return combination.template docachi2<N - 1, M - 1>( *m_dist_calc );
    }
    
  private:
    std::optional<DistanceCalculator> m_dist_calc;
  };
  
  template <typename DistanceCalculator, typename T, T N>
    struct MotherTrajectoryDistanceOfClosestApproachChi2 : public Function {
    static_assert( N >= 1, "Indices start from 1 for LoKi compatibility." );
    MotherTrajectoryDistanceOfClosestApproachChi2( std::integral_constant<T, N> ) {}
    
    void bind( TopLevelInfo& top_level ) { m_dist_calc.emplace( top_level.algorithm() ); }
    
    template <typename CombinationType>
      auto operator()( CombinationType const& combination ) const {
      return combination.template mtdocachi2<N - 1>( *m_dist_calc );
    }
    
  private:
    std::optional<DistanceCalculator> m_dist_calc;
  };
 
  //  template <typename DistanceCalculator, typename T, T N, T M>
  template <typename T, T N, T M>
  struct CosAngleBetweenDecayProducts : public Function {
    static_assert( N >= 1 && M >= 1, "Indices start from 1 for LoKi compatibility." );
    CosAngleBetweenDecayProducts( std::integral_constant<T, N>, std::integral_constant<T, M> ) {}

    template <typename CombinationType>
    auto operator()( CombinationType const& combination ) const {
      return combination.template cos_angle_prod<N - 1, M - 1>();
    }
  };
} // namespace Functors::detail

/** @namespace Functors::Combination
 *
 * TODO
 */
namespace Functors::Combination {
  template <typename DistanceCalculator = detail::DefaultDistanceCalculator_t, typename T, T N, T M>
    auto DistanceOfClosestApproach( std::integral_constant<T, N>, std::integral_constant<T, M> ) {
    // This needs a wrapper function so that `T`, `N` and `M` can be deduced
    return detail::DistanceOfClosestApproach<DistanceCalculator, T, N, M>( {}, {} );
  }

  template <typename DistanceCalculator = detail::DefaultDistanceCalculator_t, typename T, T N, T M>
    auto DistanceOfClosestApproachChi2( std::integral_constant<T, N>, std::integral_constant<T, M> ) {
    // This needs a wrapper function so that `T`, `N` and `M` can be deduced
    return detail::DistanceOfClosestApproachChi2<DistanceCalculator, T, N, M>( {}, {} );
  }
  
  template <typename DistanceCalculator = detail::DefaultDistanceCalculator_t, typename T, T N, T M>
    auto MotherTrajectoryDistanceOfClosestApproachChi2( std::integral_constant<T, N>, std::integral_constant<T, M> ) {
    // This needs a wrapper function so that `T`, `N` and `M` can be deduced
    return detail::MotherTrajectoryDistanceOfClosestApproachChi2<DistanceCalculator, T, N, M>( {}, {} );
  }
  
  template <typename T, T N, T M>
    auto CosAngleBetweenDecayProducts( std::integral_constant<T, N>, std::integral_constant<T, M> ) {
    // This needs a wrapper function so that `T`, `N` and `M` can be deduced
    return detail::CosAngleBetweenDecayProducts<T, N, M>( {}, {} );
  }

  template <typename DistanceCalculator = detail::DefaultDistanceCalculator_t>
  struct MaxDistanceOfClosestApproach : public Function {
    void bind( TopLevelInfo& top_level ) { m_dist_calc.emplace( top_level.algorithm() ); }
    template <typename CombinationType>
    auto operator()( CombinationType const& combination ) const {
      return combination.maxdoca( *m_dist_calc );
    }

  private:
    std::optional<DistanceCalculator> m_dist_calc;
  };

  template <typename DistanceCalculator = detail::DefaultDistanceCalculator_t>
  struct MaxDistanceOfClosestApproachChi2 : public Function {
    void bind( TopLevelInfo& top_level ) { m_dist_calc.emplace( top_level.algorithm() ); }
    template <typename CombinationType>
    auto operator()( CombinationType const& combination ) const {
      return combination.maxdocachi2( *m_dist_calc );
    }

  private:
    std::optional<DistanceCalculator> m_dist_calc;
  };

  template <typename DistanceCalculator = detail::DefaultDistanceCalculator_t>
  struct MaxDistanceOfClosestApproachCut : public Predicate {
    float m_thresh;
    MaxDistanceOfClosestApproachCut( float thresh ) : m_thresh( thresh ) {}

    void bind( TopLevelInfo& top_level ) { m_dist_calc.emplace( top_level.algorithm() ); }

    template <typename CombinationType>
    auto operator()( CombinationType const& combination ) const {
      return combination.maxdocacut( *m_dist_calc, m_thresh );
    }

  private:
    std::optional<DistanceCalculator> m_dist_calc;
  };

  template <typename DistanceCalculator = detail::DefaultDistanceCalculator_t>
  struct MaxDistanceOfClosestApproachChi2Cut : public Predicate {
    float m_thresh{0.f};
    MaxDistanceOfClosestApproachChi2Cut( float thresh ) : m_thresh{thresh} {}

    void bind( TopLevelInfo& top_level ) { m_dist_calc.emplace( top_level.algorithm() ); }

    template <typename CombinationType>
    auto operator()( CombinationType const& combination ) const {
      return combination.maxdocachi2cut( *m_dist_calc, m_thresh );
    }

  private:
    std::optional<DistanceCalculator> m_dist_calc;
  };

  struct Charge : public Function {
    template <typename CombinationType>
    auto operator()( CombinationType const& combination ) const {
      return combination.charge();
    }
  };
} // namespace Functors::Combination
