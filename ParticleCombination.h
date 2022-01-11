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
#include "Kernel/HeaderMapping.h"
#include "SelKernel/ParticleAccessors.h"
#include "SelKernel/Utilities.h"
#include "SelKernel/VectorOps.h"
#include "TrackKernel/TrackVertexUtils.h"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>

#include <functional>

//LoKi
/*#include "LoKi/Child.h"
#include "LoKi/Constants.h"
#include "LoKi/ILoKiSvc.h"
#include "LoKi/Interface.h"
#include "LoKi/Particles26.h"
#include "LoKi/GetTools.h"*/

namespace Sel {
  // Forward declarations so that shorthand versions of these functions can be
  // defined in Sel::detail::ParticleCombinationBase
  template <typename...>
  struct ParticleCombination;

  template <std::size_t..., typename collection_t, typename transform_t>
  auto transform( collection_t const& comb, transform_t transform );
  template <std::size_t..., typename... child_ts, typename transform_t>
  auto transform( ParticleCombination<child_ts...> const& comb, transform_t transform );

  template <typename collection_t, typename transform_t, typename reduce_t>
  auto transform_reduce( collection_t const& comb, transform_t transform, reduce_t reduce );
  template <typename... child_ts, typename transform_t, typename reduce_t>
  auto transform_reduce( ParticleCombination<child_ts...> const& comb, transform_t transform, reduce_t reduce );

  template <typename collection_t, typename transform_t, typename reduce_t>
  auto pairwise_transform_reduce( collection_t const& comb, transform_t transform, reduce_t reduce );
  template <typename... child_ts, typename transform_t, typename reduce_t>
  auto pairwise_transform_reduce( ParticleCombination<child_ts...> const& comb, transform_t transform,
                                  reduce_t reduce );
//so i think pairwise_transform_reduce extracts all the particles without input, transform
//takes ints and transforms to the corresponding particles, and transform_reduce takes no arg
//uments and extracts a single(?) particle

//forward declaration of child_t so i can use it for mtdocachi2
  /** Get a reference to the Nth child.
   *  If this is std::reference_wrapper<T>, T& is returned.
   */
  template <std::size_t N>
    child_t<N>& get() {
    return std::get<N>( m_children );
  }

  /** Get a const reference to the Nth child.
   *  If this is std::reference_wrapper<T [const]>, T const& is returned.
   */
  template <std::size_t N>
    child_t<N> const& get() const {
    return std::get<N>( m_children );
  }
 
} // namespace Sel

namespace Sel::detail {
  /** @class ParticleCombinationBase
   *  CRTP base class implementing operations that can be applied to a
   *  collection of particles. Only derived classes should be used
   *  directly.
   */
  template <typename Derived>
    struct ParticleCombinationBase {
      /** Calculate the distance of closest approach between child `N` and child
       *  `M` of the combination.
       */
      template <std::size_t N, std::size_t M, typename DistanceCalculator>
	auto doca( DistanceCalculator const& dist_calc ) const {
	return transform<N, M>( [&]( auto const& pN, auto const& pM ) { return dist_calc.particleDOCA( pN, pM ); } );
      }
      
      template <typename DistanceCalculator>
      auto maxdoca( DistanceCalculator const& dist_calc ) const {
	return pairwise_transform_reduce(
					 [&]( auto const& pA, auto const& pB ) { return dist_calc.particleDOCA( pA, pB ); }, max{} );
      }
      
      /** Calculate the distance of closest approach chi2 between child `N` and
       *  child `M` of the combination.
       */
      template <std::size_t N, std::size_t M, typename DistanceCalculator>
	auto docachi2( DistanceCalculator const& dist_calc ) const {
	return transform<N, M>( [&]( auto const& pN, auto const& pM ) { return dist_calc.particleDOCAChi2( pN, pM ); } );
      }
      
      template <std::size_t N, typename DistanceCalculator>
	auto mtdocachi2( DistanceCalculator const& dist_calc ) const {
	std::unique_ptr<LHCb::Particle> pN = child_t<N>;
	if constexpr (Sel::Utils::is_legacy_particle<pN>) {
	    LoKi::Particles::MTDOCA::result_type LoKi::Particles::MTDOCA::operator()(LoKi::Particles::MTDOCA::argument pMother ) const {
	      if ( 0 == pMother ) {
		error("Invalid particle, return 'InvalidDistance'").ignore();
		return LoKi::Constants::InvalidDistance;
	      }
	      const LHCb::Particle* pChild = LoKi::Child::child( pMother, getIndex() );
	    if ( 0 == pChild ) {
	      Error( "Invalid child particle, return 'InvalidDistance'").ignore();
	      return LoKi::Constants::InvalidDistance;
	    }
	    
	    // clone the mother and move it to the PV
	    std::unique_ptr<LHCb::Particle> tempMother( pMother->clone() );
	    //pChild is child, pMother is mother
	    const LHCb::VertexBase* aPV = bestVertex( pChild );
	    
	    //Update the position errors
	    tempMother->setReferencePoint( aPV->position() );
	    tempMother->setPosCovMatrix( aPV->covMatrix() );
	    
	    tempMother = moveMother(pMother, pN);
	    /** Calculate the distance of closest approach chi2 between child `N` and
	     *  child `M` of the combination.                                              
	     */
	    return transform<N>( [&]( auto const& pN) { return dist_calc.particleDOCAChi2( pN, tempMother ); } ); 
	    }
	  }
      }
      
    template <typename DistanceCalculator>
    auto maxdocachi2( DistanceCalculator const& dist_calc ) const {
      return pairwise_transform_reduce(
          [&]( auto const& pA, auto const& pB ) { return dist_calc.particleDOCAChi2( pA, pB ); }, max{} );
    }

    // TODO benchmark if short-circuiting actually helps
    template <typename DistanceCalculator>
    auto maxdocacut( DistanceCalculator const& dist_calc, float cut ) const {
      return maxdoca( dist_calc ) < cut;
    }
    template <typename DistanceCalculator>
    auto maxdocachi2cut( DistanceCalculator const& dist_calc, float cut ) const {
      return maxdocachi2( dist_calc ) < cut;
    }

    /** Sum the 3-momenta of the contained particles.
     */
    auto threeMomentum() const {
      return transform_reduce( []( auto const& p ) { return Sel::get::threeMomentum( p ); }, std::plus<>{} );
    }

    /** Sum the 4-momenta of the contained particles.
     */
    auto momentum() const {
      return transform_reduce( []( auto const& p ) { return Sel::get::momentum( p ); }, std::plus<>{} );
    }

    /** Get the transverse part of the vector sum of 3-momenta of the contained
     *  particles.
     */
    auto pt() const { return Sel::get::pt( this->threeMomentum() ); }

    auto charge() const {
      return transform_reduce( []( auto const& p ) { return p.charge(); }, std::plus<>{} );
    }

    auto mass() const { return Sel::get::mass( this->momentum() ); }

    /** Unfortunately this is needed in order for explicit mask functionality
     *  of the functors to work. Returns the logical AND of the loop_mask()
     *  results from all children.
     */
    auto loop_mask() const {
      return transform_reduce( []( auto const& p ) { return p.loop_mask(); }, std::logical_and<>{} );
    }

  private:
    template <std::size_t... Is, typename transform_t>
    auto transform( transform_t transform_fn ) const {
      return ::Sel::transform<Is...>( static_cast<Derived const&>( *this ), std::move( transform_fn ) );
    }

    template <typename transform_t, typename reduce_t>
    auto transform_reduce( transform_t transform_fn, reduce_t reduce ) const {
      return ::Sel::transform_reduce( static_cast<Derived const&>( *this ), std::move( transform_fn ),
                                      std::move( reduce ) );
    }

    template <typename pairwise_transform_t, typename reduce_t>
    auto pairwise_transform_reduce( pairwise_transform_t transform_fn, reduce_t reduce ) const {
      return ::Sel::pairwise_transform_reduce( static_cast<Derived const&>( *this ), std::move( transform_fn ),
                                               std::move( reduce ) );
    }

    struct max {
      template <typename... Args>
      auto operator()( Args&&... args ) const {
        using std::max;
        return max( std::forward<Args>( args )... );
      }
    };
  };

  template <typename T>
  struct remove_reference_wrapper {
    using type = T;
  };

  template <typename T>
  struct remove_reference_wrapper<std::reference_wrapper<T>> {
    using type = T;
  };

  template <typename T>
  using remove_reference_wrapper_t = typename remove_reference_wrapper<T>::type;

  template <typename T>
  using as_const_ref_t = typename std::decay_t<remove_reference_wrapper_t<std::decay_t<T>>> const&;

  template <typename... Ts>
  using first_t = std::tuple_element_t<0, std::tuple<Ts...>>;

  template <typename>
  struct mask_true_helper {};
  template <typename... Ts>
  struct mask_true_helper<std::tuple<Ts...>> {
    static constexpr auto _() {
      static_assert( ( Sel::Utils::has_static_mask_true_v<Ts> && ... ) );
      using mask_true_t = std::common_type_t<decltype( Ts::mask_true() )...>;
      return mask_true_t{first_t<Ts...>::mask_true()};
    }
  };
} // namespace Sel::detail

namespace Sel {
  /** @class  ParticleCombinationSpan
   *  @tparam iterator_t Type of the iterators to the first and
   *                     one-past-the-last member of the combination.
   *
   *  Type representing a runtime-variable-sized combination of particles. If
   *  the particles can be of heterogeneous type then the value type yielded
   *  by the combination will be a variant type.
   */
  template <typename iterator_t>
  struct ParticleCombinationSpan : public detail::ParticleCombinationBase<ParticleCombinationSpan<iterator_t>> {
    ParticleCombinationSpan( iterator_t begin, iterator_t end ) : m_begin{std::move( begin )}, m_end{end} {}
    friend iterator_t begin( ParticleCombinationSpan const& span ) { return span.m_begin; }
    friend iterator_t end( ParticleCombinationSpan const& span ) { return span.m_end; }

    static auto mask_true() { return detail::mask_true_helper<tuple_t>::_(); }

  private:
    iterator_t m_begin, m_end;
    // either the raw value type or std::variant<Ts...>
    using item_t = typename std::decay_t<detail::remove_reference_wrapper_t<decltype( *std::declval<iterator_t>() )>>;
    // either std::tuple<item_t> or std::tuple<Ts...> if item_t is a variant
    using tuple_t = typename Utils::is_variant<item_t>::tuple_type;
  };

  /** @class  ParticleCombination
   *  @tparam child_ts... Pack of concrete child [proxy] types.
   *
   *  Type representing a compile time constant sized combination of particles
   *  of known type. The given particles are stored by value, which is sensible
   *  for proxy types. In other cases then std::reference_wrapper can be used.
   *
   *  @todo Support subscript and begin/end, returning variant types?
   */
  template <typename... child_ts>
  struct ParticleCombination : public detail::ParticleCombinationBase<ParticleCombination<child_ts...>> {
    static_assert( sizeof...( child_ts ) > 1, "Combinations of fewer than 2 particles don't make sense..." );
    ParticleCombination( child_ts... children ) : m_children{std::move( children )...} {}

    /** Helper type containing all child types with wrapping by
     *  std::reference_wrapper removed.
     */
    using child_types = boost::mp11::mp_list<std::decay_t<detail::remove_reference_wrapper_t<child_ts>>...>;

    /** Get the type of the Nth child.
     */
    template <std::size_t N>
    using child_t = boost::mp11::mp_at_c<child_types, N>;

    /** Helper type containing T const& versions of 'child_types'
     */
    using child_crefs = boost::mp11::mp_transform<std::add_lvalue_reference_t,
                                                  boost::mp11::mp_transform<std::add_const_t, child_types>>;

    /** Helper for determining the result of applying a particular transform to
     *  this combination. We could consider applying std::common_type_t, but to
     *  be consistent we would want to ensure that the same behaviour was coded
     *  into ParticleCombinationSpan.
     */
    template <typename transform_t>
    struct transform_result {
      template <typename child_t>
      using result_t  = LHCb::invoke_or_visit_result_t<transform_t, child_t>;
      using results_t = boost::mp11::mp_transform<result_t, child_crefs>;
      static_assert( boost::mp11::mp_same<results_t>::value,
                     "Sel::transform() was given a functor with different return types for "
                     "different members of a Sel::ParticleCombination." );
      using type = boost::mp11::mp_first<results_t>;
    };

    /** Helper for determining the result of applying a pairwise transform to
     *  this combination. We could consider applying std::common_type_t, but to
     *  be consistent we would want to ensure the same behaviour in
     *  ParticleCombinationSpan, where the transform could be applied as a
     *  visitor on a variant type.
     */
    template <typename transform_t>
    struct pairwise_transform_result {
      template <typename child1_t, typename child2_t>
      using result_t  = LHCb::invoke_or_visit_result_t<transform_t, child1_t, child2_t>;
      using results_t = boost::mp11::mp_product<result_t, child_crefs, child_crefs>;
      static_assert( boost::mp11::mp_same<results_t>::value,
                     "Sel::pairwise_transform() was given a functor with different return types for "
                     "different members of a Sel::ParticleCombination." );
      using type = boost::mp11::mp_first<results_t>;
    };

    /** Get a reference to the Nth child.
     *  If this is std::reference_wrapper<T>, T& is returned.
     */
    /*  template <std::size_t N>
    child_t<N>& get() {
      return std::get<N>( m_children );
      }*/

    /** Get a const reference to the Nth child.
     *  If this is std::reference_wrapper<T [const]>, T const& is returned.
     */
    /*template <std::size_t N>
    child_t<N> const& get() const {
      return std::get<N>( m_children );
    }*/

    /** Unfortunately this is needed in order for the AcceptAll functor (ALL) to work.
     */
    static auto mask_true() { return detail::mask_true_helper<typename std::tuple<child_ts...>>::_(); }

  private:
    // Store children by value -- they're often proxies, and if not then the
    // calling code can use std::reference_wrapper.
    std::tuple<child_ts...> m_children;
  };

  /** Shorthand for a ParticleCombination containing N of the same child type.
   */
  template <typename child_t, std::size_t N>
  using ParticleCombinationN =
      boost::mp11::mp_apply<ParticleCombination, boost::mp11::mp_repeat_c<boost::mp11::mp_list<child_t>, N>>;

  /** Apply the given transform to the combination members given by `Is...` and
   *  return the result.
   */
  template <std::size_t... Is, typename collection_t, typename transform_t>
  auto transform( collection_t const&, transform_t ) {
    static_assert( false && sizeof...( Is ) == 42, "TODO implement this" );
  }

  /** Specialised apply-to-children for ParticleCombination, where more is
   *  known at compile time.
   */
  template <std::size_t... Is, typename... child_ts, typename transform_t>
  auto transform( ParticleCombination<child_ts...> const& comb, transform_t transform ) {
    return LHCb::invoke_or_visit( transform, comb.template get<Is>()... );
  }

  /** Apply a transform-and-reduce operation to a non-empty particle
   *  combination. This returns
   *    reduce( reduce( transform( comb[0] ), transform( comb[1] ) ),
   *            transform( comb[2] ) )
   *  and so on, for all elements of the given combination. In order to
   *  simplify type deduction in client code, it is an error to pass an empty
   *  combination. If the combination only contains one element, the result
   *  of transform( comb[0] ) is converted to the deduced return type of the
   *  reduce operation, even though reduce is not called.
   *
   *  @todo "Divide and conquer"?
   *        e.g. reduce( reduce( v0, v1 ), reduce( v2, v3 ) ) instead of the
   *        current reduce( reduce( reduce( v0, v1 ), v2 ), v3 )
   */
  template <typename combination_t, typename transform_t, typename reduce_t>
  auto transform_reduce( combination_t const& comb, transform_t transform, reduce_t reduce ) {
    using std::begin;
    using std::end;
    auto iter = begin( comb );
    auto endi = end( comb );
    if ( iter == endi ) {
      throw GaudiException{"Empty particle combination -- this should never happen", "Sel::transform",
                           StatusCode::FAILURE};
    }
    // T const& from std::reference_wrapper<T>
    using item_t    = detail::as_const_ref_t<decltype( *iter )>;
    using value_t   = decltype( Utils::invoke_or_visit( transform, std::declval<item_t>() ) );
    using reduced_t = std::invoke_result_t<reduce_t, value_t, value_t>;
    reduced_t value = Utils::invoke_or_visit( transform, item_t{*iter} );
    for ( ++iter; iter != endi; ++iter ) {
      value = std::invoke( reduce, value, Utils::invoke_or_visit( transform, item_t{*iter} ) );
    }
    return value;
  }

  /** Specialised transform-and-reduce for ParticleCombination, where more is
   *  known at compile time.
   *  @todo "Divide and conquer"?
   */
  template <typename... child_ts, typename transform_t, typename reduce_t>
  auto transform_reduce( ParticleCombination<child_ts...> const& comb, transform_t transform, reduce_t reduce ) {
    // When reading this code it's worth thinking about the case that value_t
    // is bool, reduce_t is std::plus<> and reduced_t is int.
    using value_t    = typename ParticleCombination<child_ts...>::template transform_result<transform_t>::type;
    using reduced_t  = std::invoke_result_t<reduce_t, value_t, value_t>;
    using reduced2_t = std::invoke_result_t<reduce_t, reduced_t, value_t>;
    static_assert( std::is_same_v<reduced_t, reduced2_t> );
    // Writing this out explicitly avoids an avoidable value_t -> reduced_t
    // conversion.
    auto value = std::invoke( reduce, value_t{LHCb::invoke_or_visit( transform, comb.template get<0>() )},
                              value_t{LHCb::invoke_or_visit( transform, comb.template get<1>() )} );
    static_assert( std::is_same_v<reduced_t, decltype( value )> );
    LHCb::Utils::unwind<2, sizeof...( child_ts )>( [&]( auto n ) {
      value = std::invoke( reduce, value, value_t{LHCb::invoke_or_visit( transform, comb.template get<n>() )} );
    } );
    return value;
  }

  /** Apply a transform-and-reduce operation to all pairs of particles in the
   *  given combination. It is an error to pass a 0- or 1-body combination.
   *  This effectively returns:
   *   reduce( reduce( transform( comb[0], comb[1] ), transform( comb[0], comb[2] ) ),
   *           transform( comb[1], comb[2] ) )
   *  and so on. In case the combination is 2-body, so there is only one pair
   *  of particles, the result of transform( comb[0], comb[1] ) is converted to
   *  the deduced return type of the reduce operation, even though reduce is
   *  not called.
   *
   *  @todo "Divide and conquer"?
   */
  template <typename collection_t, typename transform_t, typename reduce_t>
  auto pairwise_transform_reduce( collection_t const& comb, transform_t transform, reduce_t reduce ) {
    using std::begin;
    using std::end;
    auto iter = begin( comb );
    auto endi = end( comb );
    if ( iter == endi ) {
      throw GaudiException{"Empty collection of children -- this should never happen", "Sel::pairwise_transform",
                           StatusCode::FAILURE};
    }
    // T const& from std::reference_wrapper<T>
    using item_t     = detail::as_const_ref_t<decltype( *iter )>;
    using value_t    = decltype( Utils::invoke_or_visit( transform, std::declval<item_t>(), std::declval<item_t>() ) );
    using reduced_t  = std::invoke_result_t<reduce_t, value_t, value_t>;
    auto second_iter = iter;
    ++second_iter;
    if ( second_iter == endi ) {
      throw GaudiException{"<2-particle combination -- this should never happen", "Sel::pairwise_transform",
                           StatusCode::FAILURE};
    }
    // Calculate the first value explicitly to avoid having to pass in an
    // initial value of type reduced_t
    reduced_t value = Utils::invoke_or_visit( transform, item_t{*iter}, item_t{*second_iter} ); // (0, 1)
    do {
      ++second_iter;
      for ( ; second_iter != endi; ++second_iter ) {
        value = std::invoke( reduce, value, Utils::invoke_or_visit( transform, item_t{*iter}, item_t{*second_iter} ) );
      }
      ++iter;
      second_iter = iter;
    } while ( iter != endi );
    return value;
  }

  /** Specialised pairwise transform-and-reduce for ParticleCombination,
   *  where more is known at compile time.
   *  @todo "Divide and conquer"?
   */
  template <typename... child_ts, typename transform_t, typename reduce_t>
  auto pairwise_transform_reduce( ParticleCombination<child_ts...> const& comb, transform_t transform,
                                  reduce_t reduce ) {
    using value_t    = typename ParticleCombination<child_ts...>::template pairwise_transform_result<transform_t>::type;
    using reduced_t  = std::invoke_result_t<reduce_t, value_t, value_t>;
    using reduced2_t = std::invoke_result_t<reduce_t, reduced_t, value_t>;
    static_assert( std::is_same_v<reduced_t, reduced2_t> );
    // In principle we could avoid a result_t -> reduced_t conversion here in a
    // >2-body combination.
    reduced_t value = LHCb::invoke_or_visit( transform, comb.template get<0>(), comb.template get<1>() );
    LHCb::Utils::unwind<0, sizeof...( child_ts )>( [&]( auto n ) {
      LHCb::Utils::unwind<n ? n + 1 : 2, sizeof...( child_ts )>( [&]( auto m ) {
        value = std::invoke( reduce, value,
                             LHCb::invoke_or_visit( transform, comb.template get<n>(), comb.template get<m>() ) );
      } );
    } );
    return value;
  }
} // namespace Sel

// Register headers
template <typename iterator_t>
struct LHCb::header_map<Sel::ParticleCombinationSpan<iterator_t>> {
  static constexpr auto value = LHCb::header_map_v<iterator_t> + "SelKernel/ParticleCombination.h";
};

template <typename... child_ts>
struct LHCb::header_map<Sel::ParticleCombination<child_ts...>> {
  static constexpr auto value = ( LHCb::header_map_v<child_ts> + ... ) + "SelKernel/ParticleCombination.h";
};
