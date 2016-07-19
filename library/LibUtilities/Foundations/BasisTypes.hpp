#ifndef LIBUTILITIES_FOUNDATIONS_BASISTYPES
#define LIBUTILITIES_FOUNDATIONS_BASISTYPES

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

struct ModifiedLegendre {};
struct Lagrange {};
struct BernsteinTriangle {};
struct BernsteinTetrahedron {};

namespace traits
{
    // Default values for all types
    struct basis_traits_base {
        static const bool is_tensor_product = false;
        static const unsigned int dimension = 0;
        static const unsigned int get_total_modes(unsigned int nmodes) {
            return nmodes;
        }
        typedef Segment native_shape;
    };


    // Primary template
    template<typename... TPts>
    struct basis_traits : public basis_traits_base {
    };

    template<typename T1, typename... TBasis>
    struct basis_traits<T1, TBasis...> : public basis_traits_base {
        static const unsigned int dimension = basis_traits<T1>::dimension + basis_traits<TBasis...>::dimension;
        static_assert(dimension <= 3, "Composite basis dimension > 3.");
    };

    template<>
    struct basis_traits<ModifiedLegendre> : public basis_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
    };

    template<>
    struct basis_traits<Lagrange> : public basis_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
    };

    template<>
    struct basis_traits<BernsteinTriangle> : public basis_traits_base
    {
        static const bool is_tensor_product = false;
        static const unsigned int dimension = 2;
        typedef Triangle native_shape;
    };

    template<>
    struct basis_traits<BernsteinTetrahedron> : public basis_traits_base
    {
        static const bool is_tensor_product = false;
        static const unsigned int dimension = 3;
        typedef Tetrahedron native_shape;
    };



    struct expansion_traits_base
    {
        static const bool is_valid = false;
    };

    template<typename TShape, typename TBasis>
    struct expansion_traits : public shape_traits<TShape>,
                                 public basis_traits<TBasis>,
                                 public expansion_traits_base
    {
    };

    template<typename TBasis>
    struct expansion_traits<Segment, TBasis> :
                                 public shape_traits<Segment>,
                                 public basis_traits<TBasis>,
                                 public expansion_traits_base
    {
        static const bool is_valid = true;
    };


    template<typename TBasis>
    struct expansion_traits<Quadrilateral, TBasis> :
                                 public shape_traits<Quadrilateral>,
                                 public basis_traits<TBasis>,
                                 public expansion_traits_base
    {
        static const bool is_valid = basis_traits<TBasis>::is_tensor_product;
    };

    template<>
    struct expansion_traits<Triangle, BernsteinTriangle> :
                                public shape_traits<Triangle>,
                                public basis_traits<BernsteinTriangle>,
                                public expansion_traits_base
    {
        static const bool is_valid = true;
    };
}
}
}
}

#endif
