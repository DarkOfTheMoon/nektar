#ifndef LIBUTILITIES_FOUNDATIONS_POINTSTRAITS
#define LIBUTILITIES_FOUNDATIONS_POINTSTRAITS

#include <LibUtilities/Foundations/ShapeTypes.hpp>
#include <LibUtilities/Foundations/PointsTypes.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{
namespace traits
{
    // Default values for all types
    struct points_traits_base {
        static const bool is_tensor_product = false;
        static const unsigned int dimension = 0;
        static const unsigned int get_total_points(unsigned int npts) {
            return npts;
        }
        typedef Segment native_shape;
    };


    // Primary template
    template<typename... TPts>
    struct points_traits : public points_traits_base {
    };

    template<typename T1, typename... TPts>
    struct points_traits<T1, TPts...> : public points_traits_base {
        static const unsigned int dimension = points_traits<T1>::dimension + points_traits<TPts...>::dimension;
        static_assert(dimension <= 3, "Composite points dimension > 3.");
    };

    template<>
    struct points_traits<GaussGaussLegendre> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
        static const std::string name;
        static const int id;
        typedef Segment native_shape;
    };

    template<>
    struct points_traits<GaussRadauMLegendre> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
        typedef Segment native_shape;
    };

    template<>
    struct points_traits<GaussRadauPLegendre> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
        typedef Segment native_shape;
    };

    template<>
    struct points_traits<GaussLobattoLegendre> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
        typedef Segment native_shape;
    };

    template<>
    struct points_traits<Fourier> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
        typedef Segment native_shape;
    };

    template<>
    struct points_traits<Fekete> : public points_traits_base
    {
        static const bool is_tensor_product = false;
        static const unsigned int dimension = 2;
        static const unsigned int get_total_points(unsigned int npts) {
            return npts*(npts+1)/2;
        }
        typedef Triangle native_shape;
    };

    struct distribution_traits_base
    {
        static const bool is_valid = false;
    };

    template<typename TShape, typename TPts>
    struct distribution_traits : public shape_traits<TShape>,
                                 public points_traits<TPts>,
                                 public distribution_traits_base
    {
    };

    template<typename TPts>
    struct distribution_traits<Segment, TPts> :
                                 public shape_traits<Segment>,
                                 public points_traits<TPts>,
                                 public distribution_traits_base
    {
        static const bool is_valid = true;
    };

    template<typename TPts>
    struct distribution_traits<Quadrilateral, TPts> :
                                 public shape_traits<Quadrilateral>,
                                 public points_traits<TPts>,
                                 public distribution_traits_base
    {
        static const bool is_valid = points_traits<TPts>::is_tensor_product;
    };

    template<>
    struct distribution_traits<Triangle, Fekete> :
                                public shape_traits<Triangle>,
                                public points_traits<Fekete>,
                                public distribution_traits_base
    {
        static const bool is_valid = true;
    };
}
}
}
}

#endif
