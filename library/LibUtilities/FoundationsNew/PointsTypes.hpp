namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{
struct Segment {};
struct Triangle {};
struct Quadrilateral {};
struct Hexahedron {};
struct Tetrahedron {};
struct Prism {};
struct Pyramid {};

struct GaussGaussLegendre {};
struct GaussRadauMLegendre {};
struct GaussRadauPLegendre {};
struct GaussLobattoLegendre {};
struct Fekete {};
struct Fourier {};


namespace traits
{
    struct shape_traits_base {
            static const int dimension = 0;
            static const bool is_simplex = false;
            static constexpr char name[] = "";
    };

    // Primary shape traits template
    template<typename TShape>
    struct shape_traits : public shape_traits_base {
    };

    template<>
    struct shape_traits<Segment> : public shape_traits_base {
            static const int dimension = 1;
            static const bool is_simplex = false;
            static constexpr char name[] = "Segment";
    };

    template<>
    struct shape_traits<Triangle> : public shape_traits_base {
            static const int dimension = 2;
            static const bool is_simplex = true;
            static constexpr char name[] = "Triangle";
    };

    template<>
    struct shape_traits<Quadrilateral> : public shape_traits_base {
            static const int dimension = 2;
            static const bool is_simplex = false;
            static constexpr char name[] = "Quadrilateral";
    };

    template<>
    struct shape_traits<Tetrahedron> : public shape_traits_base {
            static const int dimension = 3;
            static const bool is_simplex = true;
            static constexpr char name[] = "Tetrahedron";
    };

    template<>
    struct shape_traits<Prism> : public shape_traits_base {
            static const int dimension = 3;
            static const bool is_simplex = true;
            static constexpr char name[] = "Prism";
    };

    template<>
    struct shape_traits<Pyramid> : public shape_traits_base {
            static const int dimension = 3;
            static const bool is_simplex = true;
            static constexpr char name[] = "Pyramid";
    };

    template<>
    struct shape_traits<Hexahedron> : public shape_traits_base {
            static const int dimension = 3;
            static const bool is_simplex = false;
            static constexpr char name[] = "Hexahedron";
    };


    // Default values for all types
    struct points_traits_base {
        static const bool is_tensor_product = false;
        static const unsigned int dimension = 0;
        static const unsigned int get_total_points(unsigned int npts) {
            return npts;
        }
    };

    // Primary template
    template<typename TPts>
    struct points_traits : public points_traits_base {
    };

    template<>
    struct points_traits<GaussGaussLegendre> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
    };

    template<>
    struct points_traits<GaussRadauMLegendre> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
    };

    template<>
    struct points_traits<GaussRadauPLegendre> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
    };

    template<>
    struct points_traits<GaussLobattoLegendre> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
    };

    template<>
    struct points_traits<Fourier> : public points_traits_base
    {
        static const bool is_tensor_product = true;
        static const unsigned int dimension = 1;
    };

    template<>
    struct points_traits<Fekete> : public points_traits_base
    {
        static const bool is_tensor_product = false;
        static const unsigned int dimension = 2;
        static const unsigned int get_total_points(unsigned int npts) {
            return npts*(npts+1)/2;
        }
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
    struct distribution_traits<Quadrilateral, TPts> :
                                 public shape_traits<Quadrilateral>,
                                 public points_traits<TPts>,
                                 public distribution_traits_base
    {
        static const bool is_valid = points_traits<TPts>::is_tensor_product;
    };

    template<>
    struct distribution_traits<Triangle, Fekete> :
                                public shape_traits<Quadrilateral>,
                                public points_traits<Fekete>,
                                public distribution_traits_base
    {
        static const bool is_valid = true;
    };
}
}
}
}
