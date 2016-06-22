#ifndef LIBUTILITIES_FOUNDATIONS_SHAPETYPES
#define LIBUTILITIES_FOUNDATIONS_SHAPETYPES

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

}
}
}
}

#endif
