#ifndef LIBUTILITIES_FOUNDATIONS_SHAPETYPES
#define LIBUTILITIES_FOUNDATIONS_SHAPETYPES

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

// Types of geometry types.
enum ShapeType
{
    eNoShapeType,
    ePoint,
    eSegment,
    eTriangle,
    eQuadrilateral,
    eTetrahedron,
    ePyramid,
    ePrism,
    eHexahedron,
    SIZE_ShapeType
};

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
            static const enum ShapeType type = eNoShapeType;
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
            static const enum ShapeType type = eSegment;
    };

    template<>
    struct shape_traits<Triangle> : public shape_traits_base {
            static const int dimension = 2;
            static const bool is_simplex = true;
            static constexpr char name[] = "Triangle";
            static const enum ShapeType type = eTriangle;
    };

    template<>
    struct shape_traits<Quadrilateral> : public shape_traits_base {
            static const int dimension = 2;
            static const bool is_simplex = false;
            static constexpr char name[] = "Quadrilateral";
            static const enum ShapeType type = eQuadrilateral;
    };

    template<>
    struct shape_traits<Tetrahedron> : public shape_traits_base {
            static const int dimension = 3;
            static const bool is_simplex = true;
            static constexpr char name[] = "Tetrahedron";
            static const enum ShapeType type = eTetrahedron;
    };

    template<>
    struct shape_traits<Prism> : public shape_traits_base {
            static const int dimension = 3;
            static const bool is_simplex = true;
            static constexpr char name[] = "Prism";
            static const enum ShapeType type = ePrism;
    };

    template<>
    struct shape_traits<Pyramid> : public shape_traits_base {
            static const int dimension = 3;
            static const bool is_simplex = true;
            static constexpr char name[] = "Pyramid";
            static const enum ShapeType type = ePyramid;
    };

    template<>
    struct shape_traits<Hexahedron> : public shape_traits_base {
            static const int dimension = 3;
            static const bool is_simplex = false;
            static constexpr char name[] = "Hexahedron";
            static const enum ShapeType type = eHexahedron;
    };

}


}
}
}

#endif
