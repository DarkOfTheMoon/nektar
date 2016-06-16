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

struct Tensor {};
struct Fekete {};


namespace traits
{
    // Default values for all types
    struct points_traits_base {
        static const bool is_tensor_product = false;
        static const bool is_valid = true;
        static const int  dimension = 0;
    };

    // Primary template
    template<typename shape, typename pointstype>
    struct points_traits : public points_traits_base {
    };


    // Tensor product types
    template<typename shape>
    struct points_traits<shape, Tensor> : public points_traits_base
    {
        static const bool is_tensor_product = true;
    };

    template<>
    struct points_traits<Quadrilateral, Fekete> : public points_traits_base
    {
        static const bool is_valid = false;
    };
}
}
}
}
