#include <LibUtilities/Foundations/ShapeTypes.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{
namespace traits
{
    // These are needed as char[] is not integer or enum type
    constexpr char shape_traits<Segment>::name[];
    constexpr char shape_traits<Triangle>::name[];
    constexpr char shape_traits<Quadrilateral>::name[];
    constexpr char shape_traits<Tetrahedron>::name[];
    constexpr char shape_traits<Prism>::name[];
    constexpr char shape_traits<Pyramid>::name[];
    constexpr char shape_traits<Hexahedron>::name[];
}
}
}
}
