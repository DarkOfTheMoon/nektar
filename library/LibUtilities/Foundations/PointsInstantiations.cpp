#include <LibUtilities/Foundations/Points.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{
static std::string points[] = {
        GetPointsFactory<NekDouble>().RegisterCreatorFunction("GaussGaussLegendre_Segment", Points<NekDouble, Segment, GaussGaussLegendre>::create),
        GetPointsFactory<NekDouble>().RegisterCreatorFunction("GaussRadauMLegendre_Segment", Points<NekDouble, Segment, GaussRadauMLegendre>::create),
};
}
}
}


