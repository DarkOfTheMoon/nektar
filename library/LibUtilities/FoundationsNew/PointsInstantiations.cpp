#include <LibUtilities/FoundationsNew/Points.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{
static std::string points[] = {
        GetPointsFactory<double>().RegisterCreatorFunction("GaussGaussLegendre_Segment", Points<double, Segment, GaussGaussLegendre>::create),
        GetPointsFactory<double>().RegisterCreatorFunction("GaussRadauMLegendre_Segment", Points<double, Segment, GaussRadauMLegendre>::create),
};
}
}
}


