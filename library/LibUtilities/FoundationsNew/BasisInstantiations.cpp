#include <LibUtilities/FoundationsNew/Basis.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{
static std::string bases[] = {
        GetBasisFactory<NekDouble>().RegisterCreatorFunction("ModifiedLegendre_GaussGaussLegendre_Segment", Basis<NekDouble, Segment, Points<NekDouble, Segment, GaussGaussLegendre>, ModifiedLegendre>::create),
        GetBasisFactory<NekDouble>().RegisterCreatorFunction("ModifiedLegendre_GaussRadauMLegendre_Segment", Basis<NekDouble, Segment, Points<NekDouble, Segment, GaussRadauMLegendre>, ModifiedLegendre>::create),
};
}
}
}
