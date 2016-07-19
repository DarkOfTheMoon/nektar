#include <LibUtilities/Foundations/Basis.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{
static std::string bases[] = {
        GetBasisFactory<NekDouble>().RegisterCreatorFunction("ModifiedLegendre_GaussGaussLegendre_Segment", Basis<NekDouble, Segment, std::tuple<GaussGaussLegendre>, std::tuple<ModifiedLegendre>>::create),
        GetBasisFactory<NekDouble>().RegisterCreatorFunction("ModifiedLegendre_GaussRadauMLegendre_Segment", Basis<NekDouble, Segment, std::tuple<GaussRadauMLegendre>, std::tuple<ModifiedLegendre>>::create),
        GetBasisFactory<NekDouble>().RegisterCreatorFunction("ModifiedLegendre,ModifiedLegendre_GaussGaussLegendre,GaussGaussLegendre_Quadrilateral", Basis<NekDouble, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>, std::tuple<ModifiedLegendre, ModifiedLegendre>>::create)
};
}
}
}
