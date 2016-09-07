#include <LibUtilities/Foundations/Basis.hpp>
#include <LibUtilities/Foundations/Basis/BasisGauss.hpp>
#include <LibUtilities/Foundations/Basis/BasisBernstein.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

BasisFactory<NekDouble>& GetBasisFactory()
{
    typedef Loki::SingletonHolder<BasisFactory<NekDouble>,
                                  Loki::CreateUsingNew,
                                  Loki::NoDestroy> Type;
    // Putting Loki::ClassLevelLockable causes an assertion!
    return Type::Instance();
}

static std::string bases[] = {
        GetBasisFactory().RegisterCreatorFunction("MOD,GGL,Segment", Basis<NekDouble, Segment, std::tuple<GaussGaussLegendre>, std::tuple<ModifiedLegendre>>::create,
                                                  "Modified basis with Gauss-Gauss-Legendre on Segment"),
        GetBasisFactory().RegisterCreatorFunction("MOD,GRM,Segment", Basis<NekDouble, Segment, std::tuple<GaussRadauMLegendre>, std::tuple<ModifiedLegendre>>::create,
                                                  "Modified basis with Gauss-RadauM-Legendre on Segment"),
        GetBasisFactory().RegisterCreatorFunction("MOD_MOD,GGL_GGL,Quadrilateral", Basis<NekDouble, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>, std::tuple<ModifiedLegendre, ModifiedLegendre>>::create,
                                                  "Modified basis with Gauss-Gauss-Legendre on Quadrilateral")
};

}
}
}
