#include <LibUtilities/Foundations/Basis.hpp>

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
        GetBasisFactory<NekDouble>().RegisterCreatorFunction("MOD_GGL_Segment", Basis<NekDouble, Segment, std::tuple<GaussGaussLegendre>, std::tuple<ModifiedLegendre>>::create),
        GetBasisFactory<NekDouble>().RegisterCreatorFunction("MOD_GRL_Segment", Basis<NekDouble, Segment, std::tuple<GaussRadauMLegendre>, std::tuple<ModifiedLegendre>>::create),
        GetBasisFactory<NekDouble>().RegisterCreatorFunction("MOD,MOD_GGL,GGL_Quadrilateral", Basis<NekDouble, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>, std::tuple<ModifiedLegendre, ModifiedLegendre>>::create)
};

}
}
}
