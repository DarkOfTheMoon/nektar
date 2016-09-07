#include <LibUtilities/Foundations/Points/PointsGauss.hpp>
#include <LibUtilities/Foundations/Points/PointsFekete.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{


PointsFactory<NekDouble>& GetPointsFactory()
{
    typedef Loki::SingletonHolder<PointsFactory<NekDouble>,
                                  Loki::CreateUsingNew,
                                  Loki::NoDestroy> Type;

    // Loki::ClassLevelLockable

    return Type::Instance();
}


static std::string points[] = {
        GetPointsFactory().RegisterCreatorFunction("GGL,Segment", Points<NekDouble, Segment, GaussGaussLegendre>::create,
                                                   "Gauss-Gauss-Legendre on Segment"),
        GetPointsFactory().RegisterCreatorFunction("GRM,Segment", Points<NekDouble, Segment, GaussRadauMLegendre>::create,
                                                   "Gauss-RadauM-Legendre on Segment"),
};
}
}
}


