#include <LibUtilities/Foundations/Points.hpp>

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
                                  Loki::NoDestroy,
                                  Loki::ClassLevelLockable> Type;
    return Type::Instance();
}


static std::string points[] = {
        GetPointsFactory().RegisterCreatorFunction("GaussGaussLegendre_Segment", Points<NekDouble, Segment, GaussGaussLegendre>::create),
        GetPointsFactory().RegisterCreatorFunction("GaussRadauMLegendre_Segment", Points<NekDouble, Segment, GaussRadauMLegendre>::create),
};
}
}
}


