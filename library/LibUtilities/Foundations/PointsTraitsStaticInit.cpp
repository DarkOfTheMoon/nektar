#include <LibUtilities/Foundations/PointsTraits.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

PointsStringToIdMap& GetPointsStringToIdMap()
{
    typedef Loki::SingletonHolder<PointsStringToIdMap,
        Loki::CreateUsingNew,
        Loki::NoDestroy,
        Loki::ClassLevelLockable> Type;
    return Type::Instance();
}

PointsIdToStringMap& GetPointsIdToStringMap()
{
    typedef Loki::SingletonHolder<PointsIdToStringMap,
        Loki::CreateUsingNew,
        Loki::NoDestroy,
        Loki::ClassLevelLockable> Type;
    return Type::Instance();
}

namespace traits
{
const std::string points_traits<GaussGaussLegendre>::name = GetPointsStringToIdMap().insert(std::pair<std::string, int>("GaussGaussLegendre", 1)).first->first;
const int points_traits<GaussGaussLegendre>::id = GetPointsIdToStringMap().insert(std::pair<int, std::string>(1, "GaussGaussLegendre")).first->first;
}

}
}
}
