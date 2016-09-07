#ifndef LIBUTILITIES_FOUNDATIONS_POINTSTYPES
#define LIBUTILITIES_FOUNDATIONS_POINTSTYPES

#include <unordered_map>
#include <loki/Singleton.h>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

struct GaussGaussLegendre {};
struct GaussRadauMLegendre {};
struct GaussRadauPLegendre {};
struct GaussLobattoLegendre {};
struct Fekete {};
struct Fourier {};

typedef std::unordered_map<std::string, int> PointsStringToIdMap;
typedef std::unordered_map<int, std::string> PointsIdToStringMap;

PointsStringToIdMap& GetPointsStringToIdMap();
PointsIdToStringMap& GetPointsIdToStringMap();

}
}
}

#endif
