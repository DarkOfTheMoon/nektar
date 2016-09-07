#include <LibUtilities/Foundations/Basis/BasisKey.h>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

LIB_UTILITIES_EXPORT std::ostream& operator<<(std::ostream& os, const BasisKey& rhs)
{
    os << "Print key" << std::endl;
    return os;
}


}
}
}
