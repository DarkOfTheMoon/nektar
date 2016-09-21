#ifndef LIBRARY_LIBUTILITIES_FOUNDATIONS_BASISKEY_H_
#define LIBRARY_LIBUTILITIES_FOUNDATIONS_BASISKEY_H_

#include <iostream>

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{


class BasisKey
{
    public:
        friend bool operator<(const BasisKey &lhs, const BasisKey &rhs)
        {
            return (lhs.m_id < rhs.m_id);
        }

    public:
        std::string m_id;
        unsigned int m_dim;
        //unsigned char m_basistype[3][3];
        unsigned int m_nummodes[3];
        BasisParamList m_params;
        PointsKey m_ptsKey;

        struct opLess
        {
            LIB_UTILITIES_EXPORT  bool operator()(const BasisKey &lhs, const BasisKey &rhs) const
            {
                return operator<(lhs, rhs);
            }
        };

};

LIB_UTILITIES_EXPORT std::ostream& operator<<(std::ostream& os, const BasisKey& rhs);

}
}
}

#endif /* LIBRARY_LIBUTILITIES_FOUNDATIONS_BASISKEY_H_ */
