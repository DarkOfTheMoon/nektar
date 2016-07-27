#ifndef LIBRARY_LIBUTILITIES_FOUNDATIONS_FOUNDATIONS_HPP_
#define LIBRARY_LIBUTILITIES_FOUNDATIONS_FOUNDATIONS_HPP_

#include <map>

// This is just a collection of useful templates, typedefs and helper classes

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

template <class T>
struct is_not_tuple { static const bool value = true; };

template <typename... Args>
struct is_not_tuple<std::tuple<Args...>> { static const bool value = false; };

template <typename... Args>
struct is_not_tuple<const std::tuple<Args...>> {
        static const bool value = false;
};

typedef double NekDouble;
typedef std::string BasisParamKey;
typedef NekDouble BasisParamValue;
typedef std::map<BasisParamKey, BasisParamValue> BasisParamList;

// 3-character identifier class
template<typename TData>
class TriArray
{
    public:
        TriArray()
        {
            m_data[0] = ' ';
            m_data[1] = ' ';
            m_data[2] = ' ';
        }

        template<int N, typename = typename std::enable_if<N==4, void>::type >
        TriArray(TData (&pSrc)[N])
        {
            m_data[0] = *pSrc;
            m_data[1] = *(pSrc+1);
            m_data[2] = *(pSrc+2);
        }
        TriArray(const TriArray& pSrc)
        {
            memcpy(m_data, pSrc.m_data, 3);
        }
        ~TriArray() {}

        TriArray<TData> &operator=(const TriArray& pSrc)
        {
            m_data[0] = pSrc.m_data[0];
            m_data[1] = pSrc.m_data[1];
            m_data[2] = pSrc.m_data[2];
            return *this;
        }

        template<int N, typename = typename std::enable_if<N==4, void>::type >
        TriArray<TData> &operator=(TData (&pSrc)[N])
        {
            m_data[0] = *pSrc;
            m_data[1] = *(pSrc+1);
            m_data[2] = *(pSrc+2);
            return *this;
        }

        bool operator==(const TriArray& pSrc) const
        {
            return (m_data[0] == pSrc.m_data[0] &&
                    m_data[1] == pSrc.m_data[1] &&
                    m_data[2] == pSrc.m_data[2]);
        }

        bool operator!=(const TriArray& pSrc) const
        {
            return !operator==(pSrc);
        }

        bool operator<(const TriArray& pSrc) const
        {
            return (m_data[0] <= pSrc.m_data[0] &&
                    m_data[1] <= pSrc.m_data[1] &&
                    m_data[2] <  pSrc.m_data[2]);
        }

        const std::string GetString() const
        {
            return std::string(m_data);
        }
    private:
        TData m_data[3];
};

}
}
}

#endif /* LIBRARY_LIBUTILITIES_FOUNDATIONS_FOUNDATIONS_HPP_ */
