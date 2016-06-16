#include <iostream>

#include "PointsTypes.hpp"

namespace Nektar
{
namespace LibUtilities
{

// Dummy Array class for the purpose of isolated development
struct OneD {};
template<typename dim, typename datatype>
class Array {};

namespace Foundations
{


//
//template<bool b, typename T = void>
//struct assert_if {
//    typedef void type;
//    static_assert(b, "That is not a valid combination, silly!");
//};
//
//template<typename T>
//struct assert_if<true, T> {
//    typedef T type;
//};


struct PointsKey
{
      std::string name;
      int npoints[3];
};


/**
 * @class PointsBase
 * @brief Points base class defining interface and data members
 */
template<typename datatype>
class PointsBase
{
    public:
        Array<OneD, datatype> GetPoints()
        {
            return m_points;
        }

    protected:
        Array<OneD, datatype> m_points;
        Array<OneD, datatype> m_weights;

        PointsBase() : m_points(), m_weights() {}
};

// Primary template for Points classes.
template<typename shape, typename pointstype, typename datatype,
        class Enable = void>
class Points
{
    public:
        Points() {}
};


// Triangular Points distribution specialisation
// This adds
template<typename pointstype, typename datatype>
class Points<Triangle, pointstype, datatype, void> : public PointsBase<datatype>
{
    public:
        Points(int np1, int np2) {
            // do something?
            PointsBase<datatype>::m_points = Array<OneD, datatype>();
        }
};

// Quadrilateral Points distribution specialisation
template<typename pointstype, typename datatype>
class Points<Quadrilateral, pointstype, datatype> : public PointsBase<datatype>
{
        static_assert(traits::points_traits<Quadrilateral, pointstype>::is_valid,
                "This combination is not valid.");

    public:
        Points(int np1, int np2) {
            // do something?
            PointsBase<datatype>::m_points = Array<OneD, datatype>();
        }
};


}
}
}
