#include <LibUtilities/FoundationsNew/Points.hpp>

using namespace Nektar::LibUtilities::Foundations;

int main () {
    PointsKey key;
    key.m_numpoints[0] = 5;
    key.m_numpoints[1] = 5;
    key.m_numpoints[2] = 1;

    Points<Quadrilateral, GaussGaussLegendre, double> P(key);

    // This one should cause a compile error
    //Points<Quadrilateral, Fekete, double> Pfail(key);

    //Points<Triangle, GaussGaussLegendre, double> P2(key);
}
