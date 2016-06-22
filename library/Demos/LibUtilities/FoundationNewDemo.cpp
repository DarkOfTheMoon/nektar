#include <iostream>
using namespace std;

#include <LibUtilities/FoundationsNew/Points.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities::Foundations;

int main () {
    PointsKey key;
    key.m_numpoints[0] = 5;
    key.m_numpoints[1] = 6;
    key.m_numpoints[2] = 7;

    //Points<double, Quadrilateral> Pnotype(key);
    Points<double, Quadrilateral, GaussGaussLegendre> P(key);
    Points<double, Quadrilateral, GaussGaussLegendre, GaussGaussLegendre> Pquad(key);

    Array<OneD, double> p1 = Pquad.GetZ(0);
    Array<OneD, double> p2 = Pquad.GetZ(1);
    cout << p1.num_elements() << endl;
    cout << p2.num_elements() << endl;

    cout << sizeof(Array<OneD, double>) << endl;
    cout << sizeof(PointsKey) << endl;
    cout << sizeof(Points<double, Quadrilateral, GaussGaussLegendre>) << endl;
    cout << sizeof(Points<double, Quadrilateral, GaussGaussLegendre, GaussGaussLegendre>) << endl;
    cout << traits::points_traits<GaussGaussLegendre,GaussGaussLegendre,GaussGaussLegendre>::dimension << endl;
    // This one should cause a compile error
    //Points<Quadrilateral, Fekete, double> Pfail(key);

    //Points<Triangle, GaussGaussLegendre, double> P2(key);
}
