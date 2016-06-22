#include <iostream>
using namespace std;

#include <LibUtilities/FoundationsNew/Points.hpp>
#include <LibUtilities/FoundationsNew/Basis.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities::Foundations;

int main () {
    PointsKey key;
    key.m_numpoints[0] = 5;
    key.m_numpoints[1] = 6;
    key.m_numpoints[2] = 7;

    //Points<double, Quadrilateral> Pnotype(key);
    Points<double, Segment, GaussGaussLegendre> P(key);
    Points<double, Quadrilateral, GaussGaussLegendre, GaussGaussLegendre> Pquad(key);
    Points<double, Hexahedron, GaussGaussLegendre, GaussGaussLegendre, GaussGaussLegendre> Phex(key);

    Array<OneD, double> p1 = Pquad.GetZ(0);
    Array<OneD, double> p2 = Pquad.GetZ(1);
    cout << p1.num_elements() << endl;
    cout << p2.num_elements() << endl;

    Array<OneD, double> h1 = Phex.GetZ(0);
    Array<OneD, double> h2 = Phex.GetZ(1);
    Array<OneD, double> h3 = Phex.GetZ(2);
    cout << h1.num_elements() << endl;
    cout << h2.num_elements() << endl;
    cout << h3.num_elements() << endl;

    PointsBase<double>* ptr1 = new Points<double, Segment, GaussGaussLegendre>(key);
    PointsBase<double>* ptr2 = new Points<double, Quadrilateral, GaussGaussLegendre, GaussGaussLegendre>(key);

    cout << ptr1->GetZ(0).num_elements() << endl;
    cout << ptr2->GetZ(0).num_elements() << endl;
    cout << ptr2->GetZ(1).num_elements() << endl;
//    cout << sizeof(Array<OneD, double>) << endl;
//    cout << sizeof(PointsKey) << endl;
//    cout << sizeof(Points<double, Quadrilateral, GaussGaussLegendre>) << endl;
//    cout << sizeof(Points<double, Quadrilateral, GaussGaussLegendre, GaussGaussLegendre>) << endl;
//    cout << traits::points_traits<GaussGaussLegendre,GaussGaussLegendre,GaussGaussLegendre>::dimension << endl;
//


    // This one should cause a compile error
    //Points<Quadrilateral, Fekete, double> Pfail(key);

    //Points<Triangle, GaussGaussLegendre, double> P2(key);
}
