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

    BasisKey bkey;
    bkey.m_nummodes[0] = 4;
    bkey.m_nummodes[1] = 5;
    bkey.m_nummodes[2] = 6;

    // Test some points keys
    Points<double, Segment,       GaussGaussLegendre> P(key);
    Points<double, Quadrilateral, GaussGaussLegendre, GaussGaussLegendre> Pquad(key);
    Points<double, Hexahedron,    GaussGaussLegendre, GaussGaussLegendre, GaussGaussLegendre> Phex(key);

    // These should fail
    //Points<double, Segment,       Fekete> Pseg_fek(key);
    //Points<double, Triangle,      GaussGaussLegendre> Ptri_gauss(key);

    // Test extraction of constituent points distributions
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

    // Test access through base class pointer
    PointsBase<double>* ptr1 = new Points<double, Segment,       GaussGaussLegendre>(key);
    PointsBase<double>* ptr2 = new Points<double, Quadrilateral, GaussGaussLegendre, GaussGaussLegendre>(key);
    cout << ptr1->GetZ(0).num_elements() << endl;
    cout << ptr2->GetZ(0).num_elements() << endl;
    cout << ptr2->GetZ(1).num_elements() << endl;

//    // Test accessing undefined points distribution
//    try {
//        cout << ptr2->GetZ(2).num_elements() << endl;
//    }
//    catch (...) {
//        cout << "Unable to access 3rd GGL distribution in a quad." << endl;
//    }

    // Test instantiation of a Basis
    Basis<double, Segment, Points<double, Segment, GaussGaussLegendre>, ModifiedLegendre> B0(bkey);

    // These should fail
    //Basis<double, Triangle, Points<double, Segment, GaussGaussLegendre>, BernsteinTriangle> Bbtri(bkey);
    //Basis<double, Triangle, Points<double, Segment, GaussGaussLegendre>, BernsteinTriangle> Bbtri(bkey);
    //Basis<double, Triangle, Points<double, Triangle, GaussGaussLegendre>, BernsteinTriangle> Bbtri(bkey);

//    cout << sizeof(Array<OneD, double>) << endl;
//    cout << sizeof(PointsKey) << endl;
    cout << "Soze of GGL points: " << sizeof(Points<double, Segment, GaussGaussLegendre>) << endl;
//    cout << sizeof(Points<double, Quadrilateral, GaussGaussLegendre, GaussGaussLegendre>) << endl;
//    cout << traits::points_traits<GaussGaussLegendre,GaussGaussLegendre,GaussGaussLegendre>::dimension << endl;
//


    // Create a points object using the factory
    PointsSharedPtr<double> facptr = GetPointsFactory<double>().CreateInstance("GaussGaussLegendre_Segment", key);
    cout << facptr->GetNumPoints() << endl;

    cout << "Size of Modified basis: " << sizeof(Basis<double, Segment, Points<double, Segment, GaussGaussLegendre>, ModifiedLegendre>) << endl;
    // This one should cause a compile error
    //Points<Quadrilateral, Fekete, double> Pfail(key);

    //Points<Triangle, GaussGaussLegendre, double> P2(key);
}
