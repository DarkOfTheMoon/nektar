///////////////////////////////////////////////////////////////////////////////
//
// File: FoundationNewDemo.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Demonstrator for new basis classes
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <type_traits>
using namespace std;

#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.hpp>
#include <LibUtilities/Foundations/Basis.hpp>
#include <LibUtilities/Foundations/Basis/BasisTemplate.hpp>
#include <LibUtilities/Foundations/Basis/BasisGauss.hpp>
#include <LibUtilities/Foundations/Basis/BasisBernstein.hpp>
#include <LibUtilities/Foundations/Interp.hpp>
#include <LibUtilities/Foundations/Points/PointsTraits.hpp>
using namespace Nektar;
using namespace Nektar::LibUtilities::Foundations;


/**
 * This program demonstrates and tests the use of the new points and basis
 * classes.
 */
int main () {
    // Traits are used to hold information about shapes, bases, point
    // distributions, etc.
    cout << "Triangle dim: " << traits::shape_traits<Triangle>::dimension << endl;

    // Create a key for instantiating a Points object. This simply holds the
    // number of points in the points distribution in up to three dimensions.
    // Unused dimensions are ignored, so instantiating a 2D points distribution
    // will ignore m_numpoints[2], for example.
    PointsKey key;
    key.m_numpoints[0] = 5;
    key.m_numpoints[1] = 6;
    key.m_numpoints[2] = 7;

    // We create a key for instantiating a Basis object. As well as the number
    // of modes in each direction, we store an associated PointsKey and an ID
    // descriptor which specified the type of Bases, Points and Shape that the
    // basis key describes. For composite bases, the types of each constituent
    // basis are separated by the underscore (_). For example:
    //    MOD_MOD,GLL_GLL,Quadrilateral
    BasisKey bkey;
    bkey.m_id = "MOD,GGL,Segment";
    bkey.m_ptsKey = key;
    bkey.m_nummodes[0] = 4;
    bkey.m_nummodes[1] = 5;
    bkey.m_nummodes[2] = 6;

    // Instantiate a Points object for a segment using Gauss-Gauss-Legendre
    // point distributions
    Points<double, Segment,       GaussGaussLegendre> P(key);

    // Instantiate a composite Points object for a quadrilateral using
    // Gauss-Gauss-Legendre point distributions in each direction
    Points<double, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>> Pquad(key);

    // Instantiate a composite Points object for a hexahedron.
    Points<double, Hexahedron,    std::tuple<GaussGaussLegendre, GaussGaussLegendre, GaussGaussLegendre>> Phex(key);

    // Points traits describe at compile-time, properties of a Points class.
    // These can be tested to ensure consistency of the template parameters.
    // For example, to ensure the dimension of the Shape parameter, matches that
    // of the provided Point distributions. The following examples, should
    // therefore fail to compile and produce a static assertion error.

    // -- A segment with a triangle Fekete electrostatic distribution
    // Points<double, Segment,       Fekete> Pseg_fek(key);

    // -- A triangle with a 1D distribution
    // Points<double, Triangle,      GaussGaussLegendre> Ptri_gauss(key);

    // -- A segment with a 2D composite distribution
    // Points<double, Segment, std::tuple<GaussGaussLegendre, GaussGaussLegendre>> Pseg_gauss_gauss(key);

    // In composite point distributions, the points in the different directions
    // can be retrieved with an index to GetZ. This function will throw a
    // std::logic_error if the index is out of range.
    Array<OneD, double> p1 = Pquad.GetZ(0);
    Array<OneD, double> p2 = Pquad.GetZ(1);
    cout << "Size in quad dim0: " << p1.num_elements() << endl;
    cout << "Size in quad dim1: " << p2.num_elements() << endl;

    Array<OneD, double> h1 = Phex.GetZ(0);
    Array<OneD, double> h2 = Phex.GetZ(1);
    Array<OneD, double> h3 = Phex.GetZ(2);
    cout << "Size in hex dim0: " << h1.num_elements() << endl;
    cout << "Size in hex dim1: " << h2.num_elements() << endl;
    cout << "Size in hex dim2: " << h3.num_elements() << endl;

    // To use the classes with the Factory pattern, these methods must be
    // accessible through the base class pointer. Here we test that methods can
    // be called successfully through the base class pointer.
    PointsBase<double>* ptr1 = new Points<double, Segment,       GaussGaussLegendre>(key);
    PointsBase<double>* ptr2 = new Points<double, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>>(key);
    cout << "Size in seg ptr dim0: " << ptr1->GetZ(0).num_elements() << endl;
    cout << "Size in quad ptr dim0: " << ptr2->GetZ(0).num_elements() << endl;
    cout << "Size in quad ptr dim1: " << ptr2->GetZ(1).num_elements() << endl;
    // We test that accessing an index out of range indeed produces the error.
    try {
        cout << ptr2->GetZ(2).num_elements() << endl;
    }
    catch (const std::logic_error& e) {
        cout << "Error message produced: " << e.what() << endl;
    }

    // A Basis is defined in a similar manner to Points. Template parameters
    // specify the data type, target shape, points distribution and basis type.
    // A one-dimensional and two-dimensional basis may be instantiated as below.
    // Again, the dimension of the (possibly composite) basis type, must match
    // that of the shape and the (possibly composite) points distribution. In
    // addition, the dimension of each constituent basis must match the
    // dimension of the corresponding points distribution.
    Basis<double, Segment, GaussGaussLegendre, ModifiedLegendre> B0(bkey);
    Basis<double, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>, std::tuple<ModifiedLegendre, ModifiedLegendre>> Bquad(bkey);

    // -- Here we try to instantiate a Bernstein triangle basis on a triangle
    //    but with a 1D point distribution, which should fail.
    // Basis<double, Triangle, GaussGaussLegendre, BernsteinTriangle> Bbtri(bkey);

    // We first print out the available classes in the Points factory
    GetPointsFactory().PrintAvailableClasses(std::cout);

    // We can then instantiate a Points object using the factory. The key is a
    // string which specified the point distribution and shape.
    PointsSharedPtr<double> facptr = GetPointsFactory().CreateInstance("GGL,Segment", key);

    // We can similarly instantiate a Basis object using the factory. The key
    // is again a string (same as what is stored in the basis key).
    GetBasisFactory().PrintAvailableClasses(std::cout);
    BasisSharedPtr<double> bfacptr = GetBasisFactory().CreateInstance("MOD,GGL,Segment", bkey);
    BasisSharedPtr<double> bfacquadptr = GetBasisFactory().CreateInstance("MOD_MOD,GGL_GGL,Quadrilateral", bkey);

    // We can access shape traits through the pointer
    cout << "Dim of quad is: " << bfacquadptr->GetShapeDimension() << endl;

    // A manager can also be used to instantiate Bases. In this case, the
    // manager uses the string stored in the basis key to actually instantiate
    // the object.
    BasisSharedPtr<double> bmanptr = LibUtilities::Foundations::GetBasisManager()[bkey];

    // Finally, we look at the memory footprint of the various classes.
    cout << "Size of points key: " << sizeof(PointsKey) << endl;
    cout << "Size of GGL points: " << sizeof(Points<double, Segment, GaussGaussLegendre>) << endl;
    cout << "Size of GGL quad pts: " << sizeof(Points<double, Quadrilateral, std::tuple<GaussGaussLegendre, GaussGaussLegendre>>) << endl;
    cout << "Size of Modified basis: " << sizeof(Basis<double, Segment, GaussGaussLegendre, ModifiedLegendre>) << endl;

}
