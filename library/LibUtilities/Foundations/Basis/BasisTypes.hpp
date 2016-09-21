///////////////////////////////////////////////////////////////////////////////
//
// File: IProductWRTBase.cpp
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
// Description: IProductWRTBase operator implementations
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBUTILITIES_FOUNDATIONS_BASISTYPES
#define LIBUTILITIES_FOUNDATIONS_BASISTYPES

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/ShapeTypes.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

struct ModifiedLegendre {};
struct Lagrange {};
struct BernsteinTriangle {};
struct BernsteinTetrahedron {};

namespace traits
{

// Default values for all types
struct basis_traits_base {
    static const unsigned int dimension = 0;
    static const unsigned int num_constituent_bases = 1;
    static const bool is_boundary_interior = false;
    static const bool is_collocation = false;
    static const unsigned int get_total_modes(unsigned int nmodes) {
        return nmodes;
    }
    typedef Segment native_shape;
};


// Primary template
template<typename... TBasis>
struct basis_traits : public basis_traits_base {
    static const unsigned int num_constituent_bases = sizeof...(TBasis);
};

template<typename T1, typename... TBasis>
struct basis_traits<T1, TBasis...> : public basis_traits_base {
    static const unsigned int dimension = basis_traits<T1>::dimension + basis_traits<TBasis...>::dimension;
    static_assert(dimension <= 3, "Composite basis dimension > 3.");
    static const bool is_boundary_interior = basis_traits<T1>::is_boundary_interior && basis_traits<TBasis...>::is_boundary_interior;
    static const bool is_collocation = basis_traits<T1>::is_collocation && basis_traits<TBasis...>::is_collocation;
};

template<>
struct basis_traits<ModifiedLegendre> : public basis_traits_base
{
    static const unsigned int dimension = 1;
    static const bool is_boundary_interior = true;
};

template<>
struct basis_traits<Lagrange> : public basis_traits_base
{
    static const unsigned int dimension = 1;
    static const bool collocation = true;
};

template<>
struct basis_traits<BernsteinTriangle> : public basis_traits_base
{
    static const unsigned int dimension = 2;
    typedef Triangle native_shape;
};

template<>
struct basis_traits<BernsteinTetrahedron> : public basis_traits_base
{
    static const unsigned int dimension = 3;
    typedef Tetrahedron native_shape;
};


template<typename TShape, typename TBasis>
struct expansion_traits_base
{
    static const bool is_valid = false;
    static const bool is_boundary_interior = std::is_same<TBasis, ModifiedLegendre>::value;
    static const int get_num_coefficients(int Na, int Nb = 1, int Nc = 1)
    {
        return 0;
    }
    static const int get_num_bnd_coefficients(int Na, int Nb = 1, int Nc = 1)
    {
        return 0;
    }
};

template<typename TShape, typename TBasis>
struct expansion_traits : public shape_traits<TShape>,
                             public basis_traits<TBasis>,
                             public expansion_traits_base<TShape, TBasis>
{
};

template<typename TBasis>
struct expansion_traits<Segment, TBasis> :
                             public shape_traits<Segment>,
                             public basis_traits<TBasis>,
                             public expansion_traits_base<Segment, TBasis>
{
    static const bool is_valid = true;
    static const int get_num_coefficients(int Na)
    {
        return Na;
    }
    static const int get_num_bnd_coefficients(int Na)
    {
        return 2;
    }
};


template<typename TBasis>
struct expansion_traits<Quadrilateral, TBasis> :
                             public shape_traits<Quadrilateral>,
                             public basis_traits<TBasis>,
                             public expansion_traits_base<Quadrilateral, TBasis>
{
    static const bool is_valid = true;
    static const int get_num_coefficients(int Na, int Nb)
    {
        ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
        return Na*Nb;
    }
    static const int get_num_bnd_coefficients(int Na, int Nb)
    {
        ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
        return 2*(Na-1) + 2*(Nb-1);
    }
};


template<typename TBasis>
struct expansion_traits<Triangle, TBasis> :
                             public shape_traits<Triangle>,
                             public basis_traits<TBasis>,
                             public expansion_traits_base<Triangle, TBasis>
{
    static const bool is_valid = false;
    static const int get_num_coefficients(int Na, int Nb)
    {
        ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
        ASSERTL1(Na <= Nb, "order in 'a' direction is higher "
                 "than order in 'b' direction");
        return Na*(Na+1)/2 + Na*(Nb-Na);
    }
    static const int get_num_bnd_coefficients(int Na, int Nb)
    {
        ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
        ASSERTL1(Na <= Nb, "order in 'a' direction is higher "
                 "than order in 'b' direction");
        return (Na-1) + 2*(Nb-1);
    }
};


template<typename TBasis>
struct expansion_traits<Hexahedron, TBasis> :
                             public shape_traits<Hexahedron>,
                             public basis_traits<TBasis>,
                             public expansion_traits_base<Hexahedron, TBasis>
{
    static const bool is_valid = true;
    static const int get_num_coefficients(int Na, int Nb, int Nc)
    {
        ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
        ASSERTL2(Nc > 1, "Order in 'c' direction must be > 1.");
        return Na*Nb*Nc;
    }
    static const int get_num_bnd_coefficients(int Na, int Nb, int Nc)
    {
        ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
        ASSERTL2(Nc > 1, "Order in 'c' direction must be > 1.");
        return 2*Na*Nb + 2*Na*Nc + 2*Nb*Nc
                - 4*(Na + Nb + Nc) + 8;
    }
};


template<typename TBasis>
struct expansion_traits<Tetrahedron, TBasis> :
                             public shape_traits<Tetrahedron>,
                             public basis_traits<TBasis>,
                             public expansion_traits_base<Tetrahedron, TBasis>
{
    static const bool is_valid = true;

    /**
     * Adds up the number of cells in a truncated Nc by Nc by Nc
     * pyramid, where the longest Na rows and longest Nb columns are
     * kept. Example: (Na, Nb, Nc) = (3, 4, 5); The number of
     * coefficients is the sum of the elements of the following
     * matrix:
     *
     * |5  4  3  2  0|
     * |4  3  2  0   |
     * |3  2  0      |
     * |0  0         |
     * |0            |
     *
     * Sum = 28 = number of tet coefficients.
     */
    static const int get_num_coefficients(int Na, int Nb, int Nc)
    {
        ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
        ASSERTL2(Nc > 1, "Order in 'c' direction must be > 1.");
        ASSERTL1(Na <= Nc, "order in 'a' direction is higher "
                 "than order in 'c' direction");
        ASSERTL1(Nb <= Nc, "order in 'b' direction is higher "
                 "than order in 'c' direction");
        int nCoef = 0;
        for (int a = 0; a < Na; ++a)
        {
            for (int b = 0; b < Nb - a; ++b)
            {
                for (int c = 0; c < Nc - a - b; ++c)
                {
                    ++nCoef;
                }
            }
        }
        return nCoef;
    }
    static const int get_num_bnd_coefficients(int Na, int Nb, int Nc)
    {
        ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
        ASSERTL2(Nc > 1, "Order in 'c' direction must be > 1.");
        ASSERTL1(Na <= Nc, "order in 'a' direction is higher "
                 "than order in 'c' direction");
        ASSERTL1(Nb <= Nc, "order in 'b' direction is higher "
                 "than order in 'c' direction");

        int nCoef =    Na*(Na+1)/2 + (Nb-Na)*Na // base
                  +    Na*(Na+1)/2 + (Nc-Na)*Na // front
                  + 2*(Nb*(Nb+1)/2 + (Nc-Nb)*Nb)// 2 other sides
                  - Na - 2*Nb - 3*Nc            // less edges
                  + 4;                          // plus vertices

        return nCoef;
    }
};


template<typename TBasis>
struct expansion_traits<Prism, TBasis> :
                             public shape_traits<Prism>,
                             public basis_traits<TBasis>,
                             public expansion_traits_base<Prism, TBasis>
{
    static const bool is_valid = true;
    static const int get_num_coefficients(int Na, int Nb, int Nc)
    {
        ASSERTL1(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL1(Nb > 1, "Order in 'b' direction must be > 1.");
        ASSERTL1(Nc > 1, "Order in 'c' direction must be > 1.");
        ASSERTL1(Na <= Nc, "Order in 'a' direction is higher "
                "than order in 'c' direction.");

        return Nb*traits::expansion_traits<Triangle, TBasis>::get_num_coefficients(Na,Nc);
    }
    static const int get_num_bnd_coefficients(int Na, int Nb, int Nc)
    {
        ASSERTL1(Na > 1, "Order in 'a' direction must be > 1.");
        ASSERTL1(Nb > 1, "Order in 'b' direction must be > 1.");
        ASSERTL1(Nc > 1, "Order in 'c' direction must be > 1.");
        ASSERTL1(Na <= Nc, "Order in 'a' direction is higher "
                "than order in 'c' direction.");

        return Na*Nb + 2*Nb*Nc              // rect faces
           + 2*( Na*(Na+1)/2 + (Nc-Na)*Na ) // tri faces
           - 2*Na - 3*Nb - 4*Nc             // less edges
           + 6;                             // plus vertices
    }
};


template<typename TBasis>
struct expansion_traits<Pyramid, TBasis> :
                             public shape_traits<Pyramid>,
                             public basis_traits<TBasis>,
                             public expansion_traits_base<Pyramid, TBasis>
{
    static const bool is_valid = true;
    static const int get_num_coefficients(int Na, int Nb, int Nc)
    {
        return 0;
    }
    static const int get_num_bnd_coefficients(int Na, int Nb, int Nc)
    {
        return 0;
    }
};


// Further specialisations
template<>
struct expansion_traits<Triangle, BernsteinTriangle> :
                            public shape_traits<Triangle>,
                            public basis_traits<BernsteinTriangle>,
                            public expansion_traits_base<Triangle, BernsteinTriangle>
{
    static const bool is_valid = true;
};

}
}
}
}

#endif
