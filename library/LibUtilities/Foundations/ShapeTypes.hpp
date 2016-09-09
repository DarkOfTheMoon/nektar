///////////////////////////////////////////////////////////////////////////////
//
// File: ShapeTraits.hpp
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
// Description: Traits for shape types
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBUTILITIES_FOUNDATIONS_SHAPETYPES
#define LIBUTILITIES_FOUNDATIONS_SHAPETYPES

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

/**
 * @brief Types of elements.
 */
enum ShapeType
{
    eNoShapeType,
    ePoint,
    eSegment,
    eTriangle,
    eQuadrilateral,
    eTetrahedron,
    ePyramid,
    ePrism,
    eHexahedron,
    SIZE_ShapeType
};

struct Point {};
struct Segment {};
struct Triangle {};
struct Quadrilateral {};
struct Tetrahedron {};
struct Pyramid {};
struct Prism {};
struct Hexahedron {};

namespace traits
{

struct shape_traits_base {
        static const int dimension = 0;
        static const bool is_simplex = false;
        static constexpr char name[] = "";
        static const enum ShapeType type = eNoShapeType;
};

// Primary shape traits template
template<typename TShape>
struct shape_traits : public shape_traits_base {
};

template<>
struct shape_traits<Segment> : public shape_traits_base {
        static const int dimension = 1;
        static const bool is_simplex = false;
        static constexpr char name[] = "Segment";
        static const enum ShapeType type = eSegment;
};

template<>
struct shape_traits<Triangle> : public shape_traits_base {
        static const int dimension = 2;
        static const bool is_simplex = true;
        static constexpr char name[] = "Triangle";
        static const enum ShapeType type = eTriangle;
};

template<>
struct shape_traits<Quadrilateral> : public shape_traits_base {
        static const int dimension = 2;
        static const bool is_simplex = false;
        static constexpr char name[] = "Quadrilateral";
        static const enum ShapeType type = eQuadrilateral;
};

template<>
struct shape_traits<Tetrahedron> : public shape_traits_base {
        static const int dimension = 3;
        static const bool is_simplex = true;
        static constexpr char name[] = "Tetrahedron";
        static const enum ShapeType type = eTetrahedron;
};

template<>
struct shape_traits<Prism> : public shape_traits_base {
        static const int dimension = 3;
        static const bool is_simplex = true;
        static constexpr char name[] = "Prism";
        static const enum ShapeType type = ePrism;
};

template<>
struct shape_traits<Pyramid> : public shape_traits_base {
        static const int dimension = 3;
        static const bool is_simplex = true;
        static constexpr char name[] = "Pyramid";
        static const enum ShapeType type = ePyramid;
};

template<>
struct shape_traits<Hexahedron> : public shape_traits_base {
        static const int dimension = 3;
        static const bool is_simplex = false;
        static constexpr char name[] = "Hexahedron";
        static const enum ShapeType type = eHexahedron;
};

}


}
}
}

#endif
