///////////////////////////////////////////////////////////////////////////////
//
// File: PointsTraitsStaticInit.cpp
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
// Description: Static initialisation of point traits
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Points/PointsTraits.hpp>

namespace Nektar
{
namespace LibUtilities
{
namespace Foundations
{

PointsStringToIdMap& GetPointsStringToIdMap()
{
    typedef Loki::SingletonHolder<PointsStringToIdMap,
        Loki::CreateUsingNew,
        Loki::NoDestroy> Type;
    return Type::Instance();
}

PointsIdToStringMap& GetPointsIdToStringMap()
{
    typedef Loki::SingletonHolder<PointsIdToStringMap,
        Loki::CreateUsingNew,
        Loki::NoDestroy> Type;
    return Type::Instance();
}

namespace traits
{
const std::string points_traits<GaussGaussLegendre>::name = GetPointsStringToIdMap().insert(std::pair<std::string, int>("GaussGaussLegendre", 1)).first->first;
const int points_traits<GaussGaussLegendre>::id = GetPointsIdToStringMap().insert(std::pair<int, std::string>(1, "GaussGaussLegendre")).first->first;
}

}
}
}
