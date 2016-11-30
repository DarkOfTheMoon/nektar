////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSurf.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: CAD object surface.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_CADSYSTEM_OCE_CADSURFOCE
#define NekMeshUtils_CADSYSTEM_OCE_CADSURFOCE

#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/CADSystem/OCE/OpenCascade.h>

namespace Nektar
{
namespace NekMeshUtils
{

class CADSurfOCE : public CADSurf
{
public:

    static CADSurfSharedPtr create()
    {
        return MemoryManager<CADSurfOCE>::AllocateSharedPtr();
    }

    static std::string key;

    CADSurfOCE()
    {
    }

    ~CADSurfOCE()
    {
    }

    void Initialise(int i, TopoDS_Shape in, std::vector<EdgeLoop> ein);

    virtual Array<OneD, NekDouble> GetBounds();
    virtual Array<OneD, NekDouble> N    (Array<OneD, NekDouble> uv);
    virtual Array<OneD, NekDouble> D1   (Array<OneD, NekDouble> uv);
    virtual Array<OneD, NekDouble> D2   (Array<OneD, NekDouble> uv);
    virtual Array<OneD, NekDouble> P    (Array<OneD, NekDouble> uv);
    virtual Array<OneD, NekDouble> locuv(Array<OneD, NekDouble> p);
    virtual NekDouble DistanceTo(Array<OneD, NekDouble> p);
    virtual void ProjectTo(Array<OneD, NekDouble> &tp,
                           Array<OneD, NekDouble> &uv);
    virtual NekDouble Curvature(Array<OneD, NekDouble> uv);

private:
    /// Function which tests the the value of uv used is within the surface
    void Test(Array<OneD, NekDouble> uv);
    /// OpenCascade object for surface.
    BRepAdaptor_Surface m_occSurface;
    /// Alternate OpenCascade object for surface. Used by reverse lookup.
    Handle(Geom_Surface) m_s;
};

typedef boost::shared_ptr<CADSurf> CADSurfSharedPtr;

}
}

#endif