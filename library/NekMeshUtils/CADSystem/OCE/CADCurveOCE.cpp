////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurve.cpp
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
//  Description: CAD object curve methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/CADSystem/OCE/CADCurveOCE.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADCurveOCE::key = GetCADCurveFactory().RegisterCreatorFunction(
    "oce", CADCurveOCE::create, "CADCurveOCE");

NekDouble CADCurveOCE::tAtArcLength(NekDouble s)
{
    NekDouble dt =
        (m_occCurve.LastParameter() - m_occCurve.FirstParameter()) / (5000);
    NekDouble t = m_occCurve.FirstParameter();

    NekDouble len = 0.0;

    while (len <= s)
    {
        gp_Pnt P1, P2;
        gp_Vec drdt1, drdt2;

        m_occCurve.D1(t, P1, drdt1);
        t += dt;
        m_occCurve.D1(t, P2, drdt2);

        len += (drdt1.Magnitude() + drdt2.Magnitude()) / 2.0 * dt;
    }

    return t - dt;
}

NekDouble CADCurveOCE::Length(NekDouble ti, NekDouble tf)
{
    Array<OneD, NekDouble> b = Bounds();
    Handle(Geom_Curve) NewCurve = new Geom_TrimmedCurve(m_c, ti, tf);
    TopoDS_Edge NewEdge = BRepBuilderAPI_MakeEdge(NewCurve);
    GProp_GProps System;
    BRepGProp::LinearProperties(NewEdge, System);
    return System.Mass() / 1000.0;
}

NekDouble CADCurveOCE::loct(Array<OneD, NekDouble> xyz)
{
    NekDouble t = 0.0;
    Array<OneD, NekDouble> b = Bounds();

    gp_Pnt loc(xyz[0]*1000.0, xyz[1]*1000.0, xyz[2]*1000.0);

    ShapeAnalysis_Curve sac;
    gp_Pnt p;
    sac.Project(m_c,loc,1e-8 ,p,t);

    if(p.Distance(loc) > 1e-5)
    {
        cerr << "large loct distance" << endl;
    }
    return t;
}

Array<OneD, NekDouble> CADCurveOCE::P(NekDouble t)
{
    Array<OneD, NekDouble> location(3);
    gp_Pnt loc = m_occCurve.Value(t);

    location[0] = loc.X();
    location[1] = loc.Y();
    location[2] = loc.Z();

    return location;
}

Array<OneD, NekDouble> CADCurveOCE::D2(NekDouble t)
{
    Array<OneD, NekDouble> out(9);
    gp_Pnt loc;
    gp_Vec d1, d2;
    m_occCurve.D2(t, loc, d1, d2);

    out[0] = loc.X();
    out[1] = loc.Y();
    out[2] = loc.Z();
    out[3] = d1.X();
    out[4] = d1.Y();
    out[5] = d1.Z();
    out[6] = d2.X();
    out[7] = d2.Y();
    out[8] = d2.Z();

    return out;
}

Array<OneD, NekDouble> CADCurveOCE::Bounds()
{
    Array<OneD, NekDouble> t(2);
    t[0] = m_occCurve.FirstParameter();
    t[1] = m_occCurve.LastParameter();

    return t;
}

Array<OneD, NekDouble> CADCurveOCE::GetMinMax()
{
    Array<OneD, NekDouble> locs(6);

    gp_Pnt start =
        BRep_Tool::Pnt(TopExp::FirstVertex(m_occEdge, Standard_True));
    gp_Pnt end = BRep_Tool::Pnt(TopExp::LastVertex(m_occEdge, Standard_True));

    locs[0] = start.X();
    locs[1] = start.Y();
    locs[2] = start.Z();
    locs[3] = end.X();
    locs[4] = end.Y();
    locs[5] = end.Z();

    return locs;
}

}
}