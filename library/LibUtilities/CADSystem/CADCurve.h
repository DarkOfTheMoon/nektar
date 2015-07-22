////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurve.h
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
//  Description: CAD object curve.
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_LIB_UTILITIES_CADSYSTEM_CADCURVE_H
#define NEKTAR_LIB_UTILITIES_CADSYSTEM_CADCURVE_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar {
namespace LibUtilities {

    /**
	 * @brief class for CAD curves.
	 *
	 * This class wraps the opencascade BRepAdaptor_Curve class for use with
     * nektar++
	 */

    class CADCurve
    {
    public:
        friend class MemoryManager<CADCurve>;

        /**
		 * @brief Defualt constructor.
		 */
        CADCurve(int i, TopoDS_Shape in);

        void Bounds(Array<OneD, NekDouble> &out);
        NekDouble Length(NekDouble ti, NekDouble tf);
        void P(NekDouble t, Array<OneD, NekDouble> &out);
        NekDouble tAtArcLength(NekDouble s);
        int GetID(){return ID;}
        void SetAdjSurf(std::vector<int> i){adjSurfs=i;}
        std::vector<int> GetAdjSurf(){return adjSurfs;}
        void GetMinMax(gp_Pnt &start, gp_Pnt &end);

    private:
        
        int ID;
        BRepAdaptor_Curve occCurve;
        std::vector<int> adjSurfs;


    };

    typedef boost::shared_ptr<CADCurve> CADCurveSharedPtr;


}
}

#endif
