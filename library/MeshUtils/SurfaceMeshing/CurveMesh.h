////////////////////////////////////////////////////////////////////////////////
//
//  File: Curvemesh.h
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
//  Description: object for individual curve meshes.
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_MESHUTILS_SURFACEMESHING_CURVEMESH_H
#define NEKTAR_MESHUTILS_SURFACEMESHING_CURVEMESH_H

#include <boost/shared_ptr.hpp>

#include <MeshUtils/MeshElem.hpp>
#include <LibUtilities/CADSystem/CADCurve.h>
#include <MeshUtils/Octree/Octree.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar {
namespace MeshUtils {

class CurveMesh
{
    public:
        friend class MemoryManager<CurveMesh>;

        /**
         * @brief default constructor
         */
        CurveMesh(bool verbose,
                  int id,
                  const LibUtilities::CADCurveSharedPtr &cad,
                  const OctreeSharedPtr &oct) :
                  m_cadcurve(cad),m_octree(oct),m_verbose(verbose),
                  m_id(id)
        {
        };

        /**
         * @brief execute meshing
         */
        void Mesh(std::map<int, MeshNodeSharedPtr> &Nodes,
                  std::map<int, MeshEdgeSharedPtr> &Edges);

        /**
         * @brief get id of first node
         */
        int GetFirstPoint(){return m_meshpoints[0];}

        /**
         * @brief get id of last node
         */
        int GetLastPoint(){return m_meshpoints.back();}

        /**
         * @brief get list of mesh nodes
         */
        std::vector<int> GetMeshPoints(){return m_meshpoints;}

        /**
         * @brief get the number of points in the curve
         */
        int GetNumPoints(){return m_meshpoints.size();}

        /**
         * @brief print report to screen
         */
        void Report();

        int SplitEdge(int a, int b,
                       std::map<int, MeshNodeSharedPtr> &Nodes,
                       std::map<int, MeshEdgeSharedPtr> &Edges);

    private:

        /**
         * @brief get node spacing sampling function
         */
        void GetSampleFunction();

        /**
         * @brief get node spacing phi function
         */
        void GetPhiFunction();

        /**
         * @brief evaluate paramter ds at curve location s
         */
        NekDouble EvaluateDS(NekDouble s);

        /**
         * @brief evaluate paramter ps at curve location s
         */
        NekDouble EvaluatePS(NekDouble s);

        /// CAD curve
        LibUtilities::CADCurveSharedPtr m_cadcurve;
        /// Octree object
        OctreeSharedPtr m_octree;
        /// length of the curve in real space
        NekDouble m_curvelength;
        /// number of sampling points used in algorithm
        int m_numSamplePoints;
        /// coords of the ends of the parametric curve
        Array<OneD, NekDouble> m_bounds;
        /// array of function ds evaluations
        std::vector<std::vector<NekDouble> > m_dst;
        /// array of function ps evaluations
        std::vector<std::vector<NekDouble> > m_ps;
        /// spacing function evaluation
        NekDouble Ae;
        /// ds
        NekDouble ds;
        /// number of edges to be made in the curve as defined by the spacing funtion
        int Ne;
        /// paramteric coordiates of the mesh nodes
        std::vector<NekDouble> meshsvalue;
        /// ids of the mesh nodes
        std::vector<int> m_meshpoints;
        /// verbosity
        bool m_verbose;
        /// id of the curvemesh
        int m_id;
};

typedef boost::shared_ptr<CurveMesh> CurveMeshSharedPtr;

}
}

#endif
