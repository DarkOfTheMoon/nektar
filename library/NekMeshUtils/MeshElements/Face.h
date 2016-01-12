////////////////////////////////////////////////////////////////////////////////
//
//  File: Face.h
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
//  Description: Mesh manipulation objects.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_MESHELEMENTS_FACE
#define NEKMESHUTILS_MESHELEMENTS_FACE

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>

namespace Nektar
{
namespace NekMeshUtils
{
    /**
     * @brief Represents a face comprised of three or more edges.
     *
     * A face is defined by a list of vertices, a list of edges joining
     * these vertices, and a list of control nodes within the interior of
     * the face, defining the shape of the face.
     */
    class Face {
    public:
        /// Create a new face.
        Face(std::vector<NodeSharedPtr> pVertexList,
             std::vector<NodeSharedPtr> pFaceNodes,
             std::vector<EdgeSharedPtr> pEdgeList,
             LibUtilities::PointsType   pCurveType)
            : m_vertexList(pVertexList),
            m_edgeList  (pEdgeList),
            m_faceNodes (pFaceNodes),
            m_curveType (pCurveType),
             m_geom    () {}

        /// Copy an existing face.
        Face(const Face& pSrc)
            : m_vertexList(pSrc.m_vertexList), m_edgeList  (pSrc.m_edgeList),
            m_faceNodes (pSrc.m_faceNodes),  m_curveType (pSrc.m_curveType),
              m_geom    (pSrc.m_geom) {}
        ~Face() {}

        /// Equality is defined by matching all vertices.
        bool operator==(Face& pSrc)
        {
            std::vector<NodeSharedPtr>::iterator it1, it2;
            for (it1 = m_vertexList.begin(); it1 != m_vertexList.end(); ++it1)
            {
                if (find(pSrc.m_vertexList.begin(), pSrc.m_vertexList.end(), *it1)
                        == pSrc.m_vertexList.end())
                {
                    return false;
                }
            }
            return true;
        }

        /// Returns the total number of nodes (vertices, edge nodes and
        /// face nodes).
        unsigned int GetNodeCount() const
        {
            unsigned int n = m_faceNodes.size();
            for (int i = 0; i < m_edgeList.size(); ++i)
            {
                n += m_edgeList[i]->GetNodeCount();
            }
            n -= m_vertexList.size();
            return n;
        }

        /// Generates a string listing the coordinates of all nodes
        /// associated with this face.
        std::string GetXmlCurveString() const
        {
            std::stringstream s;
            std::string str;

            // Treat 2D point distributions differently to 3D.
            if (m_curveType == LibUtilities::eNodalTriFekete       ||
                m_curveType == LibUtilities::eNodalTriEvenlySpaced ||
                m_curveType == LibUtilities::eNodalTriElec)
            {
                vector<NodeSharedPtr> tmp;
                int n = m_edgeList[0]->GetNodeCount();

                tmp.insert(tmp.end(), m_vertexList.begin(), m_vertexList.end());
                for (int k = 0; k < m_edgeList.size(); ++k)
                {
                    tmp.insert(tmp.end(), m_edgeList[k]->m_edgeNodes.begin(),
                               m_edgeList[k]->m_edgeNodes.end());
                    if (m_edgeList[k]->m_n1 != m_vertexList[k])
                    {
                        // If edge orientation is reversed relative to node
                        // ordering, we need to reverse order of nodes.
                        std::reverse(tmp.begin() + 3 + k*(n-2),
                                     tmp.begin() + 3 + (k+1)*(n-2));
                    }
                }
                tmp.insert(tmp.end(), m_faceNodes.begin(), m_faceNodes.end());

                for (int k = 0; k < tmp.size(); ++k) {
                    s << std::scientific << std::setprecision(8) << "    "
                      <<  tmp[k]->m_x << "  " << tmp[k]->m_y
                      << "  " << tmp[k]->m_z << "    ";
                }

                return s.str();
            }
            else
            {
                // Write out in 2D tensor product order.
                ASSERTL0(m_vertexList.size() == 4,
                         "Face nodes of tensor product only supported "
                         "for quadrilaterals.");

                int n = (int)sqrt((NekDouble)GetNodeCount());
                vector<NodeSharedPtr> tmp(n*n);

                ASSERTL0(n*n == GetNodeCount(), "Wrong number of modes?");

                // Write vertices
                tmp[0]       = m_vertexList[0];
                tmp[n-1]     = m_vertexList[1];
                tmp[n*n-1]   = m_vertexList[2];
                tmp[n*(n-1)] = m_vertexList[3];

                // Write edge-interior
                int skips[4][2] = {{0,1}, {n-1,n}, {n*n-1,-1}, {n*(n-1),-n}};
                for (int i = 0; i < 4; ++i)
                {
                    bool reverseEdge = m_edgeList[i]->m_n1 == m_vertexList[i];

                    if (!reverseEdge)
                    {
                        for (int j = 1; j < n-1; ++j)
                        {
                            tmp[skips[i][0] + j*skips[i][1]] =
                                m_edgeList[i]->m_edgeNodes[n-2-j];
                        }
                    }
                    else
                    {
                        for (int j = 1; j < n-1; ++j)
                        {
                            tmp[skips[i][0] + j*skips[i][1]] =
                                m_edgeList[i]->m_edgeNodes[j-1];
                        }
                    }
                }

                // Write interior
                for (int i = 1; i < n-1; ++i)
                {
                    for (int j = 1; j < n-1; ++j)
                    {
                        tmp[i*n+j] = m_faceNodes[(i-1)*(n-2)+(j-1)];
                    }
                }

                for (int k = 0; k < tmp.size(); ++k) {
                    s << std::scientific << std::setprecision(8) << "    "
                      <<  tmp[k]->m_x << "  " << tmp[k]->m_y
                      << "  " << tmp[k]->m_z << "    ";
                }

                return s.str();
            }
        }

        /// Generate either SpatialDomains::TriGeom or
        /// SpatialDomains::QuadGeom for this element.
        SpatialDomains::Geometry2DSharedPtr GetGeom(int coordDim)
        {
            int nEdge = m_edgeList.size();

            SpatialDomains::SegGeomSharedPtr    edges[4];
            SpatialDomains::Geometry2DSharedPtr ret;
            StdRegions::Orientation             edgeo[4];

            for (int i = 0; i < nEdge; ++i)
            {
                edges[i] = m_edgeList[i]->GetGeom(coordDim);
            }

            for (int i = 0; i < nEdge; ++i)
            {
                edgeo[i] = m_edgeList[i]->m_n1 == m_vertexList[i] ?
                    StdRegions::eForwards : StdRegions::eBackwards;
            }

            if(m_faceNodes.size() > 0)
            {
                if (nEdge == 3)
                {
                    SpatialDomains::CurveSharedPtr c =
                        MemoryManager<SpatialDomains::Curve>::
                        AllocateSharedPtr(m_id, m_curveType);

                    for(int j = 0; j < m_vertexList.size(); j++)
                    {
                        c->m_points.push_back(m_vertexList[j]->GetGeom(coordDim));
                    }
                    for(int j = 0; j < nEdge; j++)
                    {
                        vector<NodeSharedPtr> ed = m_edgeList[j]->m_edgeNodes;
                        if(edgeo[j] == StdRegions::eBackwards)
                        {
                            for(int k = ed.size()-1; k >=0; k--)
                            {
                                c->m_points.push_back(ed[k]->GetGeom(coordDim));
                            }
                        }
                        else
                        {
                            for(int k = 0; k < ed.size(); k++)
                            {
                                c->m_points.push_back(ed[k]->GetGeom(coordDim));
                            }
                        }
                    }
                    for(int j = 0; j < m_faceNodes.size(); j++)
                    {
                        c->m_points.push_back(m_faceNodes[j]->GetGeom(coordDim));
                    }

                    ret = MemoryManager<SpatialDomains::TriGeom>::
                        AllocateSharedPtr(m_id, edges, edgeo, c);
                }
                else
                {
                    SpatialDomains::CurveSharedPtr c =
                        MemoryManager<SpatialDomains::Curve>::
                        AllocateSharedPtr(m_id, m_curveType);

                    ASSERTL0(m_vertexList.size() == 4,
                             "Face nodes of tensor product only supported "
                             "for quadrilaterals.");

                    int n = (int)sqrt((NekDouble)GetNodeCount());
                    vector<NodeSharedPtr> tmp(n*n);

                    ASSERTL0(n*n == GetNodeCount(), "Wrong number of modes?");

                    // Write vertices
                    tmp[0]       = m_vertexList[0];
                    tmp[n-1]     = m_vertexList[1];
                    tmp[n*n-1]   = m_vertexList[2];
                    tmp[n*(n-1)] = m_vertexList[3];

                    // Write edge-interior
                    int skips[4][2] = {{0,1}, {n-1,n}, {n*n-1,-1}, {n*(n-1),-n}};
                    for (int i = 0; i < 4; ++i)
                    {
                        bool reverseEdge = edgeo[i] == StdRegions::eBackwards;

                        if (reverseEdge)
                        {
                            for (int j = 1; j < n-1; ++j)
                            {
                                tmp[skips[i][0] + j*skips[i][1]] =
                                    m_edgeList[i]->m_edgeNodes[n-2-j];
                            }
                        }
                        else
                        {
                            for (int j = 1; j < n-1; ++j)
                            {
                                tmp[skips[i][0] + j*skips[i][1]] =
                                    m_edgeList[i]->m_edgeNodes[j-1];
                            }
                        }
                    }

                    // Write interior
                    for (int i = 1; i < n-1; ++i)
                    {
                        for (int j = 1; j < n-1; ++j)
                        {
                            tmp[i*n+j] = m_faceNodes[(i-1)*(n-2)+(j-1)];
                        }
                    }

                    for (int k = 0; k < tmp.size(); ++k)
                    {
                        c->m_points.push_back(tmp[k]->GetGeom(coordDim));
                    }

                    ret = MemoryManager<SpatialDomains::QuadGeom>::
                        AllocateSharedPtr(m_id, edges, edgeo, c);
                }
            }
            else
            {
                if (nEdge == 3)
                {
                    ret = MemoryManager<SpatialDomains::TriGeom>::
                        AllocateSharedPtr(m_id, edges, edgeo);
                }
                else
                {
                    ret = MemoryManager<SpatialDomains::QuadGeom>::
                        AllocateSharedPtr(m_id, edges, edgeo);
                }
            }

            return ret;
        }

        /// ID of the face.
        unsigned int m_id;
        /// List of vertex nodes.
        std::vector<NodeSharedPtr> m_vertexList;
        /// List of corresponding edges.
        std::vector<EdgeSharedPtr> m_edgeList;
        /// List of face-interior nodes defining the shape of the face.
        std::vector<NodeSharedPtr> m_faceNodes;
        /// Distribution of points in this face.
        LibUtilities::PointsType   m_curveType;
        /// Element(s) which are linked to this face.
        vector<pair<ElementSharedPtr, int> > m_elLink;

        SpatialDomains::Geometry2DSharedPtr m_geom;
    };
    /// Shared pointer to a face.
    typedef boost::shared_ptr<Face> FaceSharedPtr;

    bool operator==(FaceSharedPtr const &p1, FaceSharedPtr const &p2);
    bool operator< (FaceSharedPtr const &p1, FaceSharedPtr const &p2);

    struct FaceHash : std::unary_function<FaceSharedPtr, std::size_t>
    {
        std::size_t operator()(FaceSharedPtr const& p) const
        {
            unsigned int              nVert = p->m_vertexList.size();
            std::size_t               seed  = 0;
            std::vector<unsigned int> ids(nVert);

            for (int i = 0; i < nVert; ++i)
            {
                ids[i] = p->m_vertexList[i]->m_id;
            }

            std::sort(ids.begin(), ids.end());
            boost::hash_range(seed, ids.begin(), ids.end());

            return seed;
        }
    };
    typedef boost::unordered_set<FaceSharedPtr, FaceHash> FaceSet;

}
}

#endif
