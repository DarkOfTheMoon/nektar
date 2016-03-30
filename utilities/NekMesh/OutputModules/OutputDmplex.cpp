////////////////////////////////////////////////////////////////////////////////
//
//  File0utputDmplex.cpp
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
//  Description: Nektar++ DMPlex output format.
//
////////////////////////////////////////////////////////////////////////////////

#include "OutputDmplex.h"
#include "LibUtilities/Communication/CommMpi.h"
#include "PetscEnv.h"

#include "petscdmplex.h"
#include "petscviewerhdf5.h"
#include <tinyxml.h>

using Nektar::LibUtilities::GetCommFactory;
using Nektar::LibUtilities::CommSharedPtr;
using Nektar::LibUtilities::CommMpiSharedPtr;
using Nektar::LibUtilities::CommMpi;

namespace Nektar
{
namespace Utilities
{
ModuleKey OutputDmplex::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "dmx"),
    OutputDmplex::create,
    "Writes a Nektar++ XML and DMPlex file.");

OutputDmplex::OutputDmplex(MeshSharedPtr m) : OutputXmlBase(m)
{
    // Create any options needed by the module.
}

OutputDmplex::~OutputDmplex()
{
}

#define PCALL(petscFunc, args)                                                 \
    {                                                                          \
        PetscErrorCode ierr = petscFunc args;                                  \
        if (ierr)                                                              \
        {                                                                      \
            cerr << "Error in PETSc call: function = " << #petscFunc           \
                 << ", file = " << __FILE__ << ", line = " << __LINE__         \
                 << ", code = " << ierr << endl;                               \
            abort();                                                           \
        }                                                                      \
    }

void OutputDmplex::Process()
{
    // We'll set up MPI
    CommSharedPtr comm =
        GetCommFactory().CreateInstance("ParallelMPI", 0, NULL);
    CommMpiSharedPtr mpi = boost::dynamic_pointer_cast<CommMpi>(comm);

    // And then the PETSc environment
    {
        PetscEnv penv(NULL, NULL);
        // When this variable goes out of scope, it will shut down petsc
        // Call comm->Finalize() after to deal with MPI.

        if (m_mesh->m_verbose)
        {
            cout << "OutputDmplex: Writing file..." << endl;
        }

        TiXmlDocument doc;
        TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
        doc.LinkEndChild(decl);

        TiXmlElement *root = new TiXmlElement("NEKTAR");
        doc.LinkEndChild(root);

        // Begin <GEOMETRY> section
        TiXmlElement *geomTag = new TiXmlElement("GEOMETRY");
        geomTag->SetAttribute("DIM", m_mesh->m_expDim);
        geomTag->SetAttribute("SPACE", m_mesh->m_spaceDim);
        // TODO: Add the path to the DMPlex file
        root->LinkEndChild(geomTag);

        DM plex;
        PCALL(DMPlexCreate, (mpi->GetComm(), &plex));
        PCALL(DMPlexSetDimension, (plex, m_mesh->m_spaceDim));

        // Now we need to know the size of the "chart", which is the set of
        // all vertices, edges, faces, and cells.
        unsigned vertStart = 0;
        unsigned vertCount = m_mesh->m_vertexSet.size();
        unsigned vertEnd   = vertStart + vertCount;

        unsigned edgeStart = vertEnd;
        unsigned edgeCount = m_mesh->m_edgeSet.size();
        unsigned edgeEnd   = edgeStart + edgeCount;

        unsigned faceStart = edgeEnd;
        unsigned faceCount = m_mesh->m_faceSet.size();
        unsigned faceEnd   = faceStart + faceCount;

        unsigned elemStart = faceEnd;
        unsigned elemCount = m_mesh->m_element[m_mesh->m_expDim].size();
        unsigned elemEnd   = elemStart + elemCount;

        unsigned chartSize = elemEnd;

        // Set this for the plex
        PCALL(DMPlexSetChart, (plex, 0, chartSize));

        // Annoyingly, I think we have to first set all the cone sizes.
        // Start with the elements
        {
            vector<ElementSharedPtr> &elmt =
                m_mesh->m_element[m_mesh->m_expDim];
            std::set<ElementSharedPtr>::iterator it;
            std::set<ElementSharedPtr> tmp(elmt.begin(), elmt.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                ElementSharedPtr el = *it;
                unsigned i          = el->GetId() + elemStart;
                PCALL(DMPlexSetConeSize, (plex, i, el->GetFaceCount()));
            }
        }
        // Now faces
        {
            std::set<FaceSharedPtr>::iterator it;
            std::set<FaceSharedPtr> tmp(m_mesh->m_faceSet.begin(),
                                        m_mesh->m_faceSet.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                FaceSharedPtr fa             = *it;
                unsigned i                   = fa->m_id + faceStart;
                vector<EdgeSharedPtr> &edges = fa->m_edgeList;
                PCALL(DMPlexSetConeSize, (plex, i, edges.size()));
            }
        }
        // Edges next
        {
            std::set<EdgeSharedPtr>::iterator it;
            std::set<EdgeSharedPtr> tmp(m_mesh->m_edgeSet.begin(),
                                        m_mesh->m_edgeSet.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                EdgeSharedPtr ed = *it;
                unsigned i       = ed->m_id + edgeStart;
                PCALL(DMPlexSetConeSize, (plex, i, 2));
            }
        }
        // Now build the internal state...
        PCALL(DMSetUp, (plex));

        // ...and go through setting the actual connections
        // Start with the elements
        {
            vector<ElementSharedPtr> &elmt =
                m_mesh->m_element[m_mesh->m_expDim];
            std::set<ElementSharedPtr>::iterator it;
            std::set<ElementSharedPtr> tmp(elmt.begin(), elmt.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                ElementSharedPtr el = *it;
                unsigned i          = el->GetId() + elemStart;
                PetscInt cone[el->GetFaceCount()];
                for (unsigned j = 0; j < el->GetFaceCount(); ++j)
                    cone[j]     = el->GetFace(j)->m_id + faceStart;

                PCALL(DMPlexSetCone, (plex, i, cone));
            }
        }
        // Now faces
        {
            std::set<FaceSharedPtr>::iterator it;
            std::set<FaceSharedPtr> tmp(m_mesh->m_faceSet.begin(),
                                        m_mesh->m_faceSet.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                FaceSharedPtr fa             = *it;
                unsigned i                   = fa->m_id + faceStart;
                vector<EdgeSharedPtr> &edges = fa->m_edgeList;
                PetscInt cone[edges.size()];
                for (unsigned j = 0; j < edges.size(); ++j)
                    cone[j]     = edges[j]->m_id + edgeStart;
                PCALL(DMPlexSetCone, (plex, i, cone));
            }
        }
        // Edges next
        {
            std::set<EdgeSharedPtr>::iterator it;
            std::set<EdgeSharedPtr> tmp(m_mesh->m_edgeSet.begin(),
                                        m_mesh->m_edgeSet.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                EdgeSharedPtr ed = *it;
                unsigned i       = ed->m_id + edgeStart;
                PetscInt cone[2];
                cone[0] = ed->m_n1->m_id + vertStart;
                cone[1] = ed->m_n2->m_id + vertStart;
                PCALL(DMPlexSetCone, (plex, i, cone));
            }
        }
        // Now set the coordinates
        {
            PetscSection coordSection;
            PCALL(DMGetCoordinateSection, (plex, &coordSection));
            PCALL(PetscSectionSetNumFields, (coordSection, 1));
            unsigned dim = m_mesh->m_spaceDim;
            PCALL(PetscSectionSetFieldComponents, (coordSection, 0, dim));
            PCALL(PetscSectionSetChart, (coordSection, vertStart, vertEnd));
            for (unsigned i = vertStart; i < vertEnd; i++)
            {
                PCALL(PetscSectionSetDof, (coordSection, i, dim));
                PCALL(PetscSectionSetFieldDof, (coordSection, i, 0, dim));
            }
            PCALL(PetscSectionSetUp, (coordSection));
            PetscInt coordSize;
            PCALL(PetscSectionGetStorageSize, (coordSection, &coordSize));
            Vec coordinates;
            PCALL(VecCreate, (mpi->GetComm(), &coordinates));
            PCALL(PetscObjectSetName,
                  ((PetscObject)coordinates, "coordinates"));
            PCALL(VecSetSizes, (coordinates, coordSize, PETSC_DETERMINE));
            PCALL(VecSetType, (coordinates, VECSTANDARD));
            PetscScalar *coords;
            PCALL(VecGetArray, (coordinates, &coords));

            std::set<NodeSharedPtr>::iterator it;
            std::set<NodeSharedPtr> tmp(m_mesh->m_vertexSet.begin(),
                                        m_mesh->m_vertexSet.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                NodeSharedPtr n     = *it;
                unsigned i          = n->m_id + vertStart;
                coords[i * dim + 0] = n->m_x;
                coords[i * dim + 1] = n->m_y;
                coords[i * dim + 2] = n->m_z;
            }

            PCALL(VecRestoreArray, (coordinates, &coords));
            PCALL(DMSetCoordinatesLocal, (plex, coordinates));
            PCALL(VecDestroy, (&coordinates));
        }

        // Now build all the links etc
        PCALL(DMPlexSymmetrize, (plex));
        PCALL(DMPlexStratify, (plex));

        {
            // Extract the output filename and extension
            string filename = m_config["outfile"].as<string>();
            size_t iDot     = filename.rfind('.');
            filename.replace(iDot + 1, 3, "dmp");

            // In PETSc-ese a viewer deals with IO
            PetscViewer viewer;
            // Create an HDF5 one
            //                     PCALL(PetscViewerBinaryOpen,
            //                            (mpi->GetComm(),
            //                            m_config["outfile"].as<string>().c_str(),
            //                            FILE_MODE_WRITE, &viewer));
            // HDF issues - try text
            PCALL(PetscViewerASCIIOpen,
                  (mpi->GetComm(), filename.c_str(), &viewer));
            // Write the plex out
            PCALL(PetscObjectView, ((PetscObject)plex, viewer));
            // Close the file
            PCALL(PetscViewerDestroy, (&viewer));
        }

        // Serialise the other bits of the mesh to the XML
        WriteXmlCurves(geomTag);
        WriteXmlComposites(geomTag);
        WriteXmlDomain(geomTag);

        // Extract the output filename and extension
        doc.SaveFile(m_config["outfile"].as<string>());

    } // Petsc shutdown here.
    comm->Finalise();
}
}
}
