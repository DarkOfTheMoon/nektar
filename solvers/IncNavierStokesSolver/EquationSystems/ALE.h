///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionTerm.h
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
// Description: Base class for Navier-Stokes advection term
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_SOLVERS_ALE_H
#define NEKTAR_SOLVERS_ALE_H

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Spline.h>
#include <LibUtilities/Communication/Comm.h>
#include <SolverUtils/EquationSystem.h>
#include <MultiRegions/ExpList2D.h>
//#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/MeshGraph2D.h>
//#include <SpatialDomains/HistoryPoints.h>
//#include <SpatialDomains/SpatialData.h>

//#include <MultiRegions/ContField1D.h>
//#include <MultiRegions/ContField2D.h>
//#include <MultiRegions/ContField3D.h>

//#include <MultiRegions/DisContField1D.h>
//#include <MultiRegions/DisContField2D.h>
//#include <Auxiliary/ADRBase.h>
//#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

namespace Nektar
{

    class ALE: public SolverUtils::EquationSystem
    {
    public:
        /// Default constructor.
        ALE();

        /// Constructor
        ALE(LibUtilities::SessionReaderSharedPtr &pSession);

/// Destructor
virtual ~ALE();

void SetUpFreeSurfaceSplines();

/// Compute advection term
    void DoALE(const NekDouble aii_Dt, NekDouble time);
    void DeformMeshBoundary(const NekDouble aii_Dt);
    void DeformMesh2(const NekDouble aii_Dt);
    void UpdateFields();
    void UpdateMesh();

    void CreateFreeSurfaceSpline(const NekDouble aii_Dt, NekDouble time);

        void SetMeshVelocityatBoundary(const NekDouble aii_Dt);
        void ComputeNewCoordsatFreeSurfaceBoundary(Array<OneD, Array<OneD, NekDouble> > &newcoords, const NekDouble aii_Dt);
        void UpdateMeshFreeSurfaceBoundary(const NekDouble aii_Dt);
    void ComputeMeshVelocity(const NekDouble aii_Dt);
    void ComputeVelocityMinusMeshVelocity(Array<OneD, Array<OneD, NekDouble> > &velocity);
    void SetFreeSurfaceVelocityBC();
    void WriteArrayAlongFreeSurfaceToFile(string fname, int outputcount,Array<OneD, Array<OneD, NekDouble> > &inarray);

    NekDouble LagrangeInterpolant(NekDouble x, int npts, const Array<OneD, const NekDouble>& xpts,
               const Array<OneD, const NekDouble>& funcvals);
    NekDouble LagrangePoly(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts);
    NekDouble LinearInterpolation(NekDouble x, NekDouble y, int npts, const Array<OneD, const NekDouble>& xfs,
    const Array<OneD, const NekDouble>& yfs, const Array<OneD, const NekDouble>& vfs);


protected:
    int m_phystot;
    Array<OneD, int> m_velocity;
    Array<OneD, int> m_meshvelocity;
    int m_intsteps;
/*  Array<OneD, MultiRegions::ExpListSharedPtr > m_fields;
        SpatialDomains::BoundaryConditionsSharedPtr  m_boundaryConditions;
        int m_spacedim;
       // Array<OneD, MultiRegions::ExpListSharedPtr> m_meshvelocity; */
        Array <OneD, Array<OneD, NekDouble> >m_meshcoords;
        Array <OneD, Array<OneD, NekDouble> >m_meshcoordsold;
        MultiRegions::DisContField2DSharedPtr m_DisField;
        MultiRegions::ExpList2DSharedPtr m_ExpField;
        //Array <OneD, Array<OneD, NekDouble> >m_meshvelocity;

        SpatialDomains::MeshGraph2DSharedPtr         m_mesh2D;
        Array<OneD, Array<OneD, NekDouble> >         m_meshvelwx, m_meshvelwy;
        Array<OneD, LibUtilities::CubicSplineSharedPtr> m_freesurfacesplines;

};

/// Pointer to an AdvectionTerm object.
    typedef boost::shared_ptr<ALE> ALESharedPtr;

} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

