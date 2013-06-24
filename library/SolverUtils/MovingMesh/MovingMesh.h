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


#ifndef NEKTAR_SOLVERS_MOVINGMESH_H
#define NEKTAR_SOLVERS_MOVINGMESH_H

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Foundations/CubicSpline.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <SolverUtils/UnsteadySystem.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <SpatialDomains/MeshGraph2D.h>

namespace Nektar
{
namespace SolverUtils
 {
class MovingMesh;

/// A shared pointer to an EquationSystem object
typedef boost::shared_ptr<MovingMesh> MovingMeshSharedPtr;
/// Datatype of the NekFactory used to instantiate classes derived from
/// the EquationSystem class.
typedef LibUtilities::NekFactory<
		std::string, MovingMesh,
        const LibUtilities::SessionReaderSharedPtr&,
        const SpatialDomains::MeshGraphSharedPtr&,
		const Array<OneD, MultiRegions::ExpListSharedPtr>&,
		const SpatialDomains::BoundaryConditionsSharedPtr&
    > MovingMeshFactory;
MovingMeshFactory& GetMovingMeshFactory();

static NekDouble kFreeSurfaceBCsExtrapolation[][3] = {{ 1.0,  0.0, 0.0},
                                                    { 1.5, -0.5, 0.0},
                                                    { 23.0/12.0, -16.0/12.0, 5.0/12.0}};

static NekDouble kFSBCsExtrapolation[][3] = {{ 1.0,  0.0, 0.0},
                                                    { 2.0, -1.0, 0.0},
                                                    { 3.0, -3.0, 1.0}};
    class MovingMesh : public UnsteadySystem
    {
    	public:
        /// Default constructor.
        /// Constructor
        MovingMesh(const LibUtilities::SessionReaderSharedPtr &pSession);

        virtual void v_InitObject();

        static MovingMeshSharedPtr create(const LibUtilities::SessionReaderSharedPtr& pSession,
        		const SpatialDomains::MeshGraphSharedPtr   &pGraph,
        		const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        		const SpatialDomains::BoundaryConditionsSharedPtr  &pBoundaryConditions) {
            MovingMeshSharedPtr p = MemoryManager<MovingMesh>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

		/// Destructor
		virtual ~MovingMesh();

		inline void InitObject();

		void SetUpFreeSurfaceSplines();

		/// Compute advection term
		void DoMovingMeshInitialise(const NekDouble timestep);
	    void DoMovingMesh(NekDouble timestep,NekDouble time);
	    void DoALEAxiSymmetric(const NekDouble timestep,
	    			const NekDouble time);

	    void SetMeshVelocityatBoundary(const NekDouble timestep);
	    void DeformMeshBoundary(const NekDouble timestep);
	    void DeformMesh(const NekDouble timestep);
	    void DeformMesh2(const NekDouble timestep);
	    void UpdateFields();
	    void UpdateMesh();
	    void UpdateMeshFreeSurfaceBoundary(const NekDouble timestep);

	    void ComputeNewCoordsatFreeSurfaceBoundary(Array<OneD, Array<OneD, NekDouble> > &newcoords, const NekDouble timestep);

	    void CreateFreeSurfaceFunction(const NekDouble timestep);
	    void CreateFreeSurfaceSpline(const NekDouble timestep,
				const NekDouble time);

	    void ComputeMeshVelocity(const NekDouble timestep);
	    void ComputeVelocityMinusMeshVelocity(
				Array<OneD, Array<OneD, NekDouble> > &velocity);
	    void AddFreeSurfaceTension(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray, const NekDouble tension);
	    void AddFreeSurfaceTension2(
	    		const NekDouble tension,
	    		MultiRegions::ExpListSharedPtr pPressure,
	    		Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);
	    void SetFreeSurfaceVelocityBC();
	    void SetFreeSurfaceVelocityBC2(MultiRegions::ExpListSharedPtr pPressure, NekDouble tension);
	    void SetDirichletUnzeroBCALE();
	    void SmoothNormalsAlongFreeSurface(const NekDouble timestep);
	    void WriteValuesAlongFreeSurfaceToFile(string fname, int outputcount);
	    void WriteStressAlongFreeSurfaceToFile(string fname, int outputcount,MultiRegions::ExpListSharedPtr pPressure);
	    void WriteArrayAlongFreeSurfaceToFile(string fname, int outputcount,Array<OneD, Array<OneD, NekDouble> > &inarray);
	    void SetInitialMeshSloshing();
	    void SetInitialMeshDieSwell();
	    void Subdivwtu(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
	        	                                                 Array<OneD, Array<OneD, NekDouble> > &outarray);

	    void AddwnutoFreeSurface( Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);
	    void SetEdgePhysValsToZero(int edge, Array<OneD, NekDouble> &inarray, StdRegions::StdExpansionSharedPtr &ElExp);
	    void SetValueAlongFreeSurfaceToZero(Array<OneD, Array<OneD, NekDouble> > &inarray);
	    void SetOutflowBCALE();
	    void SetPlateBCALE();
	    void WriteSwellRatioToFileFromSpline();

	protected:
	    int 										m_phystot;
	    Array<OneD, int> 							m_velocity;
	    Array<OneD, int> 							m_meshvelocity;
	    int 										m_intsteps;
	/*    NekDouble m_reynolds;
	    NekDouble m_weber;


        Array<OneD, MultiRegions::ExpListSharedPtr > m_fields;
        SpatialDomains::BoundaryConditionsSharedPtr  m_boundaryConditions;

        int m_spacedim;
       // Array<OneD, MultiRegions::ExpListSharedPtr> m_meshvelocity; */
	    LibUtilities::SessionReaderSharedPtr        m_session;
	    int m_spacedim;

        Array <OneD, Array<OneD, NekDouble> >		m_meshcoords;
        Array <OneD, Array<OneD, NekDouble> >		m_meshcoordsold;
        MultiRegions::DisContField2DSharedPtr 		m_DisField;
        MultiRegions::ExpList2DSharedPtr 			m_ExpField;
        //Array <OneD, Array<OneD, NekDouble> >		m_meshvelocity;

        SpatialDomains::MeshGraph2DSharedPtr         m_mesh2D;
        Array<OneD, Array<OneD, NekDouble> >         m_meshvelwx, m_meshvelwy;
        Array<OneD, Array<OneD, NekDouble> >         m_FBCwx, m_FBCwy;
        Array<OneD, Array<OneD, NekDouble> >         m_FBCu;
        Array<OneD, LibUtilities::CubicSplineSharedPtr> m_freesurfacesplines;

        //Array<OneD, bool>                               m_checkIfSystemSingular;
        //MultiRegions::GlobalSysSolnType 				m_solnType;

	};

    inline void MovingMesh::InitObject()
    {
        v_InitObject();
    }
 }
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H


