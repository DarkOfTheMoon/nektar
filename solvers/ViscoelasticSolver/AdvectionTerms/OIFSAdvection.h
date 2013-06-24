///////////////////////////////////////////////////////////////////////////////
//
// File LinearisedAdvection.h
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
// Description: TBA
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_OIFSAdvection_H
#define NEKTAR_SOLVERS_OIFSAdvection_H

#include <ViscoelasticSolver/AdvectionTerms/AdvectionTerm.h>
//#include <Auxiliary/ADRBase.h>

//#define TIMING
//#ifdef TIMING
//#include <time.h>
//#include <sys/time.h>
//#endif


namespace Nektar
{


    class OIFSAdvection: public AdvectionTerm

    {
    public:
    	friend class MemoryManager<OIFSAdvection>;

		/// Creates an instance of this class
		static AdvectionTermSharedPtr create(
								const LibUtilities::SessionReaderSharedPtr& pSession,
								const SpatialDomains::MeshGraphSharedPtr& pGraph) {
			AdvectionTermSharedPtr p = MemoryManager<OIFSAdvection>::AllocateSharedPtr(pSession, pGraph);
			p->InitObject();
			return p;
		}
		/// Name of class
		static std::string className;

	protected:

        OIFSAdvection(
                const LibUtilities::SessionReaderSharedPtr&        pSession,
                const SpatialDomains::MeshGraphSharedPtr&          pGraph);


        virtual ~OIFSAdvection();

		//Virtual function for the evaluation of the advective terms
		virtual void v_DoAdvection(
								   Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
								   const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
								   Array<OneD, Array<OneD, NekDouble> > &pOutarray,
								   Array<OneD, NekDouble> &pWk);

		void DoInitialise(Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
				LibUtilities::TimeIntegrationMethod IntMethod,
		            								 const Array<OneD, const Array<OneD, NekDouble> > &y_0,
		            								 const NekDouble &timestep,
		            								 const NekDouble &time);

        void DoAdvection(Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        										LibUtilities::TimeIntegrationSolutionSharedPtr &solvector,
        										const NekDouble timestep, const NekDouble time);

	protected:

		int   m_nConvectiveFields;  /// Number of fields to be convected;
		int   m_numMultiSteps;
		int   m_ncalls;
		int   m_nsteps;
		int   m_RK4itmax;

        Array<OneD, int> m_velocityid; ///< int which identifies which components of m_fields contains the velocity (u,v,w);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_velocity;
        Array<OneD, NekDouble> m_time;
        Array<OneD, NekDouble> m_timestep;
        NekDouble m_dt;

        /// The time integration method to use.
        LibUtilities::TimeIntegrationMethod             m_timeIntMethod;
        /// The time integration scheme operators to use.
        LibUtilities::TimeIntegrationSchemeOperators    m_ode;

        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> m_intScheme;
        Array<OneD, LibUtilities::TimeIntegrationSolutionSharedPtr> m_utilde;

        /// Array holding velocity field
        Array<OneD, MultiRegions::ExpListSharedPtr>     m_fields;



        //Function for the evaluation of the linearised advective terms
        void ComputeAdvectionTerm(
                Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                const Array<OneD, Array<OneD, NekDouble> > &pV,
                const Array<OneD, const NekDouble> &pU,
                Array<OneD, NekDouble> &pOutarray,
                int pVelocityComponent,
                Array<OneD, NekDouble> &pWk);

        void DoOdeRhs(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                      Array<OneD,  Array<OneD, NekDouble> > &outarray,
                      const NekDouble time);

        void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                          Array<OneD,  Array<OneD, NekDouble> > &outarray,
                          const NekDouble time);

        void WeakAdvectionGreensDivergenceForm(
                          const Array<OneD, Array<OneD, NekDouble> > &F,
                          Array<OneD, NekDouble> &outarray);

        void WeakDGAdvection(
                          const Array<OneD, Array<OneD, NekDouble> >& InField,
                          Array<OneD, Array<OneD, NekDouble> >& OutField,
                          bool NumericalFluxIncludesNormal,
                          bool InFieldIsInPhysSpace);

        // DG Advection routines
        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux);

	};

    /// A shared pointer to an EquationSystem object
	typedef boost::shared_ptr<OIFSAdvection> OIFSAdvectionSharedPtr;


} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H

