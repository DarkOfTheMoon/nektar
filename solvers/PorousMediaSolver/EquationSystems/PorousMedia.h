///////////////////////////////////////////////////////////////////////////////
//
// File PorousMedia.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_POROUSMEDIA_H
#define NEKTAR_SOLVERS_POROUSMEDIA_H

#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/AdvectionSystem.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <PorousMediaSolver/EquationSystems/DarcyTerm.h>
#include <PorousMediaSolver/EquationSystems/PressureTerm/Extrapolate.h>
#include <SolverUtils/Forcing/Forcing.h>

namespace Nektar
{     

    enum EquationType
    {
        eNoEquationType,
        eUnsteadyPorousMedia,
        eEquationTypeSize
    };
    
    // Keep this consistent with the enums in EquationType.
    const std::string kEquationTypeStr[] = 
    {
        "NoType",
        "UnsteadyPorousMedia",
    };

    enum DarcyTermMethod
    {
        eExplicit,
        eExplicitSpatiallyVarying,
        eImplicitIsotropic,
        eImplicitAnisotropic,
        SIZE_DarcyTerm        //!< Length of enum list
    };

    // Keep this consistent with the enums in EquationType.
    const std::string DarcyTermMethodStr[] = 
    {
        "Explicit",
        "ExplicitSpatiallyVarying",
        "ImplicitIsotropic",
        "ImplicitAnisotropic"
    };


    enum AdvectionForm
    {
        eNoAdvectionForm,
        eConvective,
        eNonConservative,
        eLinearised,
        eAdjoint,
        eSkewSymmetric,
        eNoAdvection,
        eAdvectionFormSize
    };
    
    // Keep this consistent with the enums in EquationType.
    const std::string kAdvectionFormStr[] = 
    {
        "NoType",
        "Convective",
        "NonConservative",
        "Linearised",
        "Adjoint",
        "SkewSymmetric"
        "NoAdvection"
        };
	
    /**
     * \brief This class is the base class for Navier Stokes problems
     *
     */
    
    //    class PorousMedia: public SolverUtils::EquationSystem
    class PorousMedia: public SolverUtils::AdvectionSystem
    {
    public:           
        // Destructor
        virtual ~PorousMedia();

        virtual void v_InitObject();

        void AddForcing(const SolverUtils::ForcingSharedPtr& pForce);
        
    protected: 

        /**
         * Constructor.
         */
        PorousMedia(const LibUtilities::SessionReaderSharedPtr& pSession);

        LibUtilities::TimeIntegrationSolutionSharedPtr  m_integrationSoln;

        /// Advection term
        //AdvectionTermSharedPtr m_advObject;

        /// Forcing terms
        std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

        ExtrapolateSharedPtr m_extrapolation;

        /// Number of fields to be convected; 
        int   m_nConvectiveFields;  

        /// int which identifies which components of m_fields contains the velocity (u,v,w);
        Array<OneD, int> m_velocity; 
 

        /// Pointer to field holding pressure field
        MultiRegions::ExpListSharedPtr m_pressure;  
        
        NekDouble   m_kinvis;        ///< Kinematic viscosity

        Array<OneD, NekDouble>  m_perm;  ///< Permeability matrix
        Array<OneD, Array<OneD, NekDouble> >  m_spatialperm;  ///< Permeability matrix
        Array<OneD, NekDouble>  m_perm_inv;  ///< inverse Permeability matrix

        Array<OneD, NekDouble>  m_darcy_fac;  ///< inverse Permeability matrix

        /// Variable permeability
        StdRegions::VarCoeffMap m_varperm;

        int         m_steadyStateSteps; ///< Check for steady state at step interval
        NekDouble   m_steadyStateTol; ///< Tolerance to which steady state should be evaluated at

        EquationType  m_equationType;  ///< equation type;
        DarcyTermMethod m_darcyType;  ///< equation type;

        int m_intSteps;  ///< Number of time integration steps AND  Order of extrapolation for pressure boundary conditions.         

        EquationType GetEquationType(void)
        {
            return m_equationType;
        }

        virtual int v_GetForceDimension()=0;

        void EvaluateAdvectionTerms(const Array<OneD, 
                                    const Array<OneD, NekDouble> > &inarray, 
                                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                    Array<OneD, NekDouble> &wk = NullNekDouble1DArray);
		
        //time dependent boundary conditions updating
        void SetBoundaryConditions(NekDouble time);

        virtual void v_PrintSummary(std::ostream &out)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

        virtual void v_DoInitialise(void)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }
	
	DarcyTermSharedPtr m_darcyEvaluation;
	
    private: 
    };
    
    typedef boost::shared_ptr<PorousMedia> PorousMediaSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_POROUSMEDIA_H
