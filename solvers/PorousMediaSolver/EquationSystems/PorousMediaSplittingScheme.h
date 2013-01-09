///////////////////////////////////////////////////////////////////////////////
//
// File BrinkmanSplittingScheme.h
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
// Description: PorousMedia's equation routine 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_BRINKMANSPLITTINGSCHEME_H
#define NEKTAR_SOLVERS_BRINKMANSPLITTINGSCHEME_H

#include <PorousMediaSolver/EquationSystems/PorousMedia.h>

namespace Nektar
{     
    
    static NekDouble kHighOrderBCsExtrapolation[][3] = {{ 1.0,  0.0, 0.0},
    { 2.0, -1.0, 0.0},
    { 3.0, -3.0, 1.0}};
    
    /**
     * \brief This class is the base class for the Velocity Correction Scheme
     *
     */
    
    class PorousMediaSplittingScheme: public PorousMedia
    {
    public:           
        
        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession) {
            SolverUtils::EquationSystemSharedPtr p = MemoryManager<PorousMediaSplittingScheme>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;
        
            
        /**
         * Constructor.
         * \param 
         * 
         */
        PorousMediaSplittingScheme(const LibUtilities::SessionReaderSharedPtr& pSession);
        
        virtual ~PorousMediaSplittingScheme();
            
        virtual void v_InitObject();
        
        void EvaluatePressureBCs(
                const Array<OneD, const Array< OneD,  NekDouble> > &fields, 
                const Array<OneD, const Array< OneD,  NekDouble> > &N, 
                const NekDouble Aii_Dt = NekConstants::kNekUnsetDouble);
        
        void SetUpPressureForcing(
            const Array<OneD, const Array<OneD, NekDouble> > &fields, 
                Array<OneD, Array<OneD, NekDouble> > &Forcing, 
            const NekDouble aii_Dt);
        
        void SetUpViscousForcing(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            Array<OneD, Array<OneD, NekDouble> > &Forcing, 
            const NekDouble aii_Dt);
        
        void SolveUnsteadyStokesSystem(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            Array<OneD, Array<OneD, NekDouble> > &outarray, 
            const NekDouble time,
            const NekDouble a_iixDt);
        
        void EvaluateAdvection_SetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            Array<OneD, Array<OneD, NekDouble> > &outarray, 
            const NekDouble time);
        
    protected:
        
    private: 
        int m_pressureCalls;
        
        Array<OneD, int> m_pressureBCtoElmtID;  // Id of element to which pressure  boundary condition belongs
        Array<OneD, int> m_pressureBCtoTraceID; // Id of edge (2D) or face (3D) to which pressure boundary condition belongs
        
        Array<OneD, Array<OneD, NekDouble> >  m_pressureHBCs; //< Storage for current and previous levels of high order pressure boundary conditions. 
        
        Array<OneD, Array<OneD, int> > m_HBC;  //data structure to old all the information regarding High order pressure BCs
        
        int m_HBCnumber;                       // number of elemental expansion where a boundary is of High order type
        
        StdRegions::StdExpansionSharedPtr m_elmt; // general standard element used to deaal with HOPBC calculations
        
        Array<OneD, NekDouble> m_wavenumber;
        
        Array<OneD, NekDouble> m_beta;

        /// Variable diffusivity
        StdRegions::VarCoeffMap m_vardiff;
        
        /**  \brief This function evaluates the normal Neumann pressure boundary
         *  condition for the velocity correction scheme at the current time
         *  level which requires as input the non-linear terms at the current
         *  time level
         *
         *   \f[ \frac{\partial p}{\partial n}^n = \left [ {\bf N( u)}^n -
         *      kinvis \nabla \times \nabla {\bf u}^n \right] \cdot {\bf n} \f]
         *                                                                           
         * where \f$ {\bf n}\f$ is the unit outward normal along the edge, \f$
         * {\bf u} \f$ is the velocity field, and \f$ {\bf N(u)}\f$ are the
         * non-linear terms in the momentum equation.
         */
        
        void CalcPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields, 
            const Array<OneD, const Array<OneD, NekDouble> >  &N);
        
        void CalcPressureBCs2D(
            const Array<OneD, const Array<OneD, NekDouble> > &fields, 
            const Array<OneD, const Array<OneD, NekDouble> >  &N);
        
        void CalcPressureBCs3D(
            const Array<OneD, const Array<OneD, NekDouble> > &fields, 
            const Array<OneD, const Array<OneD, NekDouble> >  &N);
        
        // Virtual functions 
        virtual void v_PrintSummary(std::ostream &out);
        
        virtual void v_DoSolve(void);
        
        virtual void v_TransCoeffToPhys(void);
        
        virtual void v_TransPhysToCoeff(void);
        
        virtual void v_DoInitialise(void);
        
        virtual Array<OneD, bool> v_GetSystemSingularChecks();
    };
    
    
    typedef boost::shared_ptr<PorousMediaSplittingScheme> PorousMediaSplittingSchemeSharedPtr;
    
} //end of namespace

#endif //BRINKMANSPLITTINGSCHEME_H
