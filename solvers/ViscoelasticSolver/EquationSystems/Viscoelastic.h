///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokes.h
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
// Description: Viscoelastic constitutive equations
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VISCOELASTIC_VISCOELASTIC_H
#define NEKTAR_SOLVERS_VISCOELASTIC_VISCOELASTIC_H

#include<ViscoelasticSolver/EquationSystems/IncNavierStokes.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
	static NekDouble kStressDivergenceExtrapolation[][3] = {{ 1.0,  0.0, 0.0},
                                                    { 2.0, -1.0, 0.0},
                                                    { 3.0, -3.0, 1.0}};
	enum StressTreatment
	{
		eNoStressTreatment,
		eExplicit,
		eImplicit,
		eStressTreatmentSize
	};

	// Keep this consistent with the enums in EquationType.
	const std::string kStressTreatmentStr[] =
	{
		"NoType",
		"Explicit",
		"Implicit"
	};

    class Viscoelastic: public IncNavierStokes
    {
    public:

        /**
         * Constructor.
         * \param
         *
         *
         */
    	Viscoelastic(const LibUtilities::SessionReaderSharedPtr &pSession);

    	virtual void v_InitObject();

    	void SetInitialConditionsStress(NekDouble initialtime);

    	void AddStressDivergence(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    	                                                     Array<OneD, Array<OneD, NekDouble> > &outarray);
    	void SubGGTDivergence(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    	                                                          Array<OneD, Array<OneD, NekDouble> > &outarray);
    	void SubGDivergence(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    	                                                          Array<OneD, Array<OneD, NekDouble> > &outarray);
        void AddStressDivergenceWeakFormExplicit(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                                      Array<OneD, Array<OneD, NekDouble> > &forcing);
        void CalcStressDivergenceWeakForm(const Array<OneD, const Array<OneD, NekDouble> > &stressin,
        												  Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);
        void AddStressTimesNormalToVelocityNeumannBC(const Array<OneD, const Array<OneD, NekDouble> > &stressin,
        											  Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);

        void AddStressDivergenceWeakFormExplicit2(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                                      Array<OneD, Array<OneD, NekDouble> > &forcing);

        void AddStressDivergenceWeakFormExplicit3(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                                      Array<OneD, Array<OneD, NekDouble> > &forcing);

    //	void SetBoundaryConditions(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
      //          Array<OneD, Array<OneD, NekDouble> > &outarray);

    	void EvaluateRhsViscoelasticEquation(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    	                Array<OneD, Array<OneD, NekDouble> > &outarray);

    	void EvaluateViscoelasticEquationImplicitely(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    	    		    	                Array<OneD, Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt);

    	void CalcSourceTerm(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &outarray);

    	void CalcSourceTermImplicit(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    			                Array<OneD, Array<OneD, NekDouble> > &outarray);

    	void StressAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &outarray, int vstart);

    	void CalcWeakDGSourceTerm(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);

    	void StressWeakDGAdvection(const Array<OneD, const Array<OneD, NekDouble> > &stressinarray,
    			                Array<OneD, Array<OneD, NekDouble> > &Weakoutarray,int vstart);

    	void WeakDGAdvection( Array<OneD, Array<OneD, NekDouble> > &velocity,
    	                      const Array<OneD, Array<OneD, NekDouble> >& InField,
    	                      Array<OneD, Array<OneD, NekDouble> >& OutField,
    	                      bool NumericalFluxIncludesNormal,
    	                      bool InFieldIsInPhysSpace,
    	                      int nvariables,int vstart);

        // DG Advection routines
        void GetFluxVector(Array<OneD, Array<OneD, NekDouble> > &velocity, const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &flux, int vstart);

        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &velocity,
        		Array<OneD, Array<OneD, NekDouble> > &physfield,
        		Array<OneD, Array<OneD, NekDouble> > &numflux,int vstart);

        void WeakAdvectionGreensDivergenceForm(
                        const Array<OneD, Array<OneD, NekDouble> > &F,
                        Array<OneD, NekDouble> &outarray, int vstart);

        void WeakPenaltyforScalar(int var,
                     	 Array<OneD, NekDouble> &physfield,
                           Array<OneD, NekDouble> &penaltyflux,
                           NekDouble time);

    	void SetBoundaryConditionsStress(const Array<OneD, const Array<OneD, NekDouble> > &physarray, NekDouble time);

    	void SetStressWallBoundaryConditions(const Array<OneD, const Array<OneD, NekDouble> > &inarray,const NekDouble aii_Dt);
    	void SetInflowBoundaryConditionsGiesekus(const Array<OneD, const Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt);
    	void Determinante(const NekDouble phystot, const Array<OneD, const Array<OneD, NekDouble> > &M, Array<OneD, NekDouble> &Det);
    	void StressMatrixInverseOldB(const NekDouble phystot, const Array<OneD, const Array<OneD, NekDouble> > &M, Array<OneD, Array<OneD, NekDouble> > &Mout);
    	void StressMatrixInversetimesRHS(const NekDouble phystot, Array<OneD, Array<OneD, NekDouble> > &M,
    			 Array<OneD, Array<OneD,NekDouble> > &RHS, Array<OneD, Array<OneD, NekDouble> > &tau);
    	void SetUpMatrixandRHSOldB(const NekDouble phystot, const Array<OneD, const Array<OneD, NekDouble> > &inarray, Array<OneD, Array<OneD, NekDouble> > &D,
    			                Array<OneD, Array<OneD, NekDouble> > &M,Array<OneD, Array<OneD, NekDouble> > &RHS,const NekDouble aii_Dt);
    	void VelocityGradientTensor(Array<OneD, Array<OneD, NekDouble> > &outarray,Array<OneD, Array<OneD, NekDouble> > &L);
    	void AddGiesekusTermtoRHS(const NekDouble phystot, const Array<OneD, const Array<OneD, NekDouble> > &stressin,
    				Array<OneD, Array<OneD, NekDouble> > &RHS,Array<OneD, Array<OneD, NekDouble> > &RHSpN);
    	void CalcSourceTermImplicitOldB(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    		                Array<OneD, Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt);
    	void CalcSourceTermImplicitGiesekus(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    			                Array<OneD, Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt);

    	void SetStressWallBoundaryConditionsOldB(const Array<OneD, const Array<OneD, NekDouble> > &inarray,const NekDouble aii_Dt);

    	void SetStressSymmetryBoundaryConditionsOldB(const Array<OneD, const Array<OneD, NekDouble> > &inarray,const Array<OneD, const Array<OneD, NekDouble> > &outarray,const NekDouble aii_Dt);

		/* functions to compute waters analytical solution for time dependent poiseuille flow for an Oldroyd-B fluid */
		void SetBoundaryWatersTAU(int bcRegion, int cnt, NekDouble time);

		void SetBoundaryWatersU2(int bcRegion, int cnt, NekDouble time);

		NekDouble CalcDirichletWatersU(NekDouble y, NekDouble time);

		void CalcDirichletWatersTAU(NekDouble y, NekDouble &tauxx, NekDouble &tauxy, NekDouble time);

		void CalcDirichletWatersTAUCompConstants(NekDouble y, NekDouble &Cxx, NekDouble &Cxy);

		void SetDirichletWatersSteadyState(NekDouble y, NekDouble &tauxx, NekDouble &tauxy);

		void WriteWaterstoFile(int nchk,NekDouble x, NekDouble y, NekDouble xtol, NekDouble ytol);

		void WriteMinMaxtoFile(const Array<OneD, const Array<OneD, NekDouble> > &physfield, int vstart, int nvariables, string fname, bool IsFromMFields);

		void WriteMinMaxtoFileHeader(int vstart, int nvariables, string fname);

		void CalcDrag(const Array<OneD, const Array<OneD, NekDouble> > &inarray, NekDouble &sumdrag, NekDouble &sumlift);
		void WriteDragToFile(int count,const Array<OneD, const Array<OneD, NekDouble> > &inarray);
		void WriteDragtoFileHeader(string fname);
		void PrintStressAlongBoundary(const Array<OneD, const Array<OneD, NekDouble> > &inarray, string fname, int count);
		void SetOutarrayToZeroAlongWall(Array<OneD, Array<OneD, NekDouble> > &inarray);
		void PrintValueAlongBoundary(const Array<OneD, const Array<OneD, NekDouble> > &inarray, string fname, int count);
		void PrintValueAtCylinder(const Array<OneD, const Array<OneD, NekDouble> > &inarray, string fname);
		void PrintValueAtCylinderTrace( Array<OneD, NekDouble>  &inarray, string fname);
		void SetOutarrayToZeroAlongCylinder(Array<OneD, Array<OneD, NekDouble> > &inarray);
		void SolveOlroydBPseudoSpectral(Array<OneD, Array<OneD, NekDouble> > &fieldsit,
		        Array<OneD, Array<OneD, NekDouble> > &fieldsn,
		        const NekDouble aii_Dt);

		/* Functions for DEVSS scheme and Computation of Velocity Gradient Tensor G
		 *
		 */
        void AddGGTDivergenceWeakForm(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);
        void AddGGTTimesNormalToVelocityNeumannBC(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);
        void CalcGL2Projection(Array<OneD, Array<OneD, NekDouble> > &outarray);
        void AddGDivergenceWeakForm(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);
        void AddGTimesNormalToVelocityNeumannBC(Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);

        void UpdateG();

        void SetFixedCylinderVelocityField();

	    inline NekDouble GetGNexp(NekDouble t,NekDouble fac, NekDouble a, NekDouble p , NekDouble b, NekDouble q)
		{
			return (fac*( (a*exp(p*t)) + (b*exp(q*t))) );
		}

		inline NekDouble GetGNcossin(NekDouble t,NekDouble fac, NekDouble a, NekDouble p , NekDouble b, NekDouble q)
		{
			return (fac*(a*cos(p*t) + b*sin(q*t)));
		}

		inline NekDouble GetGN0exp(NekDouble fac, NekDouble a, NekDouble b)
		{
			return (fac*(a+b) );
		}

		inline NekDouble GetGN0cossin(NekDouble fac, NekDouble a)
		{
			return (fac*a);
		}

    protected:

        // int   m_nConvectiveFields;  /// Number of fields to be convected;

        /// Array holding all dependent variables.
    	Array<OneD, int> m_stress; ///< int which identifies which components of m_fields contain the stress tensor components t_xx, t_xy, t_yy
    	Array<OneD, int> m_Gindex; ///< int which identifies which components of m_fields contain the stress tensor components G11, G12, G21, G22
    	/// Array G11, G12, G21, G22.
    	Array<OneD, MultiRegions::ExpListSharedPtr> m_G;
    	Array<OneD, Array<OneD, NekDouble> > m_stressold;
    	Array<OneD, Array<OneD, NekDouble> >  m_stressdivergenceweaku; //< Storage for current and previous levels of stress divergence in weak form
    	Array<OneD, Array<OneD, NekDouble> >  m_stressdivergenceweakv; //< Storage for current and previous levels of stress divergence in weak form

        NekDouble     m_relaxationtime;                ///< Relaxation time lambda
        NekDouble     m_weissenberg;                ///< Weissenberg number
        NekDouble     m_alpha;                ///< anisotropy parameter for Giesekus model
        NekDouble     m_xi;                ///< shear thinning parameter PTT model
        NekDouble     m_eps;                ///< elongational parameter PTT model
        NekDouble     m_lptt;                ///< calculate linear PTT model instead of standard expptt
        NekDouble     m_gamma;
        NekDouble     m_theta;
        // NekDouble     m_pvis; 						///< Polymeric Viscosity eta_p

        Array<OneD, MultiRegions::ExpListSharedPtr> Derived_field;
        StressTreatment m_stressTreatment;

    };

    typedef boost::shared_ptr<Viscoelastic> ViscoelasticSharedPtr;

} //end of namespace



#endif /* VISCOELASTIC_H_ */
