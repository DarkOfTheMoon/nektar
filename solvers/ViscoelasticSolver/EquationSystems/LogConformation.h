///////////////////////////////////////////////////////////////////////////////
//
// File LogConformation.h
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
// Description: Viscoelastic constitutive equations in log conformation formulation
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERS_VISCOELASTIC_LOGCONFORMATION_H_
#define NEKTAR_SOLVERS_VISCOELASTIC_LOGCONFORMATION_H_

//#define OUTPUT

#include <ViscoelasticSolver/EquationSystems/Viscoelastic.h>

namespace Nektar
{
    class LogConformation: public Viscoelastic
    {
    public:

        /**
         * Default constructor.
         *
         */
    	LogConformation();


        /**
         * Constructor.
         * \param
         *
         *
         */
    	LogConformation(const LibUtilities::SessionReaderSharedPtr &pSession);

    	virtual void v_InitObject();

	    void SetInitialConditionsLogc(NekDouble initialtime);

    	//void AddStressDivergence(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    	                     //                                Array<OneD, Array<OneD, NekDouble> > &outarray);

    protected:

        // int   m_nConvectiveFields;  /// Number of fields to be convected;

        /// Array holding all dependent variables.
    	Array<OneD, int> m_logc; ///< int which identifies which components of m_fields contain the stress tensor components t_xx, t_xy, t_yy

  /*  private:
        EquationType  m_equationType;  ///< equation type;
        AdvectionForm m_advectionForm; ///< Form of advection terms. */


	//void SetBoundaryConditions(NekDouble time);
    	/*!
    	 Computes the analytical solution of the
    	 spectral decomposition of an symmetric positive definite tensor T in 2D.
    	 @param[in] phystot number of physical points for which the spectral decomposition is computed
    	 @param[in] T in terms of the tensor components T[0]=Txx, T[1]=Txy, T[2]=Tyy
    	 @param[out]  EigVal  give back the Eigenvalues EigValue[0]=lambda1 and EigValue[1]=lambda2
    	 @param[out]  EigVec   contains the information about the eigenvector tensor R in terms of
							   EigVec[0]=e1^2, EigVec[1]=e1*e2, EigVec[2]=e2^2
    	 */

    	  void CompSpectralDecomp2D(int phystot, const Array<OneD, Array<OneD, NekDouble> > &T, Array<OneD, Array<OneD, NekDouble> > &EigVal,
                               Array<OneD, Array<OneD, NekDouble> > &EigVec);

      	/*!
      	 Computes the transformation of a 2D symmetric positive definite tensor from its principal frame
      	 into the global frame:
      	 T = R*Tp*R^t (R: matrix of eigenvectors, Tp: tensor T in prin)
      	 @param[in] phystot number of physical points
      	 @param[in] Tp tensor in its principal frame in terms of the tensor components Tp[0]=Tpxx, Tp[1]=Tpxy, Tp[2]=Tpyy
      	 @param[in] EigVec Eigenvector tensor R in terms of EigVec[0]=e1^2, EigVec[1]=e1*e2, EigVec[2]=e2^2
      	 @param[out]  T  tensor components T[0]=Txx, T[1]=Txy, T[2]=Tyy in the global frame
      	 */

      	  void PrincipalToGlobalTrans(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Tp, const Array<OneD, Array<OneD, NekDouble> > &EigVec,
                                 Array<OneD, Array<OneD, NekDouble> > &T);

		/*!
		 Computes the transformation of a 2D symmetric positive definite tensor from its principal frame
		 into the global frame:
		 T = R*Tp*R^t (R: matrix of eigenvectors, Tp: tensor T in prin)
		 @param[in] phystot number of physical points
		 @param[in] T tensor components T[0]=Txx, T[1]=Txy, T[2]=Tyy in the global frame
		 @param[in] EigVec Eigenvector tensor R in terms of EigVec[0]=e1^2, EigVec[1]=e1*e2, EigVec[2]=e2^2
		 @param[out]  Tp  tensor components Tp[0]=Tpxx, Tp[1]=Tpxy, Tp[2]=Tpyy in its principal frame
		 */

		  void GlobalToPrincipalTrans(int phystot, const Array<OneD, Array<OneD, NekDouble> > &T, const Array<OneD, Array<OneD, NekDouble> > &EigVec,
							   Array<OneD, Array<OneD, NekDouble> > &Tp);

		  void ElasticToConformationTensor(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Elastic,
										Array<OneD, Array<OneD, NekDouble> > &C);
		  void ElasticToLogConformationTensor(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Elastic,
				  Array<OneD, Array<OneD, NekDouble> > &Logc);
		  void LogConformationTensorToElastic(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Logc,
				  Array<OneD, Array<OneD, NekDouble> > &Elastic);
		  void LogConformationToConformationTensor(int phystot, const Array<OneD, Array<OneD, NekDouble> > &Logc,
				  Array<OneD, Array<OneD, NekDouble> > &C);
		  void SetBoundaryConditionsLogc(NekDouble time);
		  void SetBoundaryInflowLogcDG(int bcRegion, int cnt, NekDouble time);
		  void SetBoundaryInflowLogc(int bcRegion, int cnt, NekDouble time);

		  void EvaluateRhsViscoelasticEquationLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		  			    	                Array<OneD, Array<OneD, NekDouble> > &outarray);
		  void CalcWeakDGSourceTermLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		  	            Array<OneD, Array<OneD, NekDouble> > &Weakoutarray, Array<OneD, Array<OneD, NekDouble> > &outarray);
		  void CalcSourceTermLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
		  		                Array<OneD, Array<OneD, NekDouble> > &outarray);
		  void LogcToElasticL2Projection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
				  Array<OneD, Array<OneD, NekDouble> > &outarray);
	      void AddStressDivergenceLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
	                                                     Array<OneD, Array<OneD, NekDouble> > &outarray);
	      void AddStressDivergenceWeakFormLogc(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
	      													  Array<OneD, Array<OneD, NekDouble> > &forcing);
	      void CalcElasticToLogConformationTensorforCriticalPoints(int id, Array<OneD, Array<OneD, NekDouble> > &Logc);

        // Virtual functions
        /* virtual void v_AdvanceInTime(int nsteps)
        {
            ASSERTL0(false,"This method is not defined in this class");
        } */

    };

    typedef boost::shared_ptr<LogConformation> LogConformationSharedPtr;

} //end of namespace

#endif /* NEKTAR_SOLVERS_LOGCONFORMATION_H_ */
