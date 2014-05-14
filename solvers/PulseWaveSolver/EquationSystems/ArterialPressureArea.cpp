///////////////////////////////////////////////////////////////////////////////
//
// File ArterialPressureArea.cpp
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
// Description: ROuflow class
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/ArterialPressureArea.h>
#define PI 3.14159265359
namespace Nektar
{

    std::string ArterialPressureArea::className
    = GetPressureAreaFactory().RegisterCreatorFunction(
        "Arterial",
        ArterialPressureArea::create,
        "Pressure area relationship for the arterial system");

    /**
     *
     */
    ArterialPressureArea::ArterialPressureArea(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel, 
                                               const LibUtilities::SessionReaderSharedPtr pSession,
                                               const int nDomains)
        : PulseWavePressureArea(pVessel,pSession,nDomains)
    {
        m_session->LoadParameter("rho", m_rho, 0.5);

        m_beta = Array<OneD, Array<OneD, NekDouble> >(nDomains);
        m_beta_trace = Array<OneD, Array<OneD, NekDouble> >(nDomains);
    }

    /**
     *
     */
    ArterialPressureArea::~ArterialPressureArea()
    {

    }

    void ArterialPressureArea::v_ReadParameters(int omega, int nqTrace, MultiRegions::ExpListSharedPtr &field)
    {
        int nq=m_vessels[2*omega]->GetNpoints();
        m_beta[omega] = Array<OneD, NekDouble>(nq);
        EvaluateFunction(m_vessels,m_session,"beta", m_beta[omega],"MaterialProperties",0.0,omega,nq);

        m_beta_trace[omega] = Array<OneD, NekDouble>(nqTrace);
        field->ExtractTracePhys(m_beta[omega],m_beta_trace[omega]);
    }

    void ArterialPressureArea::v_GetPacons(int omega, Array<OneD, Array<OneD, NekDouble> > &pacons_trace, int nqTrace)
    {
        int nparams=1;
        pacons_trace = Array<OneD, Array<OneD, NekDouble> >(nparams);

        pacons_trace[0] = Array<OneD, NekDouble>(nqTrace);
        Vmath::Vcopy(nqTrace,m_beta_trace[omega],1,pacons_trace[0],1);
    }




    /**
     *  Contains the definition of the pressure-area relationship. To implement a new pressure-area relationship
     *  change the mathematical definitions here and also in getdpda. Also consider changing: getAfromPressure
     *  and getCharIntegral; currently use numerical methods but it is preferable to enter analyticall expressions
     *  if they can be derived for the chosen p-A relation. Furthermore ensure that the definition of the lower
     *  integration boundary A_0 is updated in getAfromPressure, getAfromChars and getCharIntegral if using the
     *  numerical methods to ensure it refers to the correct element of pacons; typically can be assumed to be
     *  equal to the passive area.
     *
     *  The p-A relation cannot be a function of u and can contain any number of parameters as specified by
     *  P-A_num_Params in the input file which should then be specified for each subdomain using:
     *  <FUNCTION NAME="P-A_Vals">
     *  <E VAR="Val[Domain ID][Parameter ID]" VALUE="" />
     *  </FUNCTION>
     */
    
    void ArterialPressureArea::v_getPressure(NekDouble A, Array<OneD, NekDouble> pacons, NekDouble &pressure)
    {
        
        /*
         * pacons is an array which holds, for the point for which at which we are evaluating this function,
         * all of the constants found in the pressure-area relation. They are in the same order as the input file.
         */
        
        //pacons[0] = p_ext
        //pacons[1] = E
        //pacons[2] = A_0
        //pacons[3] = h0
        //pacons[4] = nue
        
        NekDouble Passive = 0.0;
        NekDouble Beta_Const = 0.0;
        
        Beta_Const = ((sqrt(PI)*pacons[3]*pacons[1])/((1-pow(pacons[4],2))*pacons[2]));
        
        Passive = (Beta_Const*(sqrt(A)-sqrt(pacons[2])));
        
        pressure = pacons[0] + Passive;
        
        //cout<<"Working?"<<pressure<<endl;
        
    }
    
    /**
     *  Contains the definition of dp/da - needs to be updated when chaning p-A relation!
     */
    void ArterialPressureArea::v_getdpda(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dpda)
    {
        NekDouble Passive = 0.0;
        NekDouble Beta_Const = 0.0;
        
        Beta_Const = ((sqrt(PI)*pacons[3]*pacons[1])/((1-pow(pacons[4],2))*pacons[2]));
        Passive = (Beta_Const*0.5*(1/sqrt(A)));
        
        dpda = Passive;
        
    }
    
    /**
     *  Calculates A in terms of pressure for for any p-A relation using a Newton iteration. If using a different
     *  p-A relation the definition of A_0 may need modifying.
     */
    void ArterialPressureArea::v_getAfromPressure(NekDouble pressure, Array<OneD, NekDouble> pacons, NekDouble &A)
    {
        
        //A = (pressure - pacons[2])/pacons[0] + sqrt(pacons[1]);
        //A = A*A;
        
        NekDouble fa = 0.0;
        NekDouble dpda = 0.0;
        NekDouble p_calc = 0.0;
        NekDouble delta_A_calc = 0.0;
        NekDouble A_0 = 0.0;
        
        // Initialise to A_0
        A_0=PI*pow(pacons[2],2)/40;
        A=A_0;
        
        int proceed = 1;
        int iter = 0;
        int MAX_ITER = 200;
        
        // Tolerances for the algorithm
        NekDouble Tol = 1.0e-10;
        
        while ((proceed) && (iter < MAX_ITER))
        {
            iter =iter+1;
            
            getPressure(A, pacons, p_calc);
            fa = p_calc-pressure;
            getdpda(A, 0, pacons, dpda);
            
            delta_A_calc=fa/dpda;
            A = A - delta_A_calc;
            
            if (sqrt(delta_A_calc*delta_A_calc) < Tol)
                proceed = 0;
        }
        
    }
    
    /**
     *  Calculates the two characteristic variables for any p-A relation, currently using numerical
     *  integration. Unlike lambda and the derivitives we calculate both in the same function to avoid having
     *  to do the numerical integration twice.
     */
    void ArterialPressureArea::v_getW(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &W1, NekDouble &W2)
    {
        NekDouble integralTerm=0.0;
        getCharIntegral(A, u, pacons, integralTerm);
        W1 = u + integralTerm;
        W2 = u - integralTerm;
    }
    
    /**
     *  Contains the definition of A in terms of W_1 and W_2 for any p-A relation though definition of A_0 may
     *  need updating. Solved for using a Newton iteration and numerical integration and numerical integration.
     */
    void ArterialPressureArea::v_getAfromChars(NekDouble W_1, NekDouble W_2, Array<OneD, NekDouble> pacons, NekDouble &A)
    {
        
        //    A = (W_1 - W_2)*sqrt(m_rho/(32*pacons[0]))+sqrt(sqrt(pacons[1]));
        //    A = A*A;
        //    A = A*A;
        
        NekDouble fa = 0.0;
        NekDouble dfa = 0.0;
        NekDouble dpda = 0.0;
        NekDouble integralTerm = 0.0;
        NekDouble delta_A_calc = 0.0;
        NekDouble A_0 = 0.0;
        
        // Initialise to A_0
        A_0=PI*pow(pacons[2],2)/40;
        A=A_0;
        
        int proceed = 1;
        int iter = 0;
        int MAX_ITER = 200;
        
        // Tolerances for the algorithm
        NekDouble Tol = 1.0e-10;
        
        while ((proceed) && (iter < MAX_ITER))
        {
            iter =iter+1;
            
            getCharIntegral(A, 0, pacons, integralTerm);
            fa = integralTerm - 0.5*(W_1-W_2);
            getdpda(A, 0, pacons, dpda);
            dfa=sqrt(dpda/(m_rho*A));
            
            delta_A_calc=fa/dfa;
            A = A - delta_A_calc;
            
            if (sqrt(delta_A_calc*delta_A_calc) < Tol)
                proceed = 0;
        }
        
    }
    
    /**
     *  Contains the definition of u in terms of W_1 and W_2 for any p-A relation
     */
    void ArterialPressureArea::v_getufromChars(NekDouble W_1, NekDouble W_2, Array<OneD, NekDouble> pacons, NekDouble &u)
    {
        
        u = 0.5*(W_1+W_2);
        
    }
    
    /**
     *  Contains the definition of the the first eigenvalue lambda_1 for any p-A relation
     */
    void ArterialPressureArea::v_getlambda1(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &lambda_1)
    {
        NekDouble dpda=0.0;
        getdpda(A, u, pacons, dpda);
        
        lambda_1 = u + sqrt((A/m_rho)*dpda);
        
    }
    
    /**
     *  Contains the definition of the the second eigenvalue lambda_2 for any p-A relation
     */
    void ArterialPressureArea::v_getlambda2(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &lambda_2)
    {
        NekDouble dpda=0.0;
        getdpda(A, u, pacons, dpda);
        
        lambda_2 = u - sqrt((A/m_rho)*dpda);
    }
    
    
    /**
     *  Contains the definition of dp/du for any p-A relation
     */
    void ArterialPressureArea::v_getdpdu(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dpdu)
    {
        
        dpdu = 0;
        
    }
    
    /**
     *  Contains the definition of dW1/da for any p-A relation
     */
    void ArterialPressureArea::v_getdW1da(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW1da)
    {
        NekDouble dpda=0.0;
        getdpda(A, u, pacons, dpda);
        
        dW1da = sqrt(dpda/(A*m_rho));
        
    }
    
    /**
     *  Contains the definition of dW1/du for any p-A relation
     */
    void ArterialPressureArea::v_getdW1du(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW1du)
    {
        
        dW1du = 1;
        
    }
    
    /**
     *  Contains the definition of dW2/da for any p-A relation
     */
    void ArterialPressureArea::v_getdW2da(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW2da)
    {
        
        NekDouble dpda=0.0;
        getdpda(A, u, pacons, dpda);
        
        dW2da = -1.0*sqrt(dpda/(A*m_rho));
        
    }
    
    /**
     *  Contains the definition of dW2/du for any p-A relation
     */
    void ArterialPressureArea::v_getdW2du(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &dW2du)
    {
        
        dW2du = 1;
        
    }
    
    /**
     *  Evalautes the integral \int^A_{A_{0}}{ \sqrt{\frac{1}{\rho A} \frac{\partial p\left(A,t\right)}{\partial A}} dA}
     *  which is found in the general definition of the characteristic variables.
     *  Currenly does this numerically but if p-A relation allows suggest the analytical expression is entered instead.
     *  This current implementation is suitable for any p-A relation, though A_0 may need updating as well as num of
     *  trapeziums.
     */
    void ArterialPressureArea::v_getCharIntegral(NekDouble A, NekDouble u, Array<OneD, NekDouble> pacons, NekDouble &intergralTerm)
    {
        
        int numTrapeziums = 50;
        NekDouble A_0 = 0.1*PI*pow(pacons[2],2)/4; //For active contractions can be less than static area
        NekDouble A_n = A_0;
        NekDouble dA = (A - A_n)/numTrapeziums;
        NekDouble dpda_n = 0.0;
        
        intergralTerm = 0.0;
        
        for (int i = 2; i < numTrapeziums+1; i++)
        {
            A_n = A_n +dA;
            getdpda(A_n, u, pacons, dpda_n);
            intergralTerm = intergralTerm + sqrt(dpda_n/(m_rho*A_n));
        }
        
        intergralTerm=intergralTerm*2.0;
        getdpda(A_0, u, pacons, dpda_n);
        intergralTerm=intergralTerm+sqrt(dpda_n/(m_rho*A_0));
        getdpda(A, u, pacons, dpda_n);
        intergralTerm=intergralTerm+sqrt(dpda_n/(m_rho*A));
        intergralTerm = intergralTerm*0.5*(A-A_0)/numTrapeziums;
        
    }
    
    /**
     *  Solves an Ax=b system using the linsys utility. Note that the A matrix must be passed as a 1D buffer array of
     *  rows which is then used to form a NekMatrix
     */
    void ArterialPressureArea::v_solveAxb(int rows, const Array<OneD, NekDouble> &matrix_buf, const Array<OneD, NekDouble> &b_buf, Array<OneD, NekDouble> &x)
    {
        
        //Cast A as NekMatrix and b as NekVector
        boost::shared_ptr<NekMatrix<double> > N(new NekMatrix<double>(rows, rows, matrix_buf,eFULL));
        boost::shared_ptr<NekVector<double> > b(new NekVector<double>(rows, b_buf));
        
        //Initialise Linear System
        LinearSystem linsys(N);
        
        //Solve
        NekVector<double> result = linsys.Solve(b);
        result = linsys.SolveTranspose(b);
        
        //Convert the resulting NekVector to an array of NekDouble's for further manipulations later
        for (int i = 0; i < rows; i++)
        {
            x[i]=result[i];
        }
        
    }

}
