////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessNitricOxide.cpp
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
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EWPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes endothelial nitric oxides.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <sstream>
using namespace std;

#include "ProcessNitricOxide.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessNitricOxide::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "NO"), 
                ProcessNitricOxide::create, "Computes endothelial cell nitric-oxide.");

        ProcessNitricOxide::ProcessNitricOxide(FieldSharedPtr f) : ProcessModule(f)
        {

        }
        ProcessNitricOxide::~ProcessNitricOxide()
        {
        }
        
        void ProcessNitricOxide::Process(po::variables_map &vm)
        {
            if (m_f->m_verbose)
            {
                cout << "ProcessNitricOxide: Calculating eNOS..." << endl;
            }

            int i,j,k;
            int nfields = m_f->m_fielddef[0]->m_fields.size();
            int NumHomogeneousDir = m_f->m_fielddef[0]->m_numHomogeneousDir;

            // Set up Expansion information to use mode order from field
            m_f->m_graph->SetExpansions(m_f->m_fielddef);
            
            //Set up expansions, and extract data. 
            m_f->m_exp.resize(nfields);
            m_f->m_exp[0] = m_f->SetUpFirstExpList(NumHomogeneousDir,true);
            
            cout << "Test 3" << endl;
            for(j = 1; j < nfields; ++j)
            {
                m_f->m_exp[j] = m_f->AppendExpList(NumHomogeneousDir);
            }
            
            for (j = 0; j < nfields; ++j)
            {
                for (int k = 0; k < m_f->m_data.size(); ++k)
                {
                    m_f->m_exp[j]->ExtractDataToCoeffs(
                        m_f->m_fielddef[k],
                        m_f->m_data[k],
                        m_f->m_fielddef[k]->m_fields[j],
                        m_f->m_exp[j]->UpdateCoeffs());
                }
                m_f->m_exp[j]->BwdTrans(m_f->m_exp[j]->GetCoeffs(),
                                                   m_f->m_exp[j]->UpdatePhys());
            }

 
            cout << "Test 4" << endl;
            int spacedim   = m_f->m_graph->GetSpaceDimension();
            if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
                (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
            {
                spacedim = 3;
            }


	        int npoints = m_f->m_exp[0]->GetNpoints();
            int nout = 2;	  
	        std::cout << npoints << endl;

            Array<OneD, Array<OneD, NekDouble> > outfield(nout);
            Array<OneD, NekDouble> wss(npoints);
            
            for (i = 0; i < nout; ++i)
            {
                outfield[i] = Array<OneD, NekDouble>(npoints);
                Vmath::Zero(npoints, outfield[i],1);
            }


            cout << "Test 5" << endl;
            //Get WSS - Assumes the input is from multishear
            wss = m_f->m_exp[0]->GetPhys();
            
            // Global Variables 

            int n = 4;
            NekDouble dt = 0.001;
            NekDouble Ca0 = 0.313;
            NekDouble C_s0 = 8840.;
            NekDouble Bt = 375.;
            NekDouble kb = 9.38;
            NekDouble phi = 0.9;
            

            NekDouble steps = 1000;

            Array<OneD, NekDouble> k_1(n), k_2(n), k_3(n), k_4(n);
            Array<OneD, NekDouble> wi(n);
            Array<OneD, NekDouble> w(n,0.0),dw(n,0.0);
            
   
        
            //RK 4
            
            for (i=0;i<npoints;++i)
            {

                //Initialize Solution
                wi[0] = 0.0;
                wi[1] = Ca0*(Ca0 + Bt + kb)/(Ca0 + kb);
                wi[2] = C_s0;  
                wi[3] = 0.0;    

                for (j=0;j<steps;++j)
                {
                    dw = Rx(wi,n,wss[i]);
                    for (k=0;k<n;++k)
                    {
                        k_1[k] = dw[k]*dt;
                        w[k] = wi[k] + 0.5*k_1[k];
                    }

                    dw = Rx(w,n,wss[i]);
                    for (k=0;k<n;++k)
                    {
                        k_2[k] = dw[k]*dt;
                        w[k] = wi[k] + 0.5*k_2[k];
                    }

                    dw = Rx(w,n,wss[i]);
                    for (k=0;k<n;++k)
                    {
                        k_3[k] = dw[k]*dt;
                        w[k] = wi[k] + 0.5*k_3[k];
                    }
                
                    dw = Rx(w,n,wss[i]);
                    for (k=0;k<n;++k)
                    {
                        k_4[k] = dw[k]*dt;
                        w[k] = wi[k] + 1./6.*(k_1[k] + 2.*k_2[k] + 2.*k_3[k] + k_4[k]);
                    }

                    for(k=0;k<n;++k)
                    {
                        wi[k] = w[k];
                    }

                }
                cout << i << " "<< npoints << endl; 
                outfield[1][i] = w[3];
            }

         
            for (i = 0;i < npoints; ++i)
            {
                outfield[0][i] = wss[i];
            }

            m_f->m_exp.resize(nout);
            m_f->m_fielddef = m_f->m_fielddef;
            m_f->m_exp[0] = m_f->SetUpFirstExpList(m_f->m_fielddef[0]->m_numHomogeneousDir,true);
            
            for(i = 1; i < nout; ++i)
            {
                m_f->m_exp[i] = m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
            }                
            
            m_f->m_fielddef[0]->m_fields.resize(nout);
            m_f->m_fielddef[0]->m_fields[0] = "TAWSS";
            m_f->m_fielddef[0]->m_fields[1] = "NO";
                       
            for(i = 0; i < nout; ++i)
            {
	    	cout << i << endl;
                m_f->m_exp[i]->FwdTrans(outfield[i],
                                        m_f->m_exp[i]->UpdateCoeffs());
                m_f->m_exp[i]->BwdTrans(m_f->m_exp[i]->GetCoeffs(), 
                                        m_f->m_exp[i]->UpdatePhys());
            }


            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = m_f->m_exp[0]->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
            
            for( i = 0; i < nout; ++i)
            {
                for ( j = 0; j < FieldDef.size(); ++j)
                {
                    FieldDef[j]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[i]);
                    m_f->m_exp[i]->AppendFieldData(FieldDef[j], FieldData[j]);
                }
            }
            
            m_f->m_fielddef = FieldDef;
            m_f->m_data     = FieldData;

	    }


        Array<OneD, NekDouble> ProcessNitricOxide::Rx(Array<OneD, NekDouble> &W, int n, NekDouble Tau)
        {
            //Reaction Parameter ND
            NekDouble k1 = 8.53;
            NekDouble k2 = 100.;
            NekDouble k3 = 3320.;
            NekDouble k4 = 7810.;
            NekDouble k5 = 1.6e-5;
            NekDouble k6 = 100.*500./0.32;
            NekDouble k7 = 300.*500.;
            NekDouble k8 = 38600.;
            NekDouble kb = 9.38;
            NekDouble krel = .64*500.;
            NekDouble kdis = 0.09*500./0.32;
            NekDouble mhu2 = 0.0167*500.;
            NekDouble qmax = 27500.;
            NekDouble gmax = 0.06*500./0.32;
            NekDouble kCCE = 1.28e-4;
            NekDouble C_b0 = 3.9/0.32;
            NekDouble Ca0 = 0.313;
            NekDouble C_s0 = 8840.;
            NekDouble C_ex = 4690.;
            NekDouble Bt = 375.;
            NekDouble KCICR = 0.0;
            NekDouble K1 = 0.0;
            NekDouble K2 = 0.625;
            NekDouble K3 = 0.469;
            NekDouble K4 = 1.0;
            NekDouble K5 = 0.45/0.32;
            NekDouble Kc = 0.26;
            NekDouble beta = 2.63;
            NekDouble W0 = 111.;
            NekDouble fe = 0.0134;
            NekDouble alpha = 2.;
            NekDouble omega = 0.5;
            NekDouble epsilon = 0.1;
            NekDouble Vr = 3.5;
            NekDouble phi = 0.9;

            NekDouble Shear,FShear;

            int i;
            Array<OneD, NekDouble> dW(n);
            
            
        	W[1] = 0.5*(W[1] - Bt - kb + sqrt(pow((W[1] - Bt - kb),2) + 4.*kb*W[1]));

            Shear = fe*W0*pow((epsilon*Tau + sqrt(16.*pow(beta,2) + pow(epsilon,2)*pow(Tau,2)) - 4.*beta),2)/(epsilon*Tau + sqrt(16.*pow(beta,2) + pow(epsilon,2)*pow(Tau,2)));
            FShear = 1./(1.+alpha*exp(-Shear));

            dW[0] = k1*phi/(Kc + phi)*W[1]/(K1 + W[1]) - k2*W[0];
            dW[1] = k3*W[2]*pow((W[0]/(K2+W[0])),3) - k4*pow((W[1]/(K3 + W[1])),2) + kCCE*(C_s0-W[2])*(C_ex-W[1]) + k5*pow(W[2],2) + omega*qmax*FShear - k8*W[1]/(K4 + W[1]); 
            dW[2] = -Vr*(k3*pow((W[0]/(K2+W[0])),3)*W[2] - k4*pow((W[1]/(K3 + W[1])),2) + k5*pow(W[2],2));
            dW[3] = kdis*W[1]/(K5 + W[1]) - mhu2*W[3] + (1.-omega)*gmax*FShear;
            

            return dW;

        }
    }
}           


