///////////////////////////////////////////////////////////////////////////////
//
// File DriverUnsteadyBaseFlow.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SolverUtils/DriverUnsteadyBaseFlow.h>
#include <SolverUtils/AdvectionSystem.h>
namespace Nektar
{
    namespace SolverUtils
    {
        string DriverUnsteadyBaseFlow::className = GetDriverFactory().RegisterCreatorFunction("UnsteadyBaseFlow", DriverUnsteadyBaseFlow::create);
        string DriverUnsteadyBaseFlow::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","UnsteadyBaseFlow",0);

        /**
	 *
         */
        DriverUnsteadyBaseFlow::DriverUnsteadyBaseFlow(const LibUtilities::SessionReaderSharedPtr pSession)
            : Driver(pSession)
        {
        }
    
    
        /**
         *
         */
        DriverUnsteadyBaseFlow:: ~DriverUnsteadyBaseFlow()
        {
        }
    
    
        /**
         *
         */
        void DriverUnsteadyBaseFlow::v_InitObject(ostream &out)
        {
            Driver::v_InitObject(out);
        }
    
    
        void DriverUnsteadyBaseFlow::v_Execute(ostream &out)
        
        {
            time_t starttime, endtime;
            NekDouble CPUtime;
            NekDouble NumVar = m_equ[0]->GetNvariables();
            m_dt = m_equ[0]->GetTimeStep();
            
            
            string advName;
            m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
            m_session->LoadParameter ("numCheck",   m_numCheck,   1.0);
            
            //AdvectionSystemSharedPtr A = m_equ[0]->as<AdvectionSystem>();
            AdvectionSharedPtr A = SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
            NumVar = m_equ[0]->UpdateFields()[0]->GetCoordim(0);
            if (m_session->GetSolverInfo("EqType") == "AdjointEulerCFE" ||
                m_session->GetSolverInfo("EqType") == "AdjointNavierStokesCFE")
            {
                // Number of variables for the compressible equations
                NumVar += 2;
            }
            if(m_session->DefinesSolverInfo("HOMOGENEOUS"))
            {
                if (m_session->GetSolverInfo("HOMOGENEOUS") == "1D")
                {
                    NumVar += 1;
                }
            }
            
            Array<OneD, MultiRegions::ExpListSharedPtr> fields
            = m_equ[0]->UpdateFields();
            
            time(&starttime);
            
            for (int i = 0; i < m_numCheck; ++i)
            {
                m_equ[0]->PrintSummary(out);
                
                m_equ[0]->DoInitialise();
                
                m_equ[0]->DoSolve();
                
                m_equ[0]->Checkpoint_Output(i+1);
                
            }
            time(&endtime);
            
            if (m_comm->GetRank() == 0)
            {
                CPUtime = difftime(endtime, starttime);
                cout << "-------------------------------------------" << endl;
                cout << "Total Computation Time = " << CPUtime << "s" << endl;
                cout << "-------------------------------------------" << endl;
            }
                
                // Evaluate and output computation time and solution accuracy.
                // The specific format of the error output is essential for the
                // regression tests to work.
                // Evaluate L2 Error
            for(int i = 0; i < m_equ[0]->GetNvariables(); ++i)
            {
                    Array<OneD, NekDouble> exactsoln(m_equ[0]->GetTotPoints(), 0.0);
                    
                    // Evaluate "ExactSolution" function, or zero array
                m_equ[0]->EvaluateExactSolution(i, exactsoln,
                                                    m_equ[0]->GetFinalTime());
                    
                NekDouble vL2Error   = m_equ[0]->L2Error  (i, exactsoln);
                NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln);
                    
                if (m_comm->GetRank() == 0)
                {
                        out << "L 2 error (variable " << m_equ[0]->GetVariable(i)
                        << ") : " << vL2Error << endl;
                        out << "L inf error (variable " << m_equ[0]->GetVariable(i)
                    << ") : " << vLinfError << endl;
                }
            }
            
        }
        void  DriverUnsteadyBaseFlow::v_ImportFldUnsteadyBase(
                          std::string pInfile,
                          Array<OneD, Array<OneD, NekDouble> > q0,
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                          LibUtilities::SessionReaderSharedPtr        pSession)
        {
            int nvar = pSession->GetVariables().size();
            int s;
            
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;
            int nqtot = q0[0].num_elements();
            
            Array<OneD, NekDouble> tmp_coeff(pFields[0]->GetNcoeffs(), 0.0);
            
            //Get Homogeneous
            LibUtilities::FieldIOSharedPtr fld =
            MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(
                                                           pSession->GetComm());
            fld->Import(pInfile, FieldDef, FieldData);
            
            
            if(pSession->DefinesSolverInfo("HOMOGENEOUS"))
            {
                std::string HomoStr = pSession->GetSolverInfo("HOMOGENEOUS");
            }
            
            // copy FieldData into m_fields
            for(int j = 0; j < nvar; ++j)
            {
                for(int i = 0; i < FieldDef.size(); ++i)
                {
                    if((pSession->DefinesSolverInfo("HOMOGENEOUS") &&
                        (pSession->GetSolverInfo("HOMOGENEOUS")=="HOMOGENEOUS1D" ||
                         pSession->GetSolverInfo("HOMOGENEOUS")=="1D" ||
                         pSession->GetSolverInfo("HOMOGENEOUS")=="Homo1D")) &&
                       m_MultipleModes==false)
                    {
                        // w-component must be ignored and set to zero.
                        if (j != nvar - 2)
                        {
                            // p component (it is 4th variable of the 3D and corresponds 3nd variable of 2D)
                            s = (j == nvar - 1) ? 2 : j;
                            
                            //extraction of the 2D
                            pFields[j]->ExtractDataToCoeffs(
                                                FieldDef[i],
                                                FieldData[i],
                                                FieldDef[i]->m_fields[s],
                                                tmp_coeff);
                            
                        }
                        
                        //Put zero on higher modes
                        int ncplane = (pFields[0]->GetNcoeffs()) / m_npointsZ;
                        
                        if (m_npointsZ > 2)
                        {
                            Vmath::Zero(ncplane*(m_npointsZ-2),
                                        &tmp_coeff[2*ncplane], 1);
                        }
                    }
                    // 2D cases and Homogeneous1D Base Flows
                    else
                    {
                        bool flag = FieldDef[i]->m_fields[j] ==
                        pSession->GetVariable(j);
                        
                        ASSERTL0(flag, (std::string("Order of ") + pInfile
                                        + std::string(" data and that defined in "
                                                      "m_boundaryconditions differs")).c_str());
                        
                        pFields[j]->ExtractDataToCoeffs(FieldDef[i],
                                                        FieldData[i],
                                                        FieldDef[i]->m_fields[j],
                                                        tmp_coeff);
                    }
                }
                
                if(m_SingleMode || m_HalfMode)
                {
                    //pFields[j]->SetWaveSpace(true);
                    
                    pFields[j]->GetPlane(0)->BwdTrans(tmp_coeff, q0[j]);
                    
                    if(m_SingleMode)
                    {
                        //copy the bwd into the second plane for single Mode Analysis
                        int ncplane=(pFields[0]->GetNpoints())/m_npointsZ;
                        Vmath::Vcopy(ncplane,&q0[j][0],1,&q0[j][ncplane],1);
                    }
                }
                else
                {
                    pFields[j]->BwdTrans(tmp_coeff, q0[j]);
                }
            }
        }
    }
}

