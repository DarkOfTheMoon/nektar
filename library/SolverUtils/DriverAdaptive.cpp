///////////////////////////////////////////////////////////////////////////////
//
// File DriverAdaptive.cpp
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
// Description: Driver class for adaptive solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SolverUtils/DriverAdaptive.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>
//#include <GlobalMapping/Mapping.h>

namespace Nektar
{
    namespace SolverUtils
    {
        string DriverAdaptive::className = GetDriverFactory().RegisterCreatorFunction("Adaptive", DriverAdaptive::create);
        string DriverAdaptive::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","Adaptive",0);

        /**
	 *
         */
        DriverAdaptive::DriverAdaptive(const LibUtilities::SessionReaderSharedPtr pSession)
            : Driver(pSession)
        {
        }
    
    
        /**
         *
         */
        DriverAdaptive:: ~DriverAdaptive()
        {
        }
    
    
        /**
         *
         */
        void DriverAdaptive::v_InitObject(ostream &out)
        {
            Driver::v_InitObject(out);
            
        }
    
    
        void DriverAdaptive::v_Execute(ostream &out)        
        {
            time_t starttime, endtime;
            NekDouble CPUtime;
            
            m_equ[0]->PrintSummary(out);

            // First run using original order
            time(&starttime);
            m_equ[0]->DoInitialise();
            // Obtain initial time in case a restart was used
            NekDouble startTime = m_equ[0]->GetFinalTime();
            m_equ[0]->DoSolve();
            
            // Get Information            
            bool isHomogeneous1D;
            int nRuns, minP, maxP, sensorVar;
            NekDouble lowerTol, upperTol;
            m_session->LoadParameter("NumRuns",  nRuns, 1);
            m_session->LoadParameter("AdaptiveMinModes",  minP, 4);
            m_session->LoadParameter("AdaptiveMaxModes",  maxP, 12);
            m_session->LoadParameter("AdaptiveLowerTolerance",  lowerTol, 1e-8);
            m_session->LoadParameter("AdaptiveUpperTolerance",  upperTol, 1e-6);
            m_session->LoadParameter("AdaptiveSensorVariable",  sensorVar, 0);
            m_session->MatchSolverInfo("Homogeneous", "1D",
                                       isHomogeneous1D, false);
            
            // Get number of elements and planes
            int nExp, nPlanes;
            if (isHomogeneous1D)
            {
                nExp = m_equ[0]->UpdateFields()[0]->GetPlane(0)->GetExpSize();
                nPlanes = m_equ[0]->UpdateFields()[0]->GetZIDs().num_elements();
            }
            else
            {
                nExp = m_equ[0]->UpdateFields()[0]->GetExpSize();
                nPlanes = 1;
            }
            
            int nFields = m_equ[0]->UpdateFields().num_elements();
            
            NekDouble period  = m_session->GetParameter("TimeStep")* 
                                m_session->GetParameter("NumSteps");
            
            Array< OneD, NekDouble> coeffs;
            Array< OneD, NekDouble> phys; 
            Array< OneD, NekDouble> physReduced;
            Array< OneD, NekDouble> tmpArray;
            //GlobalMapping::MappingSharedPtr mapping;
            // Adaptive loop
            LocalRegions:: ExpansionSharedPtr Exp;
            for( int i = 1; i < nRuns; i++ )
            {
                // Get field expansions
                Array<OneD, MultiRegions::ExpListSharedPtr> fields = 
                                                    m_equ[0]->UpdateFields();

                // Determine the change to be applied in the order
                map<int, int> deltaP;
                int offset = 0;
                for ( int n = 0; n < nExp; n++)
                {
                    offset = fields[sensorVar]->GetPhys_Offset(n);

                    Exp = fields[sensorVar]->GetExp(n);

                    int P = Exp->GetBasis(0)->GetNumModes();
                    int Q = Exp->GetBasis(1)->GetNumModes();

                    int qa = Exp->GetBasis(0)->GetNumPoints();
                    int qb = Exp->GetBasis(1)->GetNumPoints();

                    // Declare orthogonal basis.
                    LibUtilities::PointsKey pa(qa,Exp->
                                                GetBasis(0)->GetPointsType());
                    LibUtilities::PointsKey pb(qb,Exp->
                                                GetBasis(1)->GetPointsType());
                    
                    StdRegions::StdExpansionSharedPtr OrthoExp;
                    switch (Exp->GetGeom()->GetShapeType())
                    {
                        case LibUtilities::eQuadrilateral:
                        {
                            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A,P-1,pa);
                            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A,Q-1,pb);
                            OrthoExp = MemoryManager<StdRegions::StdQuadExp>::
                                            AllocateSharedPtr(Ba,Bb);
                        }
                        break;
                        case LibUtilities::eTriangle:
                        {
                            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A,P-1,pa);
                            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B,Q-1,pb);
                            OrthoExp = MemoryManager<StdRegions::StdTriExp>::
                                            AllocateSharedPtr(Ba,Bb);
                        }
                        break;
                        default:
                            ASSERTL0(false, "Shape not supported.");
                            break;
                    }
                    int nq = OrthoExp->GetTotPoints();

                    NekDouble error = 0;
                    NekDouble fac = 0;
                    NekDouble tmp = 0;

                    coeffs = Array< OneD, NekDouble> (OrthoExp->GetNcoeffs());
                    physReduced = Array< OneD, NekDouble> (OrthoExp->GetTotPoints());
                    tmpArray = Array< OneD, NekDouble> (OrthoExp->GetTotPoints(), 0.0);

                    // Refinement based only on one variable
                    for ( int plane = 0; plane < nPlanes; plane++)
                    {
                        if (isHomogeneous1D)
                        {
                            phys = fields[sensorVar]->GetPlane(plane)->GetPhys()+offset;
                        }
                        else
                        {
                            phys = fields[sensorVar]->GetPhys()+offset;
                        }
                        // Project solution to lower order
                        OrthoExp->FwdTrans(phys,coeffs);
                        OrthoExp->BwdTrans(coeffs,physReduced);

                        // Calculate error =||phys-physReduced||^2 / ||phys||^2
                        Vmath::Vsub( nq, phys, 1, physReduced, 1, tmpArray, 1);
                        Vmath::Vmul( nq, tmpArray, 1, tmpArray, 1, tmpArray, 1);
                        tmp = Exp->Integral(tmpArray);

                        Vmath::Vmul( nq, phys, 1, phys, 1, tmpArray, 1);
                        fac = Exp->Integral(tmpArray);

                        tmp = abs(tmp/fac);

                        if (tmp != tmp)
                        {
                            ASSERTL0( false,
                                "Adaptive procedure encountered Nan value.");
                        }

                        // Get maximum value along planes
                        error = ( tmp > error) ? tmp : error;
                    }
                    
                    // Combine planes from different processes
                    m_session->GetComm()->GetColumnComm()->AllReduce
                                        (error, LibUtilities::ReduceMax);
                    
                    // Override tolerances if function is defined
                    if (m_session->DefinesFunction("AdaptiveLowerTolerance"))
                    {
                        int nq = Exp->GetTotPoints();
                        // Obtain points from the the element
                        Array<OneD,NekDouble>  xc0,xc1,xc2;
                        xc0 = Array<OneD,NekDouble>(nq,0.0);
                        xc1 = Array<OneD,NekDouble>(nq,0.0);
                        xc2 = Array<OneD,NekDouble>(nq,0.0);
                        Exp->GetCoords(xc0,xc1,xc2);

                        // Evaluate function from session file
                        Array< OneD, NekDouble> tolerance(nq,0.0);
                        LibUtilities::EquationSharedPtr ffunc = 
                                m_session->GetFunction("AdaptiveLowerTolerance", 0);
                        ffunc->Evaluate(xc0,xc1,xc2,tolerance);
                        lowerTol = Vmath::Vsum(nq, tolerance, 1)/nq;
                    }
                    if (m_session->DefinesFunction("AdaptiveUpperTolerance"))
                    {
                        int nq = Exp->GetTotPoints();
                        // Obtain points from the the element
                        Array<OneD,NekDouble>  xc0,xc1,xc2;
                        xc0 = Array<OneD,NekDouble>(nq,0.0);
                        xc1 = Array<OneD,NekDouble>(nq,0.0);
                        xc2 = Array<OneD,NekDouble>(nq,0.0);
                        Exp->GetCoords(xc0,xc1,xc2);

                        // Evaluate function from session file
                        Array< OneD, NekDouble> tolerance(nq,0.0);
                        LibUtilities::EquationSharedPtr ffunc = m_session->GetFunction("AdaptiveUpperTolerance", 0);
                        ffunc->Evaluate(xc0,xc1,xc2,tolerance);
                        upperTol = Vmath::Vsum(nq, tolerance, 1)/nq;
                    }

                    // Determine what to do with the polynomial order
                    if ( (error > upperTol) &&
                         (P < maxP)   )
                    {
                        deltaP[Exp->GetGeom()->GetGlobalID()] = 1;
                    }
                    else if( (error < lowerTol) &&
                              P > minP )
                    {
                        deltaP[Exp->GetGeom()->GetGlobalID()] = -1;
                    }
                    else
                    {
                        deltaP[Exp->GetGeom()->GetGlobalID()] =  0;
                    }
                }

                // Get mapping (not in master yet...)
                // mapping = GlobalMapping::Mapping::Load(m_session, 
                //                                      m_equ[0]->UpdateFields()); 
                
                // Write new expansion section to the session reader
                ReplaceExpansion(fields, deltaP); 
                // Reset GlobalLinSys Manager to avoid using too much memory
                LibUtilities::NekManager<MultiRegions::GlobalLinSysKey, 
                                    MultiRegions::GlobalLinSys>::ClearManager
                                    (std::string("GlobalLinSys"));
                
                // Initialise driver again
                Driver::v_InitObject(out);

                // Initialise equation
                m_equ[0]->SetTime(startTime + i*period);
                m_equ[0]->SetBoundaryConditions(startTime + i*period);
                m_equ[0]->SetInitialConditions(startTime + i*period, false);
                
                // Project solution to new expansion 
                //    (ExtractCoeffsToCoeffs would be simpler, but did not work)
                std::vector<std::string> fieldNames = m_session->GetVariables();
                std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                            = fields[0]->GetFieldDefinitions();
                std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
                for (int j = 0; j < FieldDef.size(); ++j)
                {
                    for (int n = 0; n<nFields; n++)
                    {
                        FieldDef[j]->m_fields.push_back(fieldNames[n]);
                        fields[n]->AppendFieldData(FieldDef[j], FieldData[j]);
                    }
                }                                
                for (int n = 0; n<nFields; n++)
                {                   
                    for (int j = 0; j < FieldDef.size(); ++j)
                    {
                        m_equ[0]->UpdateFields()[n]->ExtractDataToCoeffs(
                                FieldDef[j], FieldData[j],
                                fieldNames[n], 
                                m_equ[0]->UpdateFields()[n]->UpdateCoeffs());
                    }
                    m_equ[0]->UpdateFields()[n]->BwdTrans_IterPerExp(
                                    m_equ[0]->UpdateFields()[n]->GetCoeffs(),
                                    m_equ[0]->UpdateFields()[n]->UpdatePhys());
                }
                
                // Update mapping
                //mapping->ReplaceField(m_equ[0]->UpdateFields());
                
                // Solve equation                
                m_equ[0]->DoSolve();   
            }

            time(&endtime);

            m_equ[0]->Output();
        
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
        
        void DriverAdaptive::ReplaceExpansion(
                        Array<OneD, MultiRegions::ExpListSharedPtr>& fields,
                        map<int, int> deltaP) 
        {
            int nExp, nDim;      
            
            // Get field definitions
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefs
                        = fields[0]->GetFieldDefinitions();
            int expDim = 2;
            if (fielddefs[0]->m_numHomogeneousDir == 1)
            {
                nDim = 3;
            }
            else
            {
                nDim = 2;
            }
            // Add variables to field definition
            for(int i = 0; i < fielddefs.size(); ++i)
            {
                for (int j = 0; j < fields.num_elements(); ++j)
                {
                    fielddefs[i]->m_fields.push_back(m_session->GetVariable(j));
                }
            }         
            
            // Get tinyxml objects
            TiXmlElement* exp_tag = m_session->GetElement("NEKTAR/EXPANSIONS");
            // Clear current expansions
            exp_tag->Clear();
            
            // Write new expansion information
            for (int f = 0; f < fielddefs.size(); ++f)
            {
                nExp = fielddefs[f]->m_elementIDs.size();
                //---------------------------------------------
                // Write ELEMENTS
                TiXmlElement * elemTag = new TiXmlElement("ELEMENTS");
                exp_tag->LinkEndChild(elemTag);

                // Write FIELDS
                std::string fieldsString;
                {
                    std::stringstream fieldsStringStream;
                    bool first = true;
                    for (std::vector<int>::size_type i = 0; i
                            < fielddefs[f]->m_fields.size(); i++)
                    {
                        if (!first)
                            fieldsStringStream << ",";
                        fieldsStringStream << fielddefs[f]->m_fields[i];
                        first = false;
                    }
                    fieldsString = fieldsStringStream.str();
                }
                elemTag->SetAttribute("FIELDS", fieldsString);

                // Write SHAPE
                std::string shapeString;
                {
                    std::stringstream shapeStringStream;
                    shapeStringStream << LibUtilities::ShapeTypeMap[fielddefs[f]->m_shapeType];
                    if(fielddefs[f]->m_numHomogeneousDir == 1)
                    {
                        shapeStringStream << "-HomogenousExp1D";
                    }
                    else if (fielddefs[f]->m_numHomogeneousDir == 2)
                    {
                        shapeStringStream << "-HomogenousExp2D";
                    }

                    shapeString = shapeStringStream.str();
                }
                elemTag->SetAttribute("SHAPE", shapeString);

                // Write BASIS
                std::string basisString;
                {
                    std::stringstream basisStringStream;
                    bool first = true;
                    for (std::vector<LibUtilities::BasisType>::size_type i = 0; i < fielddefs[f]->m_basis.size(); i++)
                    {
                        if (!first)
                            basisStringStream << ",";
                        basisStringStream
                        << LibUtilities::BasisTypeMap[fielddefs[f]->m_basis[i]];
                        first = false;
                    }
                    basisString = basisStringStream.str();
                }
                elemTag->SetAttribute("BASIS", basisString);

                // Write homogeneuous length details
                if(fielddefs[f]->m_numHomogeneousDir)
                {
                    std::string homoLenString;
                    {
                        std::stringstream homoLenStringStream;
                        bool first = true;
                        for (int i = 0; i < fielddefs[f]->m_numHomogeneousDir; ++i)
                        {
                            if (!first)
                                homoLenStringStream << ",";
                            homoLenStringStream
                            << fielddefs[f]->m_homogeneousLengths[i];
                            first = false;
                        }
                        homoLenString = homoLenStringStream.str();
                    }
                    elemTag->SetAttribute("HOMOGENEOUSLENGTHS", homoLenString);
                }
				
                // Write homogeneuous planes/lines details
                if(fielddefs[f]->m_numHomogeneousDir)
                {
                    if(fielddefs[f]->m_homogeneousYIDs.size() > 0)
                    {
                        std::string homoYIDsString;
                        {
                            std::stringstream homoYIDsStringStream;
                            bool first = true;
                            for(int i = 0; i < fielddefs[f]->m_homogeneousYIDs.size(); i++)
                            {
                                if (!first)
                                    homoYIDsStringStream << ",";
                                homoYIDsStringStream << fielddefs[f]->m_homogeneousYIDs[i];
                                first = false;
                            }
                            homoYIDsString = homoYIDsStringStream.str();
                        }
                        elemTag->SetAttribute("HOMOGENEOUSYIDS", homoYIDsString);
                    }
                    
                    if(fielddefs[f]->m_homogeneousZIDs.size() > 0)
                    {
                        std::string homoZIDsString;
                        {
                            std::stringstream homoZIDsStringStream;
                            bool first = true;
                            for(int i = 0; i < fielddefs[f]->m_homogeneousZIDs.size(); i++)
                            {
                                if (!first)
                                    homoZIDsStringStream << ",";
                                homoZIDsStringStream << fielddefs[f]->m_homogeneousZIDs[i];
                                first = false;
                            }
                            homoZIDsString = homoZIDsStringStream.str();
                        }
                        elemTag->SetAttribute("HOMOGENEOUSZIDS", homoZIDsString);
                    }
                }
                
                // Write NUMMODESPERDIR
                std::string numModesString;
                {
                    std::stringstream numModesStringStream;

                    numModesStringStream << "MIXORDER:";
                    bool first = true;
                    int eId;
                    Array<OneD, int> order(nDim, 0);
                    for (int n = 0 ; n < nExp; n++)
                    {
                        eId = fielddefs[f]->m_elementIDs[n];
                        for (int i =0; i<expDim; i++)
                        {
                            order[i] = deltaP[eId];
                        }
                        for (int i =0; i<nDim; i++)
                        {
                            if (fielddefs[f]->m_uniOrder)
                            {
                                order[i] += fielddefs[f]->m_numModes[i];
                            }
                            else
                            {
                                order[i] += fielddefs[f]->m_numModes[n*nDim + i];
                            }
                            if (!first)
                                numModesStringStream << ",";
                            numModesStringStream << order[i];
                            first = false;                           
                        }
                    }
                    numModesString = numModesStringStream.str();
                }
                elemTag->SetAttribute("NUMMODESPERDIR", numModesString);

                // Write ID
                // Should ideally look at ways of compressing this stream
                // if just sequential;
                std::string idString;
                {
                    std::stringstream idStringStream;
                    GenerateSeqString(fielddefs[f]->m_elementIDs,idString);
                }
                elemTag->SetAttribute("ID", idString);
            }
        }
            
        void DriverAdaptive::GenerateSeqString(const std::vector<unsigned int> &elmtids,
                                      std::string &idString)
        {
            std::stringstream idStringStream;
            bool setdash = true;
            unsigned int endval;

            idStringStream << elmtids[0];
            for (int i = 1; i < elmtids.size(); ++i)
            {
                if(elmtids[i] == elmtids[i-1]+1)
                {
                    if(setdash)
                    {
                        idStringStream << "-";
                        setdash = false;
                    }

                    if(i == elmtids.size()-1) // last element
                    {
                        idStringStream << elmtids[i];
                    }
                    else
                    {
                        endval = elmtids[i];
                    }
                }
                else
                {
                    if(setdash == false) // finish off previous dash sequence
                    {
                        idStringStream << endval;
                        setdash = true;
                    }


                    idStringStream << "," << elmtids[i];
                }
            }
            idString = idStringStream.str();
        }
        
    }
}

