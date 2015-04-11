////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessMultiShear.cpp
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
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes tawss, osi, transwss, afi, cfi fields.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <sstream>
using namespace std;

#include "ProcessMultiShear.h"

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
        //    m_config["N"] = ConfigOption(false,"1","Number of chk or fld files");
            m_config["fromfld"] = ConfigOption(false, "NotSet",
                                               "First fld file. First underscore flags position of id in name.");

            ASSERTL0(m_config["fromfld"].as<string>().compare("NotSet") != 0,
                     "Need to specify fromfld=file.fld ");

           m_f->m_fldToBnd = false;
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

            //May not need any of this
            int nstart, i, j, nfields;
            NekDouble nfld =  1 //m_config["N"].as<NekDouble>();
            string fromfld, basename, endname, nstartStr;
            stringstream filename;
            vector<string> infiles(nfld);
            vector< boost::shared_ptr<Field> > m_fromField(nfld);



//            StdRegions::StdExpansion2DSharedPtr bc;
//            StdRegions::StdExpansionSharedPtr elmt;

            // Set up list of input fld files. 
            fromfld = m_config["fromfld"].as<string>();
            basename = fromfld.substr(0, fromfld.find_first_of("_")+1);
            filename << fromfld.substr(fromfld.find_first_of("_")+1, fromfld.size());
            filename >> nstart;
            filename.str("");
            filename << nstart;
            filename >> nstartStr;
            filename.str("");
            endname = fromfld.substr(fromfld.find(nstartStr)+nstartStr.size(), fromfld.size());
       
   
            for (i=0; i<nfld; ++i)
            {
                stringstream filename;
                filename << basename << i+nstart << endname;
                filename >> infiles[i];       
                cout << infiles[i]<<endl;
            }


            for ( i = 0; i<nfld; ++i)
            {
                m_fromField[i] = boost::shared_ptr<Field>(new Field());
                m_fromField[i]->m_session = m_f->m_session;
                m_fromField[i]->m_graph = m_f->m_graph;
                m_fromField[i]->m_fld = MemoryManager<LibUtilities::FieldIO>
                    ::AllocateSharedPtr(m_fromField[0]->m_session->GetComm());
            }
            
            //Import all fld files. 
            for (i=0; i<nfld; ++i)
            {
                if(m_f->m_exp.size())
                {
                    // Set up ElementGIDs in case of parallel processing
                    Array<OneD,int> ElementGIDs(m_f->m_exp[0]->GetExpSize());
                    for (j = 0; j < m_f->m_exp[0]->GetExpSize(); ++j)
			    {
                        ElementGIDs[j] = m_f->m_exp[0]->GetExp(j)->GetGeom()->GetGlobalID();
                    }
                    m_fromField[i]->m_fld->Import(infiles[i],m_fromField[i]->m_fielddef,
                                                  m_fromField[i]->m_data,
                                                  LibUtilities::NullFieldMetaDataMap,
                                                  ElementGIDs);
                }
                else
                {
                    m_fromField[i]->m_fld->Import(infiles[i],m_fromField[i]->m_fielddef,
                                                  m_fromField[i]->m_data,
                                                  LibUtilities::NullFieldMetaDataMap);
                }
                
                nfields    = m_fromField[i]->m_fielddef[0]->m_fields.size();
                int NumHomogeneousDir = m_fromField[i]->m_fielddef[0]->m_numHomogeneousDir;
   
                // Set up Expansion information to use mode order from field
                m_fromField[i]->m_graph->SetExpansions(m_fromField[i]->m_fielddef);
                
                //Set up expansions, and extract data. 
                m_fromField[i]->m_exp.resize(nfields);
                m_fromField[i]->m_exp[0] = m_fromField[i]->SetUpFirstExpList(NumHomogeneousDir,true);
                
                for(j = 1; j < nfields; ++j)
                {
                    m_fromField[i]->m_exp[j] = m_f->AppendExpList(NumHomogeneousDir);
                }
                
                for (j = 0; j < nfields; ++j)
                {
                    for (int k = 0; k < m_fromField[i]->m_data.size(); ++k)
                    {
                        m_fromField[i]->m_exp[j]->ExtractDataToCoeffs(
                            m_fromField[i]->m_fielddef[k],
                            m_fromField[i]->m_data[k],
                            m_fromField[i]->m_fielddef[k]->m_fields[j],
                            m_fromField[i]->m_exp[j]->UpdateCoeffs());
                    }
                    m_fromField[i]->m_exp[j]->BwdTrans(m_fromField[i]->m_exp[j]->GetCoeffs(),
                                                       m_fromField[i]->m_exp[j]->UpdatePhys());
                }
            }            


