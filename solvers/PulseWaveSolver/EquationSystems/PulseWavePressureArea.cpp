///////////////////////////////////////////////////////////////////////////////
//
// File PulseWavePressureArea.cpp
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
// Description: PulseWavePressureArea definition
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/PulseWavePressureArea.h>
#include <loki/Singleton.h>

namespace Nektar
{

    /**
     * @class PulseWavePressureArea
     *
     */
  PulseWavePressureArea::PulseWavePressureArea(
      Array<OneD, MultiRegions::ExpListSharedPtr> &pVessel, 
      const LibUtilities::SessionReaderSharedPtr &pSession)
    : m_vessels(pVessel),
      m_session(pSession)
    {
    }

    PulseWavePressureArea::~PulseWavePressureArea()
    {
    }

        
    /**
     *
     */
    PressureAreaFactory& GetPressureAreaFactory()
    {
        typedef Loki::SingletonHolder<PressureAreaFactory,
                                      Loki::CreateUsingNew,
                                      Loki::NoDestroy > Type;
        return Type::Instance();
    }

    void PulseWavePressureArea::EvaluateFunction(
        Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
        LibUtilities::SessionReaderSharedPtr              pSession,
        std::string                                       pFieldName,
        Array<OneD, NekDouble>&                           pArray,
        const std::string&                                pFunctionName,
        NekDouble                                         pTime,
        const int domain)
    {
        ASSERTL0(pSession->DefinesFunction(pFunctionName),
                 "Function '" + pFunctionName + "' does not exist.");

        unsigned int nq = pFields[0]->GetNpoints();
        if (pArray.num_elements() != nq)
        {
            pArray = Array<OneD, NekDouble> (nq);
        }

        LibUtilities::FunctionType vType;
        vType = pSession->GetFunctionType(pFunctionName, pFieldName);
        if (vType == LibUtilities::eFunctionTypeExpression)
        {
            Array<OneD, NekDouble> x0(nq);
            Array<OneD, NekDouble> x1(nq);
            Array<OneD, NekDouble> x2(nq);

            pFields[0]->GetCoords(x0, x1, x2);
            LibUtilities::EquationSharedPtr ffunc =
                pSession->GetFunction(pFunctionName, pFieldName,domain);

            ffunc->Evaluate(x0, x1, x2, pTime, pArray);
        }
        else if (vType == LibUtilities::eFunctionTypeFile)
        {
            std::string filename = pSession->GetFunctionFilename(
                pFunctionName,
                pFieldName,
                domain);

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;
            Array<OneD, NekDouble> vCoeffs(pFields[0]->GetNcoeffs());
            Vmath::Zero(vCoeffs.num_elements(), vCoeffs, 1);

            int numexp = pFields[0]->GetExpSize();
            Array<OneD,int> ElementGIDs(numexp);

            // Define list of global element ids
            for(int i = 0; i < numexp; ++i)
            {
                ElementGIDs[i] = pFields[0]->GetExp(i)->GetGeom()->GetGlobalID();
            }

            LibUtilities::FieldIOSharedPtr fld =
                MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(m_session->GetComm());

            fld->Import(filename,FieldDef,FieldData,
                          LibUtilities::NullFieldMetaDataMap,
                          ElementGIDs);
            int idx = -1;

            // Loop over all the expansions
            for (int i = 0; i < FieldDef.size(); ++i)
            {
                // Find the index of the required field in the
                // expansion segment
                for(int j = 0; j < FieldDef[i]->m_fields.size(); ++j)
                {
                    if (FieldDef[i]->m_fields[j] == pFieldName)
                    {
                        idx = j;
                    }
                }

                if (idx >= 0)
                {
                    pFields[0]->ExtractDataToCoeffs(
                        FieldDef[i], FieldData[i],
                        FieldDef[i]->m_fields[idx], vCoeffs);
                }
                else
                {
                    cout << "Field " + pFieldName + " not found." << endl;
                }
            }


            pFields[0]->BwdTrans_IterPerExp(vCoeffs, pArray);

        }
    }

}
