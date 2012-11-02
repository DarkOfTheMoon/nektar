/*
 * GlobalLinSysXxt.h
 *
 *  Created on: 19 Oct 2012
 *      Author: cc
 */

#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSXXT_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSXXT_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

namespace Xxt
{
    struct crs_data;
}

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations

        //class AssemblyMapDG;
        class ExpList;

        class GlobalLinSysXxt : public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysXxt(
                    const GlobalLinSysKey &pKey,
                    const boost::weak_ptr<ExpList> &pExp,
                    const boost::shared_ptr<AssemblyMap>
                                                            &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysXxt();

        protected:
            struct Xxt::crs_data*       m_crsData;
            Array<OneD, unsigned int>   m_Ai;
            Array<OneD, unsigned int>   m_Aj;
            Array<OneD, double>         m_Ar;

            Array<OneD, NekDouble>      m_locToGloSignMult;

            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            virtual void v_Solve(
                    const Array<OneD, const NekDouble> &in,
                    Array<OneD, NekDouble> &out,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray);

            /// Solve the linear system for given input and output vectors.
            virtual void v_SolveLinearSystem(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const int pNumDir = 0);

            void GlobalToLocalNoSign(const Array<OneD, const NekDouble> &global,
                                           Array<OneD, NekDouble> &local,
                                   const boost::shared_ptr<AssemblyMap>
                                                          &pLocToGloMap);

            void LocalToGlobalNoSign(const Array<OneD, const NekDouble> &local,
                                           Array<OneD, NekDouble> &global,
                                   const boost::shared_ptr<AssemblyMap>
                                                          &pLocToGloMap);

        };
    }
}
#endif /* GLOBALLINSYSXXT_H_ */
