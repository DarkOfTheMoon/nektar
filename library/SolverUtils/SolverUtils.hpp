namespace Nektar 
{
    namespace SolverUtils
    {
        enum EvolutionOperatorType
        {
            eNonlinear,
            eDirect,
            eAdjoint,
            eTransientGrowth,
            eTransientAdjoint,
            eSkewSymmetric,
            eAdaptiveSFD
        };
    }
}
