SET(NEKTAR_SOLVER_POROUSMEDIA ON CACHE BOOL
    "Build the porous media solver.")

IF( NEKTAR_SOLVER_POROUSMEDIA )
    SET(PorousMediaSolverSource 
       ./EquationSystems/PorousMediaSplittingScheme.cpp
       ./EquationSystems/PorousMedia.cpp
       ./PorousMediaSolver.cpp 
       )

    ADD_SOLVER_EXECUTABLE(PorousMediaSolver solvers 
			${PorousMediaSolverSource})

ENDIF( NEKTAR_SOLVER_POROUSMEDIA )
