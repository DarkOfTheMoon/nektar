SET(ViscoelasticSolverSource 
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/Viscoelastic.cpp
       ./EquationSystems/LogConformation.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./AdvectionTerms/AdjointAdvection.cpp
       ./AdvectionTerms/AdvectionTerm.cpp
       ./AdvectionTerms/LinearisedAdvection.cpp
       ./AdvectionTerms/NavierStokesAdvection.cpp
       ./AdvectionTerms/SkewSymmetricAdvection.cpp
       ./AdvectionTerms/NoAdvection.cpp
       ./ViscoelasticSolver.cpp 
       )

IF (NEKTAR_USE_ARPACK)
    SET(ViscoelasticSolverSource ${ViscoelasticSolverSource}
        ../Auxiliary/DriverArpack.cpp)
ENDIF (NEKTAR_USE_ARPACK)

ADD_SOLVER_EXECUTABLE(ViscoelasticSolver solvers
			${ViscoelasticSolverSource})
                        
