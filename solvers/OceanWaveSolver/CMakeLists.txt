CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_OCEANWAVE
    "Build the ocean wave solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NEKTAR_SOLVER_OCEANWAVE)
    SET(OceanWaveSolverSource
       ./OceanWaveSolver.cpp
       ./EquationSystems/OceanWaveSystem.cpp
       ./EquationSystems/OceanWaveSystem3D.cpp)

    ADD_SOLVER_EXECUTABLE(OceanWaveSolver solvers
        ${OceanWaveSolverSource})
ENDIF()
