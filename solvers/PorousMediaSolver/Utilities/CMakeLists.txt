IF (NEKTAR_USE_VTK)
    SET(ImageToFldSources ImageToFld.cpp)
    ADD_SOLVER_EXECUTABLE(ImageToFld solvers ${ImageToFldSources})
ENDIF (NEKTAR_USE_VTK)
