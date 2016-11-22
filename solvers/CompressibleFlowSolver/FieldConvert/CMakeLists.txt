SET(FieldConvertSources
        ../../../utilities/FieldConvert/FieldConvert.cpp
        ./ProcessModules/ProcessNonConservativeVariables.cpp
    )       

ADD_SOLVER_EXECUTABLE(FieldConvert solvers ${FieldConvertSources})
TARGET_LINK_LIBRARIES(FieldConvert FieldUtils)

# ADD_NEKTAR_TEST(chan3D_pts)
