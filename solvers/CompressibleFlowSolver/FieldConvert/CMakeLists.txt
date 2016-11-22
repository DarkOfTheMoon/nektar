SET(FieldConvertSources
        ../../../utilities/FieldConvert/FieldConvert.cpp
        ./ProcessModules/ProcessNonConservativeVariables.cpp
    )       

ADD_SOLVER_EXECUTABLE(FieldConvertCFS solvers ${FieldConvertSources})
TARGET_LINK_LIBRARIES(FieldConvertCFS FieldUtils)

# ADD_NEKTAR_TEST(chan3D_pts)
