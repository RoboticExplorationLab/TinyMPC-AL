# function add_al_lqr_test(name)
#
# Adds a new executable test called <name>_test.
# Assumes thes source code is in a file called <name>_test.
function(add_my_test name)
  set(TEST_NAME ${name}_test)
  add_executable(${TEST_NAME}
    ${TEST_NAME}.cpp

    test_utils.h
    test_utils.c

    models/bicycle_5d.h  # we don't want to include for all
    models/bicycle_5d.c

    models/planar_quadrotor.h   # we don't want to include for all
    models/planar_quadrotor.c
    )
  target_link_libraries(${TEST_NAME}
    PRIVATE
    tinympc::tinympc
    slap::slap
    gtest::gtest
    )
  add_test(NAME ${TEST_NAME}
    COMMAND ${TEST_NAME}
    )
endfunction()

# Add tests
add_my_test(model_lti)
add_my_test(model_ltv)
# add_my_test(lqrdata)
# add_my_test(lqr_ineq_utils)
# add_my_test(forward_pass)
# add_my_test(cost)
# add_my_test(lqr_lti)
# add_my_test(lqr_lti_track)
# add_my_test(al_lqr_lti)
# add_my_test(lqr_ltv)
# add_my_test(lqr_ltv_track)
# add_my_test(al_lqr_ltv)