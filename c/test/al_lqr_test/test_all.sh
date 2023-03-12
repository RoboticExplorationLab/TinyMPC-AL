#!/bin/bash
DIR="$( dirname -- "${BASH_SOURCE[0]}"; )";
echo "Run all tests at $DIR"
build/$DIR/lqrdata_test
build/$DIR/lqr_ineq_utils_test
build/$DIR/forward_pass_test
build/$DIR/cost_test
build/$DIR/lqr_lti_test
build/$DIR/lqr_lti_track_test
build/$DIR/al_lqr_lti_test
build/$DIR/lqr_ltv_test
build/$DIR/lqr_ltv_track_test
build/$DIR/al_lqr_ltv_test

