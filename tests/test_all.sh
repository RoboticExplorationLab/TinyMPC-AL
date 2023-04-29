#!/bin/bash
DIR="./tests/src";
echo "Run all tests at $DIR"
$DIR/lqrdata_test
$DIR/lqr_ineq_utils_test
$DIR/forward_pass_test
$DIR/cost_test
$DIR/lqr_lti_test
$DIR/lqr_lti_track_test
$DIR/al_lqr_lti_test
# $DIR/lqr_ltv_test
$DIR/lqr_ltv_track_test
$DIR/al_lqr_ltv_test
