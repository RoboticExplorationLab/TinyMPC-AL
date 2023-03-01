add_executable(main_tvlqr
    riccati.h riccati.c util.h main_tvlqr.c)

add_executable(main_lqr
    riccati.h riccati.c util.h main_lqr.c)

add_executable(main_lqr_track
    riccati.h riccati.c util.h main_lqr_track.c)

add_executable(planar_quad_lqr
    riccati.h riccati.c util.h quadrotor.h planar_quad_lqr.c)

add_executable(bicycle_tvlqr
    riccati.h riccati.c util.h bicycle.h bicycle_tvlqr.c)

target_link_libraries(main_tvlqr PRIVATE slap::slap)
target_link_libraries(main_lqr PRIVATE slap::slap)
target_link_libraries(main_lqr_track PRIVATE slap::slap)
target_link_libraries(planar_quad_lqr PRIVATE slap::slap)
target_link_libraries(bicycle_tvlqr PRIVATE slap::slap)
