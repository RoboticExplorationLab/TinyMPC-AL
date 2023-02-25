add_executable(main_tvlqr
    riccati.h riccati.c util.h main_tvlqr.c)

add_executable(main_lqr
    riccati.h riccati.c util.h main_lqr.c)

add_executable(main_lqr_track
    riccati.h riccati.c util.h main_lqr_track.c)
    
target_link_libraries(main_tvlqr PRIVATE slap::slap)
target_link_libraries(main_lqr PRIVATE slap::slap)
target_link_libraries(main_lqr_track PRIVATE slap::slap)