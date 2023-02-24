add_executable(main_tvlqr
    riccati.h riccati.c main_tvlqr.c)

add_executable(main_lqr
    riccati.h riccati.c main_lqr.c)

target_link_libraries(main_tvlqr PRIVATE slap::slap)
target_link_libraries(main_lqr PRIVATE slap::slap)