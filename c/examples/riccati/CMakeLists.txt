add_executable(riccati_example
    riccati.h riccati.c riccati_main.c)

target_link_libraries(riccati_example PRIVATE slap::slap)