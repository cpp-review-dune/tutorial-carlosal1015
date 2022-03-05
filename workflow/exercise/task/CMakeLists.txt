add_executable(grid-exercise2 grid-exercise2.cc)

dune_symlink_to_source_files(FILES complex.msh
                                   ldomain.msh
                                   unitcube.msh
                                   unitsquare.msh
                                   2dexample.vtu)
