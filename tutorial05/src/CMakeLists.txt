add_executable(tutorial05 tutorial05.cc)

dune_symlink_to_source_files(FILES tutorial05.ini ldomain.msh ldomain.geo)
