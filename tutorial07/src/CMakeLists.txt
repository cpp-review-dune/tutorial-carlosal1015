add_executable("tutorial07-maxwell" tutorial07-maxwell.cc)
add_executable("tutorial07-linearacoustics" tutorial07-linearacoustics.cc)
add_executable("tutorial07-shallowwater" tutorial07-shallowwater.cc)

dune_symlink_to_source_files(FILES splitdomain.msh)
dune_symlink_to_source_files(FILES tutorial07-maxwell.ini)
dune_symlink_to_source_files(FILES tutorial07-linearacoustics.ini)
dune_symlink_to_source_files(FILES tutorial07-shallowwater.ini)
