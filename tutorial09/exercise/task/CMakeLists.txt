dune_add_generated_executable(TARGET nonlinear_poisson_dg
                              UFLFILE nonlinear_poisson_dg.ufl
                              INIFILE nonlinear_poisson_dg.ini
                              )

dune_add_generated_executable(TARGET navier_stokes
                              UFLFILE navier_stokes.ufl
                              INIFILE navier_stokes.ini
                              )

dune_symlink_to_source_files(FILES nonlinear_poisson_dg.ini navier_stokes.ini cylinder2dmesh1.msh cylinder2dmesh2.msh cylinder2dmesh3.msh cylinder2dmesh4.msh)
