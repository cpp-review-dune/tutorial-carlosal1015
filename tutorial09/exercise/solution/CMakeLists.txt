dune_add_generated_executable(TARGET solution_navier_stokes
                              UFLFILE navier_stokes.ufl
                              INIFILE solution_navier_stokes.ini
                              )

dune_add_generated_executable(TARGET solution_nonlinear_poisson_dg
                              UFLFILE nonlinear_poisson_dg.ufl
                              INIFILE solution_nonlinear_poisson_dg.ini
                              )

dune_symlink_to_source_files(FILES solution_navier_stokes.ini solution_nonlinear_poisson_dg.ini cylinder2dmesh1.msh cylinder2dmesh2.msh cylinder2dmesh3.msh cylinder2dmesh4.msh)
