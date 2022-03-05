# Be careful with reordering. The slides include the heatequation part with
# line numbers.
dune_add_generated_executable(TARGET heatequation
                              UFLFILE heatequation.ufl
                              INIFILE heatequation.ini
                              SOURCE heatequation_driver.cc
                              )

dune_symlink_to_source_files(FILES heatequation.ini)


dune_add_generated_executable(TARGET poisson
                              UFLFILE poisson.ufl
                              INIFILE poisson.ini
                              )

dune_add_generated_executable(TARGET nonlinear_poisson
                              UFLFILE nonlinear_poisson.ufl
                              INIFILE nonlinear_poisson.ini
                              )

dune_add_generated_executable(TARGET nonlinear_poisson_nitschebc
                              UFLFILE nonlinear_poisson_nitschebc.ufl
                              INIFILE nonlinear_poisson_nitschebc.ini
                              )

dune_add_generated_executable(TARGET wave_equation
                              UFLFILE wave_equation.ufl
                              INIFILE wave_equation.ini
                              )



dune_symlink_to_source_files(FILES poisson.ini nonlinear_poisson.ini nonlinear_poisson_nitschebc.ini wave_equation.ini)
