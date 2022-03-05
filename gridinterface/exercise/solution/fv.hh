#ifndef FV_HH
#define FV_HH

#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/grid/common/mcmgmapper.hh>


// template parameters: GV = gridview type, E = equation parameters, ScalarField = solution vector
template<typename GV, typename E, typename ScalarField>
class FiniteVolume
{

  // mapper for elements (codim=0)
  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

  GV grid_view;
  const E& equation;
  Mapper mapper;
  ScalarField update;

public:
  FiniteVolume(const GV& grid_view_, const E& equation_)
    : grid_view(grid_view_)
    , equation(equation_)
    , mapper(grid_view_, Dune::mcmgElementLayout())
    , update(grid_view.size(0))
  {}

  void update_concentration(ScalarField& concentration, double dt)
  {

    // initialize update vector with 0
    std::fill(update.begin(),update.end(),0.0);

    for (const auto& cell : elements(grid_view))
      {

        auto cell_index = mapper.index(cell);

        auto geometry = cell.geometry();

        // cell volume
        auto volume = geometry.volume();

        // convection term -- boundary integral
        for (const auto& intersection : intersections(grid_view,cell))
          {

            auto face_geometry = intersection.geometry();

            // get the global coordinates of the face's center
            auto center = face_geometry.center();

            // get the outer normal vector at the face center
            auto normal = intersection.centerUnitOuterNormal();

            // get the velocity vector at the face center
            auto u = equation.u(center);

            // compute face integral
            auto flux = (u * normal) * face_geometry.volume() / volume;

            // handle interior face
            if (intersection.neighbor())
              {
                auto outside_index = mapper.index(intersection.outside());

                if (flux > 0.0)           // upwind decision
                  // outflow
                  update[cell_index] -= concentration[cell_index]*flux;
                else
                  // inflow
                  update[cell_index] -= concentration[outside_index]*flux;
              }

            // handle boundary face
            if (intersection.boundary())
              {
                if (flux > 0.0)           // Inflow or outflow?
                  // outflow
                  update[cell_index] -= concentration[cell_index]*flux;
                else
                  // inflow
                  update[cell_index] -= equation.c_in(center)*flux;
              }
          }
      }

    // update solution
    for (std::size_t i = 0; i < concentration.size(); ++i)
      concentration[i] += update[i] * dt;
  }

};

#endif
