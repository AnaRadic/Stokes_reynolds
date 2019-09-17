#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#include "driver.hh"

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>

template <int dim>
using Point = Dune::FieldVector<double, dim>;


template <typename GV, typename PARAMS>
double bdry_integration(const GV &gv, double dt, double tfin,PARAMS & parameters)

{

  double time = 0.0;
  const int dim = 2;
  double integral = 0.0;
  
  // Petlja po svim elementima
  while (time <= tfin-1e-8) {
     
     for (auto const &element : elements(gv)) {
    
        double elem_integral = 0.0;
        // petlja po svim stranicama elementa
       for (auto const &side : intersections(gv, element)) {
     //  auto xg = intersections.geometry().global( coord ); //dohvaćamo koordinate

       if (side.boundary()/*&& xg[0]>0.0 )*/) // Jesmo li na granici domene?
       {
        const auto sidegeo = side.geometry();
        auto outerNormal = side.centerUnitOuterNormal();
        // Zatrazimo kvadraturnu formulu na stranici (dimenzija dim -1)
        const auto &rule =
            Dune::QuadratureRules<double, dim - 1>::rule(sidegeo.type(), tfin);
        if (rule.order() < tfin)
             std::cerr << "Integracijska formula reda " << tfin << " nije dostupna.\n";
        double result = 0.0;
        // Petlja po svim integracijskim toÄŤkama
        for (auto const &qpoint : rule) {
          double weight = qpoint.weight();
          // | det (grad g) | daje Geometry objekt
          double detjac = sidegeo.integrationElement(qpoint.position());
          Point<dim> pressure (0.0);
          result += (pressure * outerNormal) * weight * detjac;
       }

        elem_integral += result;
      }
    } // kraj petlje po svim stranicama
    integral += elem_integral;
  }
   // kraj petlje po svim elementima
    time+=dt;
}
  return integral;

}


