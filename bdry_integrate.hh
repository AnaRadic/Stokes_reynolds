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
#include "driver.hh"

template<template<typename ctype, int dim> class Function, class Grid>
double bdry_integration(Grid& grid, int p) {

   std::vector<std::string> rec0fLines;
   std::string oneLine;

   ifstream file0("obstacle.msh");
   if(file0.is_open())
   {
       while(getLine(file0,oneLine))
       {
           rec0fLines.push_back(oneLine);
       }
       file0.close();
   }
   else std::cout << "File not opened";

   //test
   for(auto it=rec0fLines.begin();it!=rec0fLines.end();it++)
   {
       std::cout <<*it << std::endl;
   }
    const int dim = Grid::dimension;
    typedef typename Grid::ctype ctype;

    // Vektorsko polje koje integriramo
    Function<ctype, dim> fptr;

    // leaf GridView
    auto gridView = grid.leafGridView();

    auto endit = gridView.template end<0>();
    auto it = gridView.template begin<0>();

    double integral = 0.0;
    // Petlja po svim elementima
    for (; it != endit; ++it) {
        // Intersection  iteratori
        auto isit_end = gridView.iend(*it);
        auto isit = gridView.ibegin(*it);

        double elem_integral = 0.0;
        for (; isit != isit_end; ++isit) { // petlja po svim stranicama elementa

            if (isit->boundary())  // Jesmo li na granici domene?
            {
                const auto igeo = isit->geometry();
                const auto gt = igeo.type();
                //std::cout << igeo.center() << std::endl;

                auto outerNormal = isit->centerUnitOuterNormal();
                // std::cout << outerNormal << std::endl;

                // Zatražimo kvadraturnu formulu (obavezno koristiti referencu)
                // Formula mora biti dimenzije dim -1!
                const auto& rule = Dune::QuadratureRules<ctype, dim - 1>::rule(gt, p);

                //if (rule.order() < p)
                  //  DUNE_THROW(Dune::Exception, "order not available");

                double result = 0.0;
                // Petlja po svim integracijskim toèkama
                auto iq = rule.begin();
                for (; iq != rule.end(); ++iq) {
                    // Na iteratoru i zovemo metode position() i weight()
                    auto fval = fptr(igeo.global(iq->position()));
                    double weight = iq->weight();
                    // | det (grad g) | daje Geometry objekt
                    double detjac = igeo.integrationElement(iq->position());
                    result += (pressure * outerNormal) * weight * detjac;
                }  // kraj petlje po svim integracijskim toèkama

                //  std::cout << result << std::endl;
                elem_integral += result;

            } // kraj iif(isit->boundary()
        }
        integral += elem_integral;
    } // kraj petlje po svim elementima

    return integral;
}
