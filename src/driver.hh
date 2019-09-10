#ifndef DRIVER_HH
#define	DRIVER_HH

#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
//#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/seqistlsolverbackend.hh>
#include <dune/pdelab/backend/simple/descriptors.hh>


#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>//dodano
//#include <dune/pdelab/backend/istlmatrixbackend.hh>
//#include <dune/pdelab/backend/istlsolverbackend.hh>
//#include <dune/pdelab/localoperator/cg_stokes.hh>

#include <dune/pdelab/localoperator/taylorhoodnavierstokes.hh>
#include <dune/pdelab/newton/newton.hh>

template<typename GV, typename IF, typename PARAMS>
void driver(const GV& gv, std::string filename, PARAMS & parameters, IF & bdry_solution,double dt, double tfin)
{
    using namespace Dune::PDELab;  // skrati imena

    static const unsigned int dim = GV::dimension;
    typedef double RF;

    Dune::Timer timer;
        std::cout << "=== Initialize:" << timer.elapsed() << std::endl;
        timer.reset();

    // Konstrukcija prostora konačnih elemenata
        typedef typename GV::Grid::ctype DF;
    const int k = 2;
    const int q = 2 * k;  // preciznost integracijske formule lokalnog operatora

    // Taylor-Hoodovi elementi -- P2 za brzinu P1 za tlak
    typedef PkLocalFiniteElementMap<GV, double, double, k>      V_FEM;  // komponenta brzine
    typedef PkLocalFiniteElementMap<GV, double, double, k - 1 > P_FEM;  // tlak
    V_FEM vFem(gv);
    P_FEM pFem(gv);

    typedef Dune::PDELab::ConformingDirichletConstraints CDC;
 //   typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none, 1> VB; SMETA MU DIMENZIJA
    typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none> VB;

    //using CDC = ConformingDirichletConstraints;//mozda maknuti
    //using VB = ISTL::VectorBackend<>;//mozda maknuti
    // Ova klasa direktno konstruira vektorske elemente u R^dim.
    // Prostor mrežnih funkcija za brzinu (vektorski):  V_h

    typedef Dune::PDELab::VectorGridFunctionSpace<GV, V_FEM, dim, VB, VB, CDC> V_GFS;//ili PGFS_V_GFS
    //using V_GFS = VectorGridFunctionSpace<GV, V_FEM, dim, VB, VB, CDC>;//MOZDA MAKNUTI
    V_GFS VGfs(gv, vFem);
    VGfs.name("velocity");


    // Prostor mrežnih funkcija za tlak (skalarni): W_h

    typedef Dune::PDELab::GridFunctionSpace<GV, P_FEM, CDC, VB> P_GFS;
    //using P_GFS = GridFunctionSpace<GV, P_FEM, CDC, VB>;//MOZDA MAKNUTI
    P_GFS pGfs(gv, pFem);
    pGfs.name("pressure");


    // Prostor V_h x W_h
    // LexicographicOrderingTag daje poredak varijabli: v1,v2,p
    typedef Dune::PDELab::CompositeGridFunctionSpace
                <VB, Dune::PDELab::LexicographicOrderingTag, V_GFS, P_GFS> GFS;
    //using GFS = CompositeGridFunctionSpace<VB, LexicographicOrderingTag, V_GFS, P_GFS> ;//MOZDA MAKNUTI
    GFS gfs(VGfs, pGfs);

    // Primjena Dirichletovih ograničenja
    typedef typename GFS::template ConstraintsContainer<RF>::Type C;
    //using C = typename GFS::template ConstraintsContainer<double>::Type;//MOZDA MAKNUTI
    C cg;
    cg.clear();


    // create Taylor-Hood constraints from boundary-type
        typedef Dune::PDELab::StokesVelocityDirichletConstraints<PARAMS>                                  ScalarVelocityConstraints;
        typedef Dune::PDELab::PowerConstraintsParameters<ScalarVelocityConstraints, dim>               VelocityConstraints;
        typedef Dune::PDELab::StokesPressureDirichletConstraints<PARAMS>                                  PressureConstraints;
        typedef Dune::PDELab::CompositeConstraintsParameters<VelocityConstraints, PressureConstraints> Constraints;

    // Određivanje Dirichletove granice. Ovdje se koriste pomoćne klase koje određuju
    // Dirichletovu granicu za svaku komponentu vektorske funkcije (v1,v2,p).
    //using ScalarVelConstraints = StokesVelocityDirichletConstraints<PARAMS>;//MOZDA MAKNUTI
    //using VelocityConstraints = PowerConstraintsParameters<ScalarVelConstraints, dim>;//MOZDA MAKNUTI
    //using PressureConstraints = StokesPressureDirichletConstraints<PARAMS>;
    //using Constraints = CompositeConstraintsParameters<VelocityConstraints, PressureConstraints>;

    ScalarVelocityConstraints scalarvelocity_constraints(parameters);
    VelocityConstraints  velocity_constraints(scalarvelocity_constraints);
    PressureConstraints  pressure_constraints(parameters);
    Constraints          bconst(velocity_constraints, pressure_constraints);

    // Odredi Dirichletova ograničenja
    Dune::PDELab::constraints(bconst, gfs, cg);  // MIJEŠANJE IMENA
//    constraints(bconst, gfs, cg);  --

    //Dune::PDELab::NavierStokesDefaultParameters< GV, RF, F, B, V, J, navier, tensor >//mozda


    const bool navier = true;
    const int qorder = 2 * k;
    typedef Dune::PDELab::TaylorHoodNavierStokes<PARAMS> SLOP;
    SLOP slop(parameters);
    typedef Dune::PDELab::NavierStokesMass<PARAMS> TLOP;
    TLOP tlop(parameters);



 //   Dune::PDELab::Simple::MatrixBackend< Container >::Pattern< Matrix, GFSV, GFSU > MBE;

 //   typedef VB::MatrixBackend MBE;
    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    MBE mbe(5);
    typedef Dune::PDELab::GridOperator<GFS,GFS,SLOP,MBE,RF,RF,RF,C,C> GO1; // CC -> C
    GO1 go1(gfs,cg,gfs,cg,slop,mbe);  // cc-> cg
    typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,RF,RF,RF,C,C> GO2;
    GO2 go2(gfs,cg,gfs,cg,tlop,mbe);
    typedef Dune::PDELab::OneStepGridOperator<GO1,GO2> GO;
    GO go(go1,go2);

    // Prostorni lokalni operator - definiran u Dune::PDELab-u.
    using LOP = TaylorHoodNavierStokes<PARAMS>;
    LOP lop(parameters, q);

    //Mrežni operator
//    using MBE = ISTL::BCRSMatrixBackend<>;
//    MBE mbe(5); // maksimalan broj ne-nul elemenat u retku (samo pretpostavka)
//    using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, double, double, double, C, C>;
//    GO go(gfs, cg, gfs, cg, lop, mbe);

    // Vektor koeficijenata i interpolacija rubnog uvjeta
    using U = typename GO::Traits::Domain;
    U x0(gfs);
    x0 = 0.0;
    Dune::PDELab::interpolate(bdry_solution, gfs, x0);
    std::cout << "=== Finished interpolation:" << timer.elapsed() << std::endl;
    timer.reset();

    interpolate(bdry_solution, gfs, x0);

    // Postavi sve stupnjeve slobode koji nisu Dirichletovi na nulu.
    //set_shifted_dofs(cg, 0.0, x0);



    // Set non constrained dofs to zero
    Dune::PDELab::set_shifted_dofs(cg, 0.0, x0);


    // Linear solver

    // Linear solver. Koristimo superLU, direktan solver za rijetke matrice
    // koji može riješiti indefinitne sustave.
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LinearSolver;
    LinearSolver ls(false);

    //  Newtonow solver
    typedef Dune::PDELab::Newton<GO,LinearSolver,U> PDESOLVER;
    PDESOLVER newton(go,x0,ls);//dodaLA X0
    newton.setReassembleThreshold(0.0);
    newton.setVerbosityLevel(2);
    newton.setMaxIterations(25);
    newton.setLineSearchMaxIterations(30);
    //  Izbor metode konačnih diferencija
    Dune::PDELab::OneStepThetaParameter<RF> method(1); // implicit, theta = 1
    //  Dune::PDELab::Alexander2Parameter<RF> method;               // koeficijenti metode
    Dune::PDELab::OneStepMethod<RF,GO,PDESOLVER,U,U> osm(method,go,newton);
    osm.setVerbosityLevel(2);

//Dune::PDELab::NavierStokesDefaultParameters< GV, RF, F, B, V, J, navier, tensor >/MOZDA
    //using LS = ISTLBackend_SEQ_SuperLU;
    //LS ls(false);



    // Riješi sustav.
    //Newton<GO, LS, U> newton(go, x0, ls);
    //newton.setReassembleThreshold(0.0);
    //newton.setVerbosityLevel(2);
    //newton.setMaxIterations(25);
    //newton.setLineSearchMaxIterations(30);
    //newton.apply();

    // Izračunaj rezidual i ispiši njegovu normu.
//    U r(gfs);
//    r = 0.;
//    go.residual(x0, r);
//    std::cout << "Final Residual: " << r.two_norm() << std::endl;

    // Ispis rješenja. NOVA VERZIJA KORISTI SE SEQUENCE WRITER
    typedef Dune::SubsamplingVTKWriter<GV>  VTKW;
    VTKW vtkwriter(gv, Dune::RefinementIntervals{k});
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x0);
    Dune::VTKSequenceWriter<GV> writer(std::make_shared<VTKW>(vtkwriter), "out");

/*
    Dune::PDELab::FilenameHelper fn(filename);
    {
        Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, 2);
        Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x0);
        vtkwriter.write(fn.getName(), Dune::VTK::ascii);
        fn.increment();
    }
*/

    double time = 0.0;
    writer.write(time);
    U x1(gfs,0.0);        // sljedeći vremenski sloj -jednokoračna metoda x1=x^{n+1}, x0=x^{n}
    while (time < tfin-1e-8) {
        // postavi novo vrijeme u BC klasu
        bdry_solution.setTime(time+dt);
        cg.clear();
        Dune::PDELab::constraints(bconst,gfs,cg);

        osm.apply(time, dt, x0, bdry_solution, x1);                          // riješi sustav

        U r(gfs);
        r = 0.;

        go.residual(x1, r);
        std::cout << "Final Residual: " << r.two_norm() << std::endl;
        // graphics  STARI KOD!
//        {
//            Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, 2);
//            Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x0);
//            vtkwriter.write(fn.getName(), Dune::VTK::ascii);
//            fn.increment();
//        }
        x0 = x1;                                             // pripremi sljedeći vremenski korak
        time += dt;
        writer.write(time);

    }


}

#endif	/* DRIVER_HH */
