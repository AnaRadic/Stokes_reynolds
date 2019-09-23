#ifndef CG_STOKES_INITIAL_HH
#define CG_STOKES_INITIAL_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/localoperator/stokesparameter.hh>
#include <dune/pdelab/common/function.hh>

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================


// Klasa koja određuje tip granice
class BCTypeParam
{
private:
    double time;
public:
    // Ova klasa daje indekse: DoNothing i VelocityDirichlet i StressNeumann
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits

  {
    typedef BC::Type RangeType;
    using BoundaryCondition = Dune::PDELab::StokesBoundaryCondition;//?
  };

  template<typename I>
  inline void evaluate (const I & intersection,
                        const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord,
                        BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::coorddimension>
        xg = intersection.geometry().global( coord );
    if( xg[0] > 10.0-1e-6 )
      y = BC::DoNothing;
    else
      y = BC::VelocityDirichlet;
  }
  template <typename T>
   void setTime(T t){
     time = t;
   }
};


// Dirichletov rubni uvjet za brzinu
template<typename GV, typename RF, int dim>
class Velocity :
  public Dune::PDELab::AnalyticGridFunctionBase<
                             Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
                             Velocity<GV,RF,dim>
                                               >
{
private:
  RF time;
  RF intensity;//?

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Velocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  Velocity(const GV & gv, double intensity_ = 1.0) : BaseT(gv), intensity(intensity_) {
    time = 0.0;
  }

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    RF r = 0;

    y[1] = 0;  // vertikalna brzina je svugdje nula


    if(x[0] < 0.0+1e-6){
      y[0] = 1;
    }
    else{
      y[0]=0;
    }
  }

  template <typename T>
  void setTime(T t){
    time = t;
  }

};




// Vektorska funkcija jednaka nuli. Ovdje je koristimo  za tlak (dim_range=1)
// i za Neumannov rubni uvjet (dim_range=dim) koji nije prisutan pa ga stavljamo na nulu.
// Rubni uvjet za tlak jednako tako ne postoji i stoga vrijednost tlaka stavljamo na nulu.
template<typename GV, typename RF, std::size_t dim_range>
class ZeroFunction :
  public Dune::PDELab::AnalyticGridFunctionBase<
                            Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim_range>,
                            ZeroFunction<GV,RF,dim_range>
                                               >,
  public Dune::PDELab::InstationaryFunctionDefaults
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim_range> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ZeroFunction> BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  ZeroFunction(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    y=0;
  }
};



#endif // BC_EXTENSION_HH
