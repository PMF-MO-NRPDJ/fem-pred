#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <array>
#include <cmath>
#include <string>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh> 

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/common/fvector.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

// Egzaktno rješenje kao funkcija globalne koordinate.
template<int dim>
double exact(Dune::FieldVector<double,dim> const & x){
      return 1.0 + 3*(x*x);
}

// Slobodni koeficijent
template<int dim>
double a0(Dune::FieldVector<double,dim> const & x){
      return 1.0 + x*x;
}

// Desna strana jednadžbe.
// Aproksimiramo desnu stranu pomoću konačnih diferencija.
template <int dim>
double f(Dune::FieldVector<double, dim> x) {
  std::array<Dune::FieldVector<double, dim>, dim> xp;
  std::array<Dune::FieldVector<double, dim>, dim> xm;
  std::array<double, dim> cp;
  std::array<double, dim> cm;
  const double h = 1E-6; // Korak diferenciranja
  const double h2 = h * h;
  for (unsigned int i = 0; i < dim; ++i) // i je koordinatni smjer
  {
    xp[i] = x;
    xm[i] = x;
    xp[i][i] += h;
    xm[i][i] -= h;
    cp[i] = exact(xp[i]);
    cm[i] = exact(xm[i]);
  }
  auto c0 = exact(x);
  // operator - Laplace aproksimiran konačnim diferencijama u točki x.
  double result = 0.0;
  for (int i = 0; i < dim; i++)
    result -= (cp[i] - 2 * c0 + cm[i]) / h2;
  // dodaj slobodni član
  result += a0(x) * c0;

  return result;
}

// GV = GridView
// Vec = vektor (=BlockVector<...>)
// Izračunaj egzaktno rješenje u vrhovima mreže i vrati izračunate vrijednosti kroz coeff.
template<class GV, class Vec>
void vertexdata (const GV& gridView, Vec & coeff)
{
  const int dim = GV::dimensionworld;
  auto const & idxSet = gridView.indexSet();
  // Dajmo vektoru koeficijenata odgovarajuću dimenziju
  coeff.resize(idxSet.size(dim));

  // Iteriramo kroz sve elemente mreže
  for (auto vertex : vertices(gridView) )
  {
         auto global = vertex.geometry().corner(0);
         // Smjesti vrijednost funkcije u centru elementa na
         // odgovarajuće mjesto u vektor.
          coeff[idxSet.index(vertex)] = exact(global);
  }
}

// Izračunaj  profilMatrice:
//  profilMatrice[i] = skup svih indeksa vrhova koji su susjedni vrhu i.
//                     Biti susjedan znači pripadati istom elementu mreže.
// GV = GridView
template<class GV>
void  izracunajProfilMatrice(GV const& gv,
                             std::vector<std::set<int> > & profilMatrice)
  {
    const int dim = GV::dimensionworld;
    // Broj vrhova mreže
    const int N = gv.size(dim);
    profilMatrice.resize(N);
    auto const & idxSet = gv.indexSet();
    // Po svim elementima mreže
    for (auto const & elem : elements(gv))
      {
        // nastavi ....        
      }
  }

// Dimenzioniraj i inicijaliziraj nulama matricu krutosti i vektor
// desne strane. Matrica krutosti se dimenzionira prema podacima u
// varijabli profilMatrice. Metoda  izracunajProfilMatrice() mora biti
// pozvana prije ove metode.
template<typename GV, typename Mat, typename Vec>
void init(GV const & gv, Mat & A, Vec & b, std::vector<std::set<int> > & profilMatrice)
{
    const int dim = GV::dimensionworld;
    // Broj vrhova mreže
    const int N = gv.size(dim);

    // dimenzije matrice A i vektora desne strane b
    A.setSize(N, N);
    A.setBuildMode(Mat::random);
    b.resize(N);

    // nastavi .....
    
    // Profil matrice je time fiksiran
    // inicijalizacija nulama
    A = 0.0;
    b = 0.0;
}

// Asembliranje matrice krutosti i vektora desne strane.
// Profil matrice mora biti definiran i matrica i vektor desne strane moraju
// biti inicijalizirani nulama.
// GV = LeafGridView
// FEM = Lokalni konačni element
// Mat = tip matrice
// Vec = tip vektora
template<typename GV, typename FEM, typename Mat, typename Vec>
void assemble(GV & gv, FEM const & fem, Mat & A, Vec & b){
    // A = matrica krutosti, b = vektor desne strane
    const int  dim = GV::dimensionworld;
    using FEMGradient = typename FEM::Traits::LocalBasisType::Traits::JacobianType;
    using FEMRange    = typename FEM::Traits::LocalBasisType::Traits::RangeType;

    auto const & idxSet = gv.indexSet();
    // Petlja po svim elementima
    for(auto const & element : elements(gv))
    {
        int basis_size = fem.localBasis().size();
        auto const & geo = element.geometry();
        // kvadraturna formula
        const auto& rule = Dune::QuadratureRules<double,dim>::rule(element.type(),2);

        for (auto const & r : rule)
        {
            auto ip = r.position(); // integration point
            // (grad G)^{-t}(xi)
            auto jacInvTra = geo.jacobianInverseTransposed(ip);
            // integracijska težina wi
            auto weight = r.weight();
            // | det (grad G)|
            auto detjac = geo.integrationElement(ip);
            // Bazne funkcije
            std::vector<FEMRange> phi(basis_size);
            fem.localBasis().evaluateFunction(ip, phi);
            // grad phi (za sve bazne funkcije phi)
            std::vector<FEMGradient> gradphi(basis_size);
            fem.localBasis().evaluateJacobian(ip, gradphi);
           
            // nastavi ....
        }
    } // Kraj petlje po elementima

    // Dirichletov rubni uvjet
    // Svaku Dirichletovu točku zamijenjujemo trivijalnim linijama.
    for(auto const & element : elements(gv))
    {
        for (auto const & is : intersections(gv,element))
        {
            // jesmo li na granici
            if ( is.boundary() )
            {
                // instanciramo referentni element
                const auto ref = Dune::referenceElement<double,dim>(element.type());

                // Obiđimo sve vrhove na stranici elementa
                auto isInd = is.indexInInside(); // indeks stranice u elementu
                for (int i=0; i < ref.size(isInd,1,dim); i++)
                {
                    // nastavi ....
                }
            }
        }
    }
}


int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 2;
  using GridType = Dune::UGGrid<dim>;
  using GridView = GridType::LeafGridView ;

  // Učitaj mrežu
  std::unique_ptr<GridType> gridptr = Dune::GmshReader<GridType>::read("src_dir/square.msh");
  int nref = 0;
  if(argc >1) nref = std::stoi(argv[1]);
  gridptr->globalRefine(nref);
  auto gv = gridptr->leafGridView();
  std::cout << " Učitavanje mreže je gotovo.\n";

  // Linearna algebra
  using Matrix = Dune::BCRSMatrix<double>;
  using Vector = Dune::BlockVector<double>;
  using FEM = Dune::LagrangeSimplexLocalFiniteElement<double, double, dim, 1>;

 
  // nastavi ...


  Vector ex; // vektor točnog rješenja
  vertexdata(gv, ex);

  Dune::VTKWriter<GridView> vtkwriter(gv);
  vtkwriter.addVertexData(x, "u");
  vtkwriter.addVertexData(ex, "exact");
  Vector diff(ex);
  diff -= x;
  vtkwriter.addVertexData(diff, "diff");
  vtkwriter.write("fem", Dune::VTK::OutputType::ascii);

  ex-=x;
  std::cout << " Greška aproksimacije (L^infty) = " <<ex.infinity_norm() << std::endl;
  return 0;
}
