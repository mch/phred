#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

int main () {
  
  boost::numeric::ublas::matrix<double> m (3, 3);
  boost::numeric::ublas::vector<double> v (3);

  for (unsigned i = 0; i < std::min (m.size1 (), v.size ()); ++ i) {
    for (unsigned j = 0; j <= i; ++ j)
      m (i, j) = 3 * i + j + 1;
    v (i) = i;
  }

  std::cout << solve (m, v, boost::numeric::ublas::lower_tag ()) 
            << std::endl;
  std::cout << solve (v, m, boost::numeric::ublas::lower_tag ()) 
            << std::endl;
}

