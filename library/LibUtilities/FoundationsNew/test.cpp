#include "Points.hpp"

using namespace Nektar::LibUtilities::Foundations;

int main () {
    Points<Quadrilateral, Tensor, double> P(5,5);

    // This one should cause a compile error
    //Points<Quadrilateral, Fekete, double> Pfail(5,5);

    Points<Triangle, Tensor, double> P2(5,5);
}
