#include <cstdint>
#include <cmath>
#include <iostream>
#include <gsl/>

double V_LJ(double r)
{
    return (4.0 * (std::pow(r,-12.0)-std::pow(r,-6.0)));
}

int main()
{
    std::cout << "Ciao ese1!" << std::endl;
    
    constexpr double eps = 5.9e-3; // eV
    constexpr double sigma = 3.18; // Ams


    
    return 0;
}