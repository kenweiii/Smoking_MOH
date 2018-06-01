#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix simforward(NumericMatrix tm, IntegerVector population, int years)
{
    //one column represents the population for a given year
    IntegerMatrix finalpop(population.length(),years);
    double thisprob;
    for(int year = 0; year < years;year++)
    {
        //iterate through each member of the population and check their new state after a year
        for(int i = 0;i < population.length();i++)
        {
            thisprob = runif(1)[0];
            switch(population[i])
            {
                case 1:
                    if(thisprob < (tm(0,3) + tm(0,1) + tm(0,0))) population[i] = 4;
                    if(thisprob < (tm(0,1) + tm(0,0))) population[i] = 2;
                    if(thisprob < tm(0,0)) population[i] = 1;
                    break;
                case 2:
                    if(thisprob < (tm(1,4) + tm(1,3) + tm(1,2) + tm(1,1))) population[i] = 5;
                    if(thisprob < (tm(1,3) + tm(1,2) + tm(1,1))) population[i] = 4;
                    if(thisprob < (tm(1,2) + tm(1,1))) population[i] = 3;
                    if(thisprob < tm(1,1)) population[i] = 2;
                    break;
                case 3:
                    if(thisprob < (tm(2,3) + tm(2,2) + tm(2,1))) population[i] = 4;
                    if(thisprob < (tm(2,2) + tm(2,1))) population[i] = 3;
                    if(thisprob < tm(2,1)) population[i] = 2;
                    break;
                case 4:
                    if(thisprob < (tm(3,4) + tm(3,3) + tm(3,2) + tm(3,1))) population[i] = 5;
                    if(thisprob < (tm(3,3) + tm(3,2) + tm(3,1))) population[i] = 4;
                    if(thisprob < (tm(3,2) + tm(3,1))) population[i] = 3;
                    if(thisprob < tm(3,1)) population[i] = 2;
                    break;
                case 5:
                    if(thisprob < (tm(4,4) + tm(4,3) + tm(4,2) + tm(4,1))) population[i] = 5;
                    if(thisprob < (tm(4,3) + tm(4,2) + tm(4,1))) population[i] = 4;
                    if(thisprob < (tm(4,2) + tm(4,1))) population[i] = 3;
                    if(thisprob < tm(4,1)) population[i] = 2;
                    break;
                default:
                    break;
            }
        }
        finalpop.column(year) = population;
    }
    return finalpop;
}