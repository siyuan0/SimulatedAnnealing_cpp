#ifndef INCLUDE_SCHWEFEL
#define INCLUDE_SCHWEFEL

#include "../../lib/core.hpp"
#include <cstdlib>
#include <ostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <utility>
#include <limits>

#define DIMENSION 6 // to change this define for number of dimensions considered

namespace Schwefel
{
    static int num_of_evaluations = 0; // track the total number of evaluations of Schwefel's function
    static thread_local std::mt19937 randomGen;

class soln
{
private:
    float x[DIMENSION]; // using array instead of vector to have it stored on stack for faster creation/access/deletion
    float f;
    float _lbound;
    float _ubound;

    float evaluateObjective()
    {   // evaluate Schwefel's function on this solution
        num_of_evaluations += 1;
        float tmp = 0;
        for(int i=0; i<DIMENSION; i++) 
        {
            if(x[i] < _lbound | x[i] > _ubound) return std::numeric_limits<float>::max(); // solution is outside constraints
            tmp -= x[i] * std::sin(std::sqrt(std::fabs(x[i])));
        }
        return tmp;
    }

public:
    soln(float lowerbound, float upperbound)
    {  // randomly generate a soln within the provided constraints
        _lbound = lowerbound;
        _ubound = upperbound;
        for(int i=0; i<DIMENSION; i++) x[i] = lowerbound +  float(std::rand()) / RAND_MAX * (upperbound-lowerbound);
        f = evaluateObjective();
    }

    soln()
    {   // default constructor
        _lbound = 0;
        _ubound = 0;
        for(int i=0; i<DIMENSION; i++) x[i] = 0;
        f = 0;
    }

    void doEval(){ f = evaluateObjective(); }

    float getEval(){ return f; }

    float getX(int i){ return x[i]; }

    void setX(int i, float val){ x[i]=val; }

    friend std::ostream& operator<< (std::ostream& stream, const soln& s)
    {   // for printing out the contents of a solution
        for(int i=0; i<DIMENSION; i++) stream << s.x[i] << ", ";
        stream << s.f;
        return stream;
    }

    std::string print()
    {   // same goal as operator<<, but slightly more formatted
        std::stringstream ss;
        ss << "x: [";
        for(int i=0; i<DIMENSION-1; i++) ss << x[i] << ", ";
        ss << x[DIMENSION-1] << "] f: " << f;
        return ss.str();
    }
};

float l2(soln& s1, soln& s2)
{  // get the l2 norm of s1-s2
    float sum = 0; 
    for(int i=0; i<DIMENSION; i++) sum += std::pow(s1.getX(i) - s2.getX(i), 2);
    return std::pow(sum, 0.5);
}

void setThreadRandomGenerator(std::mt19937 gen){randomGen = gen;}

soln getBestSoln(std::vector<soln>& population)
{   // return the best soln in the population
    int idx_best = 0;
    for(int i=0; i<population.size(); i++)
    {
        if(population[i].getEval() < population[idx_best].getEval()) idx_best = i;
    }
    return population[idx_best];
}

std::vector<soln> getInitialPopulation(int size, std::unordered_map<std::string, float>& parameters)
{   // randomly initialise initial population
    std::vector<soln> v;
    for(int i=0; i<size; i++) v.push_back(soln(
        parameters["min xi"],
        parameters["max xi"]
    ));
    return v;
}

bool endSearch(std::unordered_map<std::string, float>& parameters)
{   // end when computational budget is exceeded
    return Schwefel::num_of_evaluations > parameters["max_eval"];
}

// store the problem specific methods for the GA core to run on
static ProblemCtx<soln> problemCtx = {
    .setRandomGenerator = &setThreadRandomGenerator,
    .getRandomSolutions = &getInitialPopulation,
    .getParentIdx = &getParentIdx,
    .getChildren = &getChildren,
    .updatePopulation = &updatePopulation,
    .endSearch = &endSearch
};

} // namespace Schwefel

#endif // INCLUDE_SCHWEFEL