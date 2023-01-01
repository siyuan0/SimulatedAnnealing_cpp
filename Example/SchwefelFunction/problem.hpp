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

#define DIMENSION 2 // to change this define for number of dimensions considered

namespace Schwefel
{
    static int num_of_evaluations = 0; // track the total number of evaluations of Schwefel's function
    static std::mt19937 randomGen;

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

void setRandomGen(std::mt19937& gen)
{
    randomGen = gen;
}

float l2(soln& s1, soln& s2)
{  // get the l2 norm of s1-s2
    float sum = 0; 
    for(int i=0; i<DIMENSION; i++) sum += std::pow(s1.getX(i) - s2.getX(i), 2);
    return std::pow(sum, 0.5);
}

SA_policy<soln> initialiseRuntimeInfo(std::unordered_map<std::string, float>& parameters)
{
    soln initialMaxChange{};
    for(int i=0; i<DIMENSION; i++) initialMaxChange.setX(i, parameters["initial max change"]);
    return {
        .temperature = parameters["intial temperature"],
        .maxChange = initialMaxChange
    };
}

soln getRandomSolution(std::unordered_map<std::string, float>& parameters)
{   // return a random solution within problem constraints
    soln s{parameters["min xi"], parameters["max xi"]};
    return s;
}

soln getNewSolution(std::unordered_map<std::string, float>& parameters,
                    SA_policy<soln>& runtimeInfo, soln& currSoln)
{

}


soln getBestSoln(std::vector<soln>& population)
{   // return the best soln in the population
    int idx_best = 0;
    for(int i=0; i<population.size(); i++)
    {
        if(population[i].getEval() < population[idx_best].getEval()) idx_best = i;
    }
    return population[idx_best];
}

bool endSearch(std::unordered_map<std::string, float>& parameters)
{   // end when computational budget is exceeded
    return Schwefel::num_of_evaluations > parameters["max_eval"];
}

// store the problem specific methods for the GA core to run on
static ProblemCtx<soln> problemCtx = {
    .setRandomGenerator = &setRandomGen,
    .initRuntimeInfo = &initialiseRuntimeInfo,
    .getRandomSolution = nullptr,
    .getNewSolution = nullptr,
    .acceptProbability = nullptr,
    .updateRuntimeInfo = nullptr,
    .compareSoln = nullptr,
    .endSearch = nullptr
};

} // namespace Schwefel

#endif // INCLUDE_SCHWEFEL