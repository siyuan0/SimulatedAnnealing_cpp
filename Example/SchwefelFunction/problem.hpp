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
        f = 0;
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

SA_policy<soln> initialiseRuntimeInfo(std::unordered_map<std::string, float>& parameters)
{
    soln initialMaxChange{};
    for(int i=0; i<DIMENSION; i++) initialMaxChange.setX(i, parameters["initial max change"]);
    return {
        .temperature = parameters["initial temperature"],
        .maxChange = initialMaxChange,
        .numAcceptedCurrTemp = 0,
        .numCurrTemp = 0,
        .numTempSteps = 1
    };
}

soln getRandomSolution(std::unordered_map<std::string, float>& parameters)
{   // return a random solution within problem constraints
    soln s{parameters["min xi"], parameters["max xi"]};
    s.doEval();
    return s;
}

soln getNewSolution(std::unordered_map<std::string, float>& parameters,
                    SA_policy<soln>& runtimeInfo, soln& currSoln)
{   // generate a new solution from the current solution using:
    // x_new = x_curr + D * u
    // where D is a diagonal matrix of max change in each dimension
    // and u is a vector of random values in [-1, 1]
    soln s{parameters["min xi"], parameters["max xi"]};
    std::uniform_real_distribution<float> urand{-1.0, 1.0};
    for(int i=0; i<DIMENSION; i++)
    {
        float newxi = parameters["max xi"] + 1;
        while((newxi > parameters["max xi"]) | (newxi < parameters["min xi"]))
            newxi = currSoln.getX(i) + urand(randomGen) * runtimeInfo.maxChange.getX(i);
        s.setX(i, newxi);
    }
    s.doEval();
    return s;
}

float acceptProbability(std::unordered_map<std::string, float>& parameters, 
                        SA_policy<soln>& runtimeInfo, soln& newSoln, soln& currSoln)
{   // get the acceptance probability of newsoln given curr soln
    // better solutions are always accepted
    return std::exp(-(newSoln.getEval() - currSoln.getEval())/runtimeInfo.temperature);
}

void updateRuntimeInfo(std::unordered_map<std::string, float>& parameters, 
                       SA_policy<soln>& runtimeInfo, soln& newSoln, soln& currSoln, bool accepted)
{
    if(accepted)
    {   // new solution is accepted, so we update the max change values
        for(int i=0; i<DIMENSION; i++) runtimeInfo.maxChange.setX(i, 
            runtimeInfo.maxChange.getX(i) * (1-parameters["alpha"]) +
            parameters["alpha"] * parameters["w"] * std::abs(newSoln.getX(i) - currSoln.getX(i))
        );
        runtimeInfo.numAcceptedCurrTemp += 1;
        runtimeInfo.numCurrTemp += 1;
    }else
    {
        runtimeInfo.numCurrTemp += 1;
    }
    if((runtimeInfo.numAcceptedCurrTemp > parameters["min accepted at each temperature"]) |
       (runtimeInfo.numCurrTemp > parameters["max same temperature chain"]))
    {   // desired length of markov chain at current temperature is reached, so we move forward
        // the annealing schedule and update the temperature
        runtimeInfo.temperature *= parameters["temperature scaling"];
        runtimeInfo.numTempSteps += 1;
        runtimeInfo.numAcceptedCurrTemp = 0;
        runtimeInfo.numCurrTemp = 0;
    }
}

bool compareSoln(soln& betterSoln, soln& worseSoln)
{   // compare if the betterSoln is really more optimal than the worseSoln
    return betterSoln.getEval() < worseSoln.getEval();
}

bool endSearch(std::unordered_map<std::string, float>& parameters, SA_policy<soln>& runtimeInfo)
{   // end the algorithm if any conditions are met
    if((Schwefel::num_of_evaluations > parameters["max eval"]) |
       (runtimeInfo.numTempSteps > parameters["max temperature steps"]))
    {
        return true;
    }else
    {
        return false;
    }
}

// store the problem specific methods for the SA core to run on
static ProblemCtx<soln> problemCtx = {
    .setRandomGenerator = &setRandomGen,
    .initRuntimeInfo = &initialiseRuntimeInfo,
    .getRandomSolution = &getRandomSolution,
    .getNewSolution = &getNewSolution,
    .acceptProbability = &acceptProbability,
    .updateRuntimeInfo = &updateRuntimeInfo,
    .compareSoln = &compareSoln,
    .endSearch = &endSearch
};

} // namespace Schwefel

#endif // INCLUDE_SCHWEFEL