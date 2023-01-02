#ifndef INCLUDE_SA_CORE
#define INCLUDE_SA_CORE

#include <random>
#include <unordered_map>
#include <string>
#include <ctime>
#include <iostream>
#include <ostream>

template <typename T>
struct SA_policy
{
    float temperature;
    T maxChange; // for creating new solutions
    int numAcceptedCurrTemp; // number of solutions accepted at current temperature
    int numCurrTemp; // number of trials at current temperature
    int numTempSteps; // number of temperature changes
};

// problem specific methods (to be defined for each optimisation problem)
template <typename T>
struct ProblemCtx
{   
    void (*setRandomGenerator)(std::mt19937&) = nullptr;
    SA_policy<T> (*initRuntimeInfo)(std::unordered_map<std::string, float>&) = nullptr;
    T (*getRandomSolution)(std::unordered_map<std::string, float>&) = nullptr;
    T (*getNewSolution)(std::unordered_map<std::string, float>&, SA_policy<T>&, T&) = nullptr;
    float (*acceptProbability)(std::unordered_map<std::string, float>&, SA_policy<T>&, T&, T&) = nullptr;
    void (*updateRuntimeInfo)(std::unordered_map<std::string, float>&, SA_policy<T>&, T&, T&, bool)=nullptr;
    bool (*compareSoln)(T&, T&) = nullptr;
    bool (*endSearch)(std::unordered_map<std::string, float>&, SA_policy<T>&) = nullptr;
};

template <typename T>
class SA
{
private:
    T _currSoln;
    T _bestSoln;
    std::vector<T> _allAcceptedSolns;
    std::vector<T> _allSolns;
    std::vector<float> _annealingSchedule;
    std::vector<float> _acceptProbs;
    
    std::unordered_map<std::string, float> _parameters;
    ProblemCtx<T> _problemCtx;
    SA_policy<T> _runtimeInfo;
    std::mt19937 _randGen;
public:
    SA(ProblemCtx<T>& problemCtx, std::unordered_map<std::string, float>& parameters)
    {
        _parameters = parameters;
        _problemCtx = problemCtx;
        _runtimeInfo = {};
       _randGen.seed(std::time(NULL));
       _problemCtx.setRandomGenerator(_randGen);
    }

    void printAllToFile(const std::string fileName)
    {   // prints entire optimisation journey to file 
        // each line is: [solution], temperature, acceptProb
        std::ofstream outfile;
        outfile.open(fileName, std::ios::out|std::ios::trunc);
        for(int i=0; i<_allSolns.size(); i++)
        {
            outfile << _allSolns[i] << ", " 
                    << _annealingSchedule[i] << ", "
                    << _acceptProbs[i] << '\n';
        }
        outfile.close();
    }

    void printAcceptedToFile(const std::string fileName)
    {   // prints only accepted solutions to file 
        // each line is: [solution]
        std::ofstream outfile;
        outfile.open(fileName, std::ios::out|std::ios::trunc);
        for(T s : _allAcceptedSolns)
        {
            outfile << s << '\n';
        }
        outfile.close();
    }

    std::pair<T,T> getOptimisationResult(){ return {_currSoln, _bestSoln}; }

    SA_policy<T> getRuntimeInfo(){ return _runtimeInfo; }

    void optimise()
    {
        // prepare for optimisation
        _currSoln = _problemCtx.getRandomSolution(_parameters);
        _runtimeInfo = _problemCtx.initRuntimeInfo(_parameters);
        _allSolns.clear();
        _allAcceptedSolns.clear();
        _annealingSchedule.clear();
        _acceptProbs.clear();
        int progressCounter = 0;
        std::uniform_real_distribution<float> uniformDist{0, 1.0};
        std::cout << "intial: " << _runtimeInfo.temperature << '\n';

        // do the optimisation
        while((progressCounter < _parameters["max iterations"]) && (!_problemCtx.endSearch(_parameters, _runtimeInfo)))
        {
            T newSoln = _problemCtx.getNewSolution(_parameters, _runtimeInfo, _currSoln);
            _allSolns.push_back(_currSoln);
            _annealingSchedule.push_back(_runtimeInfo.temperature);
            _acceptProbs.push_back(_problemCtx.acceptProbability(_parameters, _runtimeInfo, newSoln, _currSoln));
            float u = uniformDist(_randGen);
            if(u < _problemCtx.acceptProbability(_parameters, _runtimeInfo, newSoln, _currSoln))
            {
                // update runtimeinfo knowing that new solution is accepted
                _problemCtx.updateRuntimeInfo(_parameters, _runtimeInfo, newSoln, _currSoln, true);

                // update archive
                _currSoln = newSoln;
                _allAcceptedSolns.push_back(newSoln);
                if(_problemCtx.compareSoln(_currSoln, _bestSoln)) _bestSoln = _currSoln;
            }else
            {
                // update runtimeinfo knowing that new solution is rejected
                _problemCtx.updateRuntimeInfo(_parameters, _runtimeInfo, newSoln, _currSoln, false);
            }
        }
    }
};

#endif // INCLUDE_SA_CORE