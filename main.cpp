#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <chrono>
#include "lib/core.hpp"
#include "third_party/nlohmann/json.hpp"
#include "Example/SchwefelFunction/problem.hpp"

int main(int argc, 
         char *argv[]) {
    if(argc<=1)
    {
        std::cout << "missing paramters.json file\n";
    }else if(argc==2)
    {
        std::ifstream f(argv[1]);
        nlohmann::json data = nlohmann::json::parse(f);
        auto jmap = data.get<std::unordered_map<std::string, float>>();

        SA<Schwefel::soln> SAinst(Schwefel::problemCtx, jmap);
        auto start = std::chrono::high_resolution_clock::now();
        SAinst.optimise();
        auto finish = std::chrono::high_resolution_clock::now();
        std::cout << "Optimisation took " << 
                std::chrono::duration_cast<std::chrono::milliseconds>(finish-start).count() << "ms\n";
        if(jmap["print results"])
        {
            std::cout << "results saved to allSolutions.txt and acceptedSolutions.txt\n";
            SAinst.printAllToFile("allSolutions.txt");
            SAinst.printAcceptedToFile("acceptedSolutions.txt");
        }
        std::cout << "number of function evaluations: " << Schwefel::num_of_evaluations << '\n';
        std::cout << "final temperature: " << SAinst.getRuntimeInfo().temperature << '\n';
        std::cout << "current solution: " << SAinst.getOptimisationResult().first.print() << '\n';
        std::cout << "best solution: " << SAinst.getOptimisationResult().second.print() << '\n';
    }else
    {
        std::cout << "too many arguments\n";
    }

    return 0;
}