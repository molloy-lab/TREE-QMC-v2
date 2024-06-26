#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include "utility.hpp"
#include "tree.hpp"

class Instance {
    public:
        Instance(int argc, char **argv);
        ~Instance();
        long long solve();
        SpeciesTree *get_solution();
        void output_solution();
    private:
        std::string input_file, output_file, normal, execute, taxa_mode;
        unsigned int refine_seed, cut_seed;
        std::vector<Tree *> input;
        Dict *dict;
        SpeciesTree *output;
        bool parse(int argc, char **argv);
        std::size_t input_trees();
        void resolve_trees();
};

extern std::ofstream subproblem_csv;
extern std::string verbose;
extern unsigned long long count[8];

#endif
