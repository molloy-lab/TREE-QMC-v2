#include "instance.hpp"

Instance::Instance(int argc, char **argv) {
    input_file = output_file = "";
    normal = "2"; execute = "0"; taxa_mode = "0";
    refine_seed = 12345; cut_seed = 1;
    dict = NULL; output = NULL;
    if (parse(argc, argv)) {
        std::cout << help_info;
    }
    else {
        if (input_file == "") {
            std::cout << "input file not found" << std::endl;
            std::cout << "use -h for more information" << std::endl;
        }
        else {
            dict = new Dict;
            if (input_trees() == 0) {
                dict->update_singletons();
                if (verbose > "0") {
                    subproblem_csv.open(output_file + "_subproblems.csv");
                    subproblem_csv << "ID,PARENT,DEPTH,SIZE,ARTIFICIAL,SUBSET";
                    if (verbose > "1") {
                        subproblem_csv << ",ENTRY,PRUNED,ZERO";
                    }
                    subproblem_csv << std::endl;
                }
                std::cout << "Input has " << input.size() << " gene trees and " << dict->size() << " taxa.\n";
                if (execute == "0") resolve_trees();
            }
        }
    }
}

Instance::~Instance() {
    for (Tree *t : input) delete t;
    delete output;
    delete dict;
}

long long Instance::solve() {
    if (input.size() == 0) {
        return -1;
    }
    else {
        srand(cut_seed);
        std::string mode = normal + execute + taxa_mode;
        auto start = std::chrono::high_resolution_clock::now();
        output = new SpeciesTree(input, dict, mode);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        return duration.count();
    }
}

SpeciesTree *Instance::get_solution() {
    return output;
}

void Instance::output_solution() {
    if (output_file == "") {
        std::cout << output->to_string() << std::endl;
    }
    else {
        std::ofstream fout(output_file);
        fout << output->to_string() << std::endl;
    }
}

bool Instance::parse(int argc, char **argv) {
    std::cout << "TREE-QMC version 2.0.0\nCOMMAND: ";
    for (int j = 0; j < argc; j ++) 
        std::cout << argv[j] << " ";
    std::cout << std::endl;
    int i = 0;
    bool help = false;
    while (i < argc) {
        std::string opt(argv[i]);
        if (opt == "-h" || opt == "--help") help = true;
        if (opt == "-i" || opt == "--input") input_file = argv[++ i];
        if (opt == "-o" || opt == "--output") output_file = argv[++ i];
        if (opt == "-n" || opt == "--normalize") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "ERROR: invalid normalize parameter!" << std::endl;
                return true;
            }
            normal = param;
        }
        if (opt == "-x" || opt == "--execution") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1") {
                std::cout << "ERROR: invalid execution parameter!" << std::endl;
                return true;
            }
            execute = param;
        }
        if (opt == "-v" || opt == "--verbose") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "ERROR: invalid verbose parameter!" << std::endl;
                return true;
            }
            verbose = param;
        }
        if (opt == "--polyseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "-1") {
                for (int j = 0; j < param.length(); j ++) {
                    if (param[j] < '0' || param[j] > '9') {
                        std::cout << "ERROR: invalid polyseed parameter!" << std::endl;
                        return true;
                    }
                }
            }
            refine_seed = std::stoi(param);
            if (refine_seed < 0) refine_seed = time(0);
        }
        if (opt == "--cutseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "-1") {
                for (int j = 0; j < param.length(); j ++) {
                    if (param[j] < '0' || param[j] > '9') {
                        std::cout << "ERROR: invalid maxcutseed parameter!" << std::endl;
                        return true;
                    }
                }
            }
            cut_seed = std::stoi(param);
            if (cut_seed < 0) cut_seed = time(0);
        }
        if (opt == "--shared") taxa_mode = "1";
        i ++;
    }
    std::cout << "normalization scheme: n" + normal << std::endl;
    if (execute == "0") 
        std::cout << "execution mode: efficient" << std::endl;
    else 
        std::cout << "execution mode: brute force" << std::endl;
    std::cout << "random seed for refinement: " << refine_seed << std::endl;
    std::cout << "random seed for max-cut: " << cut_seed << std::endl;
    if (taxa_mode == "1") 
        std::cout << "WARNING: importance values are shared across gene trees!" << std::endl;
    return help;
}

std::size_t Instance::input_trees() {
    std::ifstream fin(input_file);
    if (! fin.is_open()) {
        std::cout << "ERROR: input file " << input_file << " does not exist!" << std::endl;
        return 1;
    }
    std::string newick;
    while (std::getline(fin, newick)) {
        if (newick.find(";") == std::string::npos) break;
        Tree *t = new Tree(newick, dict);
        input.push_back(t);
    }
    fin.close();
    return 0;
}

void Instance::resolve_trees() {
    size_t total = 0;
    for (Tree *t : input) 
        total += t->resolve();
    std::cout << "Performed " << total << " refinements." << std::endl;
    if (total != 0) {
        std::ofstream fout(input_file + ".refined");
        for (Tree *t : input) 
            fout << t->to_string() << std::endl;
        fout.close();
    }
}
