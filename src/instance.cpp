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
            if (verbose > "0") {
                subproblem_csv.open(output_file + "_subproblems.csv");
                subproblem_csv << "ID,PARENT,DEPTH,SIZE,ARTIFICIAL,SUBSET";
                if (verbose > "1") {
                    subproblem_csv << ",ENTRY,PRUNED,ZERO";
                }
                subproblem_csv << std::endl;
            }
            dict = new Dict;
            input_trees();
            dict->update_singletons();
            if (execute == "0") resolve_trees();
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
        if (opt == "-n" || opt == "--normalize") normal = argv[++ i];
        if (opt == "-x" || opt == "--execution") execute = argv[++ i];
        if (opt == "-v" || opt == "--verbose") verbose = argv[++ i];
        if (opt == "--polyseed") refine_seed = std::stoi(argv[++ i]);
        if (opt == "--cutseed") cut_seed = std::stoi(argv[++ i]);
        if (opt == "--shared") taxa_mode = "1";
        i ++;
    }
    return help;
}

void Instance::input_trees() {
    std::ifstream fin(input_file);
    std::string newick;
    while (std::getline(fin, newick)) {
        if (newick.find(";") == std::string::npos) break;
        Tree *t = new Tree(newick, dict);
        input.push_back(t);
    }
    fin.close();
}

void Instance::resolve_trees() {
    size_t total = 0;
    for (Tree *t : input) 
        total += t->resolve();
    if (total != 0) {
        std::cout << total << " fake nodes" << std::endl;
        std::ofstream fout(input_file + ".refined");
        for (Tree *t : input) 
            fout << t->to_string() << std::endl;
        fout.close();
    }
}
