#include "instance.hpp"

Instance::Instance(int argc, char **argv) {
    input_file = output_file = "";
    normal = "2"; execute = "0"; verbose = "0";
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
            input_trees();
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
        std::string mode = normal + execute + verbose;
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
        i ++;
    }
    return help;
}

void Instance::input_trees() {
    srand(refine_seed);
    std::ifstream fin(input_file);
    std::string newick;
    while (std::getline(fin, newick)) {
        if (newick.find(";") == std::string::npos) break;
        Tree *t = new Tree(newick, dict);
        t->resolve();
        input.push_back(t);
    }
    fin.close();
}
