#include "instance.hpp"

std::ofstream count_csv("zero-count.csv", std::ofstream::out);
unsigned long long count[8];

int main(int argc, char **argv) {
    if (! PRUNING) count_csv << "ALL,PRUNED,ZEROES,SUBPROBLEM\n";
    Instance instance(argc, argv);
    long long time = instance.solve();
    if (! PRUNING) std::cout << count[3] << ',' << count[4] << ',' << count[5] << std::endl;
    if (time >= 0) std::cout << time << "ms" << std::endl;
    if (instance.get_solution() != NULL) 
        instance.output_solution();
    return 0;
}
