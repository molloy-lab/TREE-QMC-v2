#include "instance.hpp"

long long total_count[8];

int main(int argc, char **argv) {
    Instance instance(argc, argv);
    long long time = instance.solve();
    if (time >= 0) std::cout << time << "ms" << std::endl;
    if (instance.get_solution() != NULL) 
        instance.output_solution();
    // std::cout << total_count[0] << ' ' << total_count[1] << std::endl;
    // std::cout << (double) total_count[1] / (total_count[1] + total_count[0]) << std::endl;
    // std::cout << total_count[3] << std::endl;
    return 0;
}
