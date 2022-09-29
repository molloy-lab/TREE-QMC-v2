#include "instance.hpp"

long long total_count[8];

int main(int argc, char **argv) {
    Instance instance(argc, argv);
    long long time = instance.solve();
    if (time >= 0) std::cout << time << "ms" << std::endl;
    if (instance.get_solution() != NULL) 
        instance.output_solution();
    return 0;
}
