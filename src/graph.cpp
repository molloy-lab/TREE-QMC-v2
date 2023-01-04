#include "graph.hpp"

extern unsigned long long count[8];
extern std::ofstream count_csv;

Graph::Graph(std::vector<Tree *> trees, Taxa &subset) {
    size = subset.size();
    for (index_t i = 0; i < size; i ++) {
        index2index[subset.root_at(i)] = i;
        indices.push_back(subset.root_at(i));
    }
    graph = new weight_t**[2];
    graph[0] = Matrix::new_mat(size);
    graph[1] = Matrix::new_mat(size);
    if (! PRUNING) {
        count[0] = count[1] = count[2] = 0;
    }
    for (Tree *tree : trees) {
        std::unordered_map<index_t, index_t> valid = tree->get_indices();
        subset.weight_update(valid);
        weight_t ***subgraph = tree->build_graph(subset);
        for (index_t i = 0; i < size; i ++) {
            for (index_t j = 0; j < size; j ++) {
                if (subgraph[0][i][j] > 0 || subgraph[1][i][j] > 0) {
                    index_t i_ = index2index[subset.root_at(i)];
                    index_t j_ = index2index[subset.root_at(j)];
                    graph[0][i_][j_] += subgraph[0][i][j];
                    graph[1][i_][j_] += subgraph[1][i][j];
                }
            }
        }
        Matrix::delete_mat(subgraph[0], size);
        Matrix::delete_mat(subgraph[1], size);
        delete [] subgraph;
    }
    if (! PRUNING) {
        count_csv << count[0] << ',' << count[1] << ',' << count[2] << ',';
        count_csv << "\"" + subset.to_list() + "\"" << std::endl;
    }
}

Graph::Graph(std::unordered_map<quartet_t, weight_t> &quartets, Taxa &subset) {
    size = subset.size();
    for (index_t i = 0; i < size; i ++) {
        index2index[subset.root_at(i)] = i;
        indices.push_back(subset.root_at(i));
    }
    graph = new weight_t**[2];
    graph[0] = Matrix::new_mat(size);
    graph[1] = Matrix::new_mat(size);
    for (auto elem : quartets) {
        index_t *indices = split(elem.first);
        index_t a = index2index[indices[0]], b = index2index[indices[1]], c = index2index[indices[2]], d = index2index[indices[3]];
        weight_t w = elem.second;
        graph[1][a][b] += w; graph[1][c][d] += w; graph[1][b][a] += w; graph[1][d][c] += w;
        graph[0][a][c] += w; graph[0][a][d] += w; graph[0][b][c] += w; graph[0][b][d] += w;
        graph[0][c][a] += w; graph[0][d][a] += w; graph[0][c][b] += w; graph[0][d][b] += w;
        delete [] indices;
    }
    std::cout << subset.to_string() << std::endl;
}

Graph::~Graph() {
    Matrix::delete_mat(graph[0], size);
    Matrix::delete_mat(graph[1], size);
    delete[] graph;
}

std::string Graph::to_string() {
    return Matrix::display_mat(graph[0], size) + "\n" + Matrix::display_mat(graph[1], size);
}

weight_t Graph::get_cut(std::vector<index_t> *A, std::vector<index_t> *B) {
    weight_t positive_weight = -1.0;
    std::vector<index_t> a, b;
    weight_t lower = 0.0, upper = 6.0;
    /*
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            if (graph[1][i][j] == 0) continue;
            weight_t ratio = graph[0][i][j] / graph[1][i][j];
            if (ratio > upper) upper = ratio;
        }
    }
    */
    // std::cout << upper << std::endl;
    while (lower + 0.1 < upper) {
        weight_t alpha = (lower + upper) / 2.0;
        a.clear(); b.clear();
        weight_t weight = sdp_cut(alpha, &a, &b);
        if (weight < 0.001 || a.size() <= 1 || b.size() <= 1) {
            upper = alpha;
        }
        else {
            lower = alpha;
            positive_weight = alpha;
            *A = a;
            *B = b;
        }
    }
    /*
    if (A->size() <= 1 || B->size() <= 1) {
        std::cout << Matrix::display_mat(graph[0], size) << std::endl;
        std::cout << Matrix::display_mat(graph[1], size) << std::endl;
    }
    assert(A->size() > 1 && B->size() > 1);
    */
    return positive_weight;
}

weight_t Graph::sdp_cut(weight_t alpha, std::vector<index_t> *A, std::vector<index_t> *B) {
    std::vector<Instance::InstanceTuple> input;
    weight_t avg = 0, num = size * (size - 1) / 2;
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            weight_t temp = (graph[0][i][j] - alpha * graph[1][i][j]) / num;
            if (temp < 0) temp = - temp;
            avg += temp;
        }
    }
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            weight_t weight = (graph[0][i][j] - alpha * graph[1][i][j]) / avg;
            input.push_back(Instance::InstanceTuple(std::make_pair(i + 1, j + 1), weight));
        }
    }
    MaxCutInstance instance(input, size);
    Burer2002 heuristic(instance, -1, false, NULL);
    MaxCutSimpleSolution solution = heuristic.get_best_solution();
    std::vector<int> cut = solution.get_assignments();
    for (index_t i = 0; i < cut.size(); i ++) {
        if (cut[i] < 0) 
            A->push_back(indices[i]);
        else 
            B->push_back(indices[i]);
    }
    return solution.get_weight();
}
