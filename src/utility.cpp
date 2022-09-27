#include "utility.hpp"

weight_t **Matrix::new_mat(index_t size) {
    weight_t **m = new weight_t*[size];
    for (index_t i = 0; i < size; i ++) {
        m[i] = new weight_t[size];
        for (index_t j = 0; j < size; j ++) {
            m[i][j] = 0;
        }
    }
    return m;
}

void Matrix::delete_mat(weight_t **m, index_t size) {
    for (index_t i = 0; i < size; i ++) {
        delete[] m[i];
    }
    delete[] m;
}

std::string Matrix::display_mat(weight_t **m, index_t size) {
    std::stringstream ss;
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = 0; j < size; j ++) {
            ss << std::setw(12) << std::setprecision(6) << m[i][j];
        }
        ss << std::endl;
    }
    return ss.str();
}

weight_t Matrix::diff_mat(weight_t **m1, weight_t **m2, index_t size) {
    weight_t sum = 0;
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = 0; j < size; j ++) {
            weight_t delta = m1[i][j] - m2[i][j];
            if (delta < 0) delta = - delta;
            sum += delta;
        }
    }
    return sum;
}
