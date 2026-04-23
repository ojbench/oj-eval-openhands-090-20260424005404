#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"

class matrix {
private:
    int m, n;
    fraction **data;

public:
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    matrix(int m_, int n_) {
        if (m_ <= 0 || n_ <= 0) {
            throw matrix_error();
        }
        m = m_;
        n = n_;
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = fraction(0);
            }
        }
    }

    matrix(const matrix &obj) {
        m = obj.m;
        n = obj.n;
        if (obj.data == nullptr) {
            data = nullptr;
        } else {
            data = new fraction*[m];
            for (int i = 0; i < m; i++) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; j++) {
                    data[i][j] = obj.data[i][j];
                }
            }
        }
    }

    matrix(matrix &&obj) noexcept {
        m = obj.m;
        n = obj.n;
        data = obj.data;
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
    }

    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    matrix &operator=(const matrix &obj) {
        if (this != &obj) {
            if (data != nullptr) {
                for (int i = 0; i < m; i++) {
                    delete[] data[i];
                }
                delete[] data;
            }
            
            m = obj.m;
            n = obj.n;
            if (obj.data == nullptr) {
                data = nullptr;
            } else {
                data = new fraction*[m];
                for (int i = 0; i < m; i++) {
                    data[i] = new fraction[n];
                    for (int j = 0; j < n; j++) {
                        data[i][j] = obj.data[i][j];
                    }
                }
            }
        }
        return *this;
    }

    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n || data == nullptr) {
            throw matrix_error();
        }
        return data[i-1][j];
    }

    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m || lhs.data == nullptr || rhs.data == nullptr) {
            throw matrix_error();
        }
        
        matrix result(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; i++) {
            for (int j = 0; j < rhs.n; j++) {
                fraction sum(0);
                for (int k = 0; k < lhs.n; k++) {
                    sum = sum + lhs.data[i][k] * rhs.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }
        return result;
    }

    matrix transposition() {
        if (data == nullptr || m <= 0 || n <= 0) {
            throw matrix_error();
        }
        
        matrix result(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    fraction determination() {
        if (data == nullptr || m <= 0 || n <= 0 || m != n) {
            throw matrix_error();
        }
        
        matrix temp(*this);
        fraction det(1);
        
        for (int i = 0; i < m; i++) {
            int pivot_row = i;
            for (int j = i + 1; j < m; j++) {
                if (!(temp.data[j][i] == fraction(0)) && 
                    ((temp.data[pivot_row][i] == fraction(0)) || 
                     !(temp.data[j][i] == temp.data[pivot_row][i]))) {
                    pivot_row = j;
                }
            }
            
            if (temp.data[pivot_row][i] == fraction(0)) {
                return fraction(0);
            }
            
            if (pivot_row != i) {
                for (int j = 0; j < n; j++) {
                    fraction temp_frac = temp.data[i][j];
                    temp.data[i][j] = temp.data[pivot_row][j];
                    temp.data[pivot_row][j] = temp_frac;
                }
                det = det * fraction(-1);
            }
            
            for (int j = i + 1; j < m; j++) {
                if (!(temp.data[j][i] == fraction(0))) {
                    fraction factor = temp.data[j][i] / temp.data[i][i];
                    for (int k = i; k < n; k++) {
                        temp.data[j][k] = temp.data[j][k] - factor * temp.data[i][k];
                    }
                }
            }
        }
        
        for (int i = 0; i < m; i++) {
            det = det * temp.data[i][i];
        }
        
        return det;
    }
};

class resistive_network {
private:
    int interface_size, connection_size;
    matrix adjacency, conduction;
    int *from_nodes, *to_nodes;
    fraction *resistances;

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        interface_size = interface_size_;
        connection_size = connection_size_;
        
        from_nodes = new int[connection_size];
        to_nodes = new int[connection_size];
        resistances = new fraction[connection_size];
        
        for (int i = 0; i < connection_size; i++) {
            from_nodes[i] = from[i];
            to_nodes[i] = to[i];
            resistances[i] = resistance[i];
        }
        
        adjacency = matrix(interface_size, interface_size);
        conduction = matrix(connection_size, connection_size);
        
        for (int i = 0; i < connection_size; i++) {
            int u = from[i] - 1;
            int v = to[i] - 1;
            
            fraction conductance = fraction(1) / resistance[i];
            adjacency(u + 1, v) = conductance;
            adjacency(v + 1, u) = conductance;
        }
        
        for (int i = 0; i < connection_size; i++) {
            conduction(i + 1, i) = fraction(1) / resistance[i];
        }
    }

    ~resistive_network() {
        delete[] from_nodes;
        delete[] to_nodes;
        delete[] resistances;
    }

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        if (interface_id1 == interface_id2) {
            return fraction(0);
        }
        
        matrix laplacian(interface_size, interface_size);
        
        for (int i = 0; i < interface_size; i++) {
            fraction degree(0);
            for (int j = 0; j < interface_size; j++) {
                if (i != j) {
                    degree = degree + adjacency(i + 1, j);
                    laplacian(i + 1, j) = fraction(0) - adjacency(i + 1, j);
                }
            }
            laplacian(i + 1, i) = degree;
        }
        
        matrix reduced_laplacian(interface_size - 1, interface_size - 1);
        for (int i = 0; i < interface_size - 1; i++) {
            for (int j = 0; j < interface_size - 1; j++) {
                reduced_laplacian(i + 1, j) = laplacian(i + 1, j);
            }
        }
        
        fraction det = reduced_laplacian.determination();
        
        int i_idx = interface_id1 - 1;
        int j_idx = interface_id2 - 1;
        
        matrix minor_ii(interface_size - 2, interface_size - 2);
        int row = 0, col = 0;
        for (int r = 0; r < interface_size - 1; r++) {
            if (r == i_idx) continue;
            col = 0;
            for (int c = 0; c < interface_size - 1; c++) {
                if (c == i_idx) continue;
                minor_ii(row + 1, col) = reduced_laplacian(r + 1, c);
                col++;
            }
            row++;
        }
        
        matrix minor_jj(interface_size - 2, interface_size - 2);
        row = 0; col = 0;
        for (int r = 0; r < interface_size - 1; r++) {
            if (r == j_idx) continue;
            col = 0;
            for (int c = 0; c < interface_size - 1; c++) {
                if (c == j_idx) continue;
                minor_jj(row + 1, col) = reduced_laplacian(r + 1, c);
                col++;
            }
            row++;
        }
        
        matrix minor_ij(interface_size - 2, interface_size - 2);
        row = 0; col = 0;
        for (int r = 0; r < interface_size - 1; r++) {
            if (r == i_idx) continue;
            col = 0;
            for (int c = 0; c < interface_size - 1; c++) {
                if (c == j_idx) continue;
                minor_ij(row + 1, col) = reduced_laplacian(r + 1, c);
                col++;
            }
            row++;
        }
        
        fraction det_ii, det_jj, det_ij;
        if (interface_size > 2) {
            det_ii = minor_ii.determination();
            det_jj = minor_jj.determination();
            det_ij = minor_ij.determination();
        } else {
            det_ii = fraction(1);
            det_jj = fraction(1);
            det_ij = fraction(0);
        }
        
        int sign_ii = ((i_idx + i_idx) % 2 == 0) ? 1 : -1;
        int sign_jj = ((j_idx + j_idx) % 2 == 0) ? 1 : -1;
        int sign_ij = ((i_idx + j_idx) % 2 == 0) ? 1 : -1;
        
        fraction numerator = fraction(sign_ii) * det_ii + fraction(sign_jj) * det_jj - fraction(2 * sign_ij) * det_ij;
        
        return numerator / det;
    }

    fraction get_voltage(int id, fraction current[]) {
        matrix laplacian(interface_size, interface_size);
        
        for (int i = 0; i < interface_size; i++) {
            fraction degree(0);
            for (int j = 0; j < interface_size; j++) {
                if (i != j) {
                    degree = degree + adjacency(i + 1, j);
                    laplacian(i + 1, j) = fraction(0) - adjacency(i + 1, j);
                }
            }
            laplacian(i + 1, i) = degree;
        }
        
        matrix reduced_laplacian(interface_size - 1, interface_size - 1);
        for (int i = 0; i < interface_size - 1; i++) {
            for (int j = 0; j < interface_size - 1; j++) {
                reduced_laplacian(i + 1, j) = laplacian(i + 1, j);
            }
        }
        
        matrix current_vector(interface_size - 1, 1);
        for (int i = 0; i < interface_size - 1; i++) {
            current_vector(i + 1, 0) = current[i];
        }
        
        fraction det = reduced_laplacian.determination();
        
        if (id == interface_size) {
            return fraction(0);
        }
        
        matrix temp_matrix = reduced_laplacian;
        for (int i = 0; i < interface_size - 1; i++) {
            temp_matrix(i + 1, id - 1) = current_vector(i + 1, 0);
        }
        
        fraction numerator_det = temp_matrix.determination();
        
        return numerator_det / det;
    }

    fraction get_power(fraction voltage[]) {
        fraction power(0);
        
        for (int i = 0; i < connection_size; i++) {
            int u = from_nodes[i] - 1;
            int v = to_nodes[i] - 1;
            
            fraction voltage_diff = voltage[u] - voltage[v];
            fraction connection_power = (voltage_diff * voltage_diff) / resistances[i];
            power = power + connection_power;
        }
        
        return power;
    }
};

#endif //SRC_HPP
