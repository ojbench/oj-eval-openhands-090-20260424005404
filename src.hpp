#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"

class matrix {
private:

    // m行n列的矩阵，用动态二维数组存储，每个元素是分数类实例
    int m, n;
    fraction **data;

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************

public:

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************

    // 默认构造函数
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    // 构造函数，构建 m_*n_ 的矩阵，矩阵元素设为0。
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

    // 拷贝构造函数，构建与 obj 完全相同的矩阵。
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

    // 移动拷贝构造函数。
    matrix(matrix &&obj) noexcept {
        m = obj.m;
        n = obj.n;
        data = obj.data;
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
    }

    // 析构函数。
    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    // 重载赋值号。
    matrix &operator=(const matrix &obj) {
        if (this != &obj) {
            // 释放原有内存
            if (data != nullptr) {
                for (int i = 0; i < m; i++) {
                    delete[] data[i];
                }
                delete[] data;
            }
            
            // 复制新数据
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

    // 重载括号，返回矩阵的第i行(1-based)、第j列(0-based)的元素的引用。如果 i、j 不合法，抛出 matrix_error 错误。
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n || data == nullptr) {
            throw matrix_error();
        }
        return data[i-1][j];
    }

    // 重载乘号，返回矩阵乘法 lhs * rhs 的结果。如果 lhs 的列数与 rhs 的行数不相等，抛出 matrix_error 错误。
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

    // 返回矩阵的转置。若矩阵为空，抛出 matrix_error 错误。
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

    // 返回矩阵的行列式。建议用高斯消元实现。若矩阵不是方阵或为空，抛出 matrix_error 错误。
    fraction determination() {
        if (data == nullptr || m <= 0 || n <= 0 || m != n) {
            throw matrix_error();
        }
        
        // 创建矩阵的副本用于高斯消元
        matrix temp(*this);
        fraction det(1);
        
        for (int i = 0; i < m; i++) {
            // 找到主元
            int pivot_row = i;
            for (int j = i + 1; j < m; j++) {
                if (!(temp.data[j][i] == fraction(0)) && 
                    ((temp.data[pivot_row][i] == fraction(0)) || 
                     !(temp.data[j][i] == temp.data[pivot_row][i]))) {
                    pivot_row = j;
                }
            }
            
            // 如果主元为0，行列式为0
            if (temp.data[pivot_row][i] == fraction(0)) {
                return fraction(0);
            }
            
            // 交换行（如果需要）
            if (pivot_row != i) {
                for (int j = 0; j < n; j++) {
                    fraction temp_frac = temp.data[i][j];
                    temp.data[i][j] = temp.data[pivot_row][j];
                    temp.data[pivot_row][j] = temp_frac;
                }
                det = det * fraction(-1); // 交换行改变行列式符号
            }
            
            // 消元
            for (int j = i + 1; j < m; j++) {
                if (!(temp.data[j][i] == fraction(0))) {
                    fraction factor = temp.data[j][i] / temp.data[i][i];
                    for (int k = i; k < n; k++) {
                        temp.data[j][k] = temp.data[j][k] - factor * temp.data[i][k];
                    }
                }
            }
        }
        
        // 计算对角线元素的乘积
        for (int i = 0; i < m; i++) {
            det = det * temp.data[i][i];
        }
        
        return det;
    }
};

class resistive_network {
private:

    // 节点数量 和 接线数量
    int interface_size, connection_size;

    // 矩阵A 和 矩阵C
    matrix adjacency, conduction;

    // 存储连接信息用于功率计算
    int *from_nodes, *to_nodes;
    fraction *resistances;

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************

public:

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************

    // 设置电阻网络，构建矩阵A和C。节点数量为interface_size_，接线数量为connection_size_。
    //       对于 1<=i<=connection_size_，从节点from[i-1]到节点to[i-1]有接线，对应电阻为resistance[i-1]。
    //       保证接线使得电阻网络联通，from[i-1] < to[i-1]，resitance[i-1] > 0，均合法。
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        interface_size = interface_size_;
        connection_size = connection_size_;
        
        // 存储连接信息
        from_nodes = new int[connection_size];
        to_nodes = new int[connection_size];
        resistances = new fraction[connection_size];
        
        for (int i = 0; i < connection_size; i++) {
            from_nodes[i] = from[i];
            to_nodes[i] = to[i];
            resistances[i] = resistance[i];
        }
        
        // 初始化邻接矩阵A (n x n) 和传导矩阵C (m x m)
        adjacency = matrix(interface_size, interface_size);
        conduction = matrix(connection_size, connection_size);
        
        // 构建邻接矩阵A
        for (int i = 0; i < connection_size; i++) {
            int u = from[i] - 1; // 转换为0-based
            int v = to[i] - 1;   // 转换为0-based
            
            // 邻接矩阵：A[u][v] = A[v][u] = 1/r
            fraction conductance = fraction(1) / resistance[i];
            adjacency(u + 1, v) = conductance;
            adjacency(v + 1, u) = conductance;
        }
        
        // 构建传导矩阵C (对角矩阵，C[i][i] = 1/r_i)
        for (int i = 0; i < connection_size; i++) {
            conduction(i + 1, i) = fraction(1) / resistance[i];
        }
    }

    ~resistive_network() {
        delete[] from_nodes;
        delete[] to_nodes;
        delete[] resistances;
    }

    // 返回节点 interface_id1 和 interface_id2 (1-based)之间的等效电阻。
    //       保证 interface_id1 <= interface_id2 均合法。
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        if (interface_id1 == interface_id2) {
            return fraction(0); // 同一节点间电阻为0
        }
        
        // 构建拉普拉斯矩阵 L = D - A，其中D是度数矩阵
        matrix laplacian(interface_size, interface_size);
        
        // 计算度数矩阵并构建拉普拉斯矩阵
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
        
        // 移除最后一行和最后一列（参考节点）
        matrix reduced_laplacian(interface_size - 1, interface_size - 1);
        for (int i = 0; i < interface_size - 1; i++) {
            for (int j = 0; j < interface_size - 1; j++) {
                reduced_laplacian(i + 1, j) = laplacian(i + 1, j);
            }
        }
        
        // 使用正确的等效电阻公式
        // R_eq = L^+_ii + L^+_jj - 2*L^+_ij
        // 其中 L^+ 是拉普拉斯矩阵的伪逆
        
        // 对于小规模网络，我们使用余子式方法
        // 等效电阻 = (余子式(i,i) + 余子式(j,j) - 2*余子式(i,j)) / det(L_reduced)
        
        fraction det = reduced_laplacian.determination();
        
        // 计算所需的余子式
        int i_idx = interface_id1 - 1;
        int j_idx = interface_id2 - 1;
        
        // 余子式(i,i)
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
        
        // 余子式(j,j)
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
        
        // 余子式(i,j)
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
        
        // 考虑符号 (-1)^(i+j)
        int sign_ii = ((i_idx + i_idx) % 2 == 0) ? 1 : -1;
        int sign_jj = ((j_idx + j_idx) % 2 == 0) ? 1 : -1;
        int sign_ij = ((i_idx + j_idx) % 2 == 0) ? 1 : -1;
        
        fraction numerator = fraction(sign_ii) * det_ii + fraction(sign_jj) * det_jj - fraction(2 * sign_ij) * det_ij;
        
        return numerator / det;
    }

    // 在给定节点电流I的前提下，返回节点id(1-based)的电压。认为节点interface_size(1-based)的电压为0。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电流为 current[i-1]。
    //       保证 current 使得电阻网络有解，id < interface_size 合法。
    fraction get_voltage(int id, fraction current[]) {
        // 构建拉普拉斯矩阵 L = D - A
        matrix laplacian(interface_size, interface_size);
        
        // 计算度数矩阵并构建拉普拉斯矩阵
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
        
        // 移除最后一行和最后一列（参考节点，电压为0）
        matrix reduced_laplacian(interface_size - 1, interface_size - 1);
        for (int i = 0; i < interface_size - 1; i++) {
            for (int j = 0; j < interface_size - 1; j++) {
                reduced_laplacian(i + 1, j) = laplacian(i + 1, j);
            }
        }
        
        // 构建电流向量（移除最后一个元素）
        matrix current_vector(interface_size - 1, 1);
        for (int i = 0; i < interface_size - 1; i++) {
            current_vector(i + 1, 0) = current[i];
        }
        
        // 使用克莱姆法则求解线性系统 L_reduced * V = I
        fraction det = reduced_laplacian.determination();
        
        if (id == interface_size) {
            return fraction(0); // 参考节点电压为0
        }
        
        // 创建矩阵，将第id列替换为电流向量
        matrix temp_matrix = reduced_laplacian;
        for (int i = 0; i < interface_size - 1; i++) {
            temp_matrix(i + 1, id - 1) = current_vector(i + 1, 0);
        }
        
        fraction numerator_det = temp_matrix.determination();
        
        return numerator_det / det;
    }


    // 在给定节点电压U的前提下，返回电阻网络的功率。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电压为 voltage[i-1]。
    //       保证 voltage 合法。
    fraction get_power(fraction voltage[]) {
        // 功率 P = Σ (V_i - V_j)^2 / R_ij，对所有边(i,j)求和
        
        fraction power(0);
        
        // 遍历所有连接
        for (int i = 0; i < connection_size; i++) {
            int u = from_nodes[i] - 1; // 转换为0-based
            int v = to_nodes[i] - 1;   // 转换为0-based
            
            // 计算电压差
            fraction voltage_diff = voltage[u] - voltage[v];
            // 功率 = (电压差)^2 / 电阻
            fraction connection_power = (voltage_diff * voltage_diff) / resistances[i];
            power = power + connection_power;
        }
        
        return power;
    }
};


#endif //SRC_HPP