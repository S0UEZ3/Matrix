#pragma once
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

template <class T>
class Matrix {
public:
    Matrix() = delete;
    Matrix(const std::size_t n, const std::size_t m, const T default_value);
    Matrix(const std::size_t n, const std::size_t m);
    // Constructor for square matrix
    Matrix(const std::size_t n);
    Matrix(const Matrix&) = default;

    T& operator()(const std::size_t x, const std::size_t y);
    const T& operator()(const std::size_t x, const std::size_t y) const; // Added const overload
    Matrix operator+ (const Matrix& other) const;
    Matrix& operator++ ();
    Matrix operator++ (int);
    bool operator== (const Matrix& other) const;
    bool operator!= (const Matrix& other) const;
    Matrix& operator+= (const Matrix& other);
    Matrix operator- (const Matrix& other) const;
    Matrix& operator-= (const Matrix& other);
    Matrix operator* (const Matrix& other) const;
    Matrix& operator*= (const Matrix& other);

    std::pair<std::size_t, std::size_t> size() const;
    Matrix pow(int exp) const;
    T det() const;
    Matrix inverse() const;
    Matrix t() const;

private:
    std::size_t n, m;
    std::vector<std::vector<T>> mtx;

    friend std::ostream& operator << (std::ostream&, const Matrix&);
};

template<class T>
std::ostream& operator << (std::ostream& os, const Matrix<T>& matrix);


// -----------------------------------------------------------------//

template<class T>
Matrix<T>::Matrix(const std::size_t n, const std::size_t m, const T default_value) :
    n(n), m(m), mtx(n, std::vector<T>(m, default_value)) {}

template<class T>
Matrix<T>::Matrix(const std::size_t n, const std::size_t m) :
    n(n), m(m), mtx(n, std::vector<T>(m)) {}

template<class T>
Matrix<T>::Matrix(const std::size_t n) :
    n(n), m(n), mtx(n, std::vector<T>(n, 0)) {
    for (std::size_t i = 0; i < n; ++i) {
        mtx[i][i] = 1;
    }
}

template<class T>
T& Matrix<T>::operator()(const std::size_t x, const std::size_t y) {
    if (x >= n || y >= m || x < 0 || y < 0) {
        throw std::out_of_range("Out of bounds");
    }
    return mtx[x][y];
}

template<class T>
const T& Matrix<T>::operator()(const std::size_t x, const std::size_t y) const { // Const overload
    if (x >= n || y >= m || x < 0 || y < 0) {
        throw std::out_of_range("Out of bounds");
    }
    return mtx[x][y];
}

template<class T>
Matrix<T>& Matrix<T>::operator+= (const Matrix<T>& other) {
    if (n != other.n || m != other.m) {
        throw std::logic_error("Shapes not match");
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            mtx[i][j] += other.mtx[i][j];
        }
    }
    return *this;
};

template<class T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& other) const {
    Matrix<T> summ = *this;
    summ += other;
    return summ;
};

template<class T>
Matrix<T>& Matrix<T>::operator++ () {
    if (n != m) {
        throw std::logic_error("Not square matrix");
    }
    for (std::size_t i = 0; i < n; ++i) {
        mtx[i][i]++;
    }
    return *this;
};

template<class T>
Matrix<T> Matrix<T>::operator++ (int) {
    Matrix<T> copy = *this;
    ++(*this);
    return copy;
};

template<class T>
bool Matrix<T>::operator== (const Matrix<T>& other) const {
    if (n != other.n || m != other.m) {
        return false;
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (mtx[i][j] != other.mtx[i][j]) {
                return false;
            }
        }
    }
    return true;
};

template<class T>
bool Matrix<T>::operator!= (const Matrix<T>& other) const {
    return !((*this) == other);
};

template<class T>
Matrix<T> Matrix<T>::operator- (const Matrix<T>& other) const {
    Matrix<T> diff = *this;
    diff -= other;
    return diff;
}

template<class T>
Matrix<T>& Matrix<T>::operator-= (const Matrix<T>& other) {
    if (n != other.n || m != other.m) {
        throw std::logic_error("Shapes not match");
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            mtx[i][j] -= other.mtx[i][j];
        }
    }
    return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator* (const Matrix<T>& other) const {
    Matrix<T> product = *this;
    product *= other;
    return product;
}

template<class T>
Matrix<T>& Matrix<T>::operator*= (const Matrix<T>& other) {
    if (m != other.n) {
        throw std::logic_error("Shapes not match");
    }
    Matrix<T> result(n, other.m, 0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < other.m; ++j) {
            for (std::size_t k = 0; k < m; ++k) {
                result.mtx[i][j] += mtx[i][k] * other.mtx[k][j];
            }
        }
    }
    *this = result;
    return *this;
}

template<class T>
std::pair<std::size_t, std::size_t> Matrix<T>::size() const {
    return std::make_pair(n, m);
};

template<class T>
Matrix<T> Matrix<T>::pow(int s) const {
    if (n != m) {
        throw std::logic_error("Not square matrix");
    }
    Matrix<T> result(n);
    Matrix<T> base = *this;
    while (s > 0) {
        if (s % 2 == 1) {
            result *= base;
        }
        base *= base;
        s /= 2;
    }
    return result;
}

template<class T>
T Matrix<T>::det() const {
    if (n != m) {
        throw std::logic_error("Not square matrix");
    }
    Matrix<T> temp(*this);
    T determinant = 1.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (temp(i, i) == 0) {
            for (std::size_t j = i + 1; j < n; ++j) {
                if (temp(j, i) != 0) {
                    for (std::size_t k = 0; k < n; ++k) {
                        std::swap(temp(i, k), temp(j, k));
                    }
                    determinant = -determinant;
                    break;
                }
            }
            if (temp(i, i) == 0) {
                return 0;
            }
        }
        determinant *= temp(i, i);
        for (std::size_t j = i + 1; j < n; ++j) {
            T factor = temp(j, i) / temp(i, i);
            for (std::size_t k = i; k < n; ++k) {
                temp(j, k) -= factor * temp(i, k);
            }
        }
    }
    return determinant;
}

template<class T>
Matrix<T> Matrix<T>::inverse() const {
    if (n != m) {
        throw std::logic_error("Not a square matrix");
    }

    Matrix<T> inv(n, n);
    Matrix<T> gs(*this);

    // Initialize inv as the identity matrix
    for (std::size_t i = 0; i < n; ++i) {
        inv(i, i) = 1.0;
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (gs(i, i) == 0) {
            bool swapped = false;
            for (std::size_t j = i + 1; j < n; ++j) {
                if (gs(j, i) != 0) {
                    for (std::size_t k = 0; k < n; ++k) {
                        std::swap(gs(i, k), gs(j, k));
                        std::swap(inv(i, k), inv(j, k));
                    }
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                throw std::logic_error("Singular matrix");
            }
        }

        T diag_val = gs(i, i);
        for (std::size_t j = 0; j < n; ++j) {
            gs(i, j) /= diag_val;
            inv(i, j) /= diag_val;
        }

        for (std::size_t j = 0; j < n; ++j) {
            if (i != j) {
                T f = gs(j, i);
                for (std::size_t k = 0; k < n; ++k) {
                    gs(j, k) -= f * gs(i, k);
                    inv(j, k) -= f * inv(i, k);
                }
            }
        }
    }
    return inv;
}



template<class T>
Matrix<T> Matrix<T>::t() const {
    Matrix<T> tr(m, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            tr.mtx[j][i] = mtx[i][j];
        }
    }
    return tr;
};

template<class T>
std::ostream& operator << (std::ostream& os, const Matrix<T>& matrix) {
    os << matrix.size().first << " " << matrix.size().second << '\n';
    for (std::size_t i = 0; i < matrix.size().first; ++i) {
        for (std::size_t j = 0; j < matrix.size().second; ++j) {
            os << matrix.mtx[i][j] << " ";
        }
        os << '\n';
    }
    return os;
}