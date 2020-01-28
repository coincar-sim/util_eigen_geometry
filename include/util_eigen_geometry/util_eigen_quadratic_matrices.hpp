/*
 * Copyright (c) 2017
 * FZI Forschungszentrum Informatik, Karlsruhe, Germany (www.fzi.de)
 * KIT, Institute of Measurement and Control, Karlsruhe, Germany (www.mrt.kit.edu)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace util_eigen_quadratic_matrices {
/* Extend convenience typedefs of predefined Eigen::Matrix.
 */
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Vector6d = Eigen::Matrix<double, 6, 1>;

/** @brief Removes the i-th row and the i-th column from a symmetric matrix. Specialization for dynamic-size matrices.
 * @param[in] MxM quadratic input Matrix
 * @param[in] Index of row and column that will be removed
 * @return (M-1) x (M-1) quadratic output Matrix*/
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> removeRowCol(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const int n) {
    const int size = static_cast<int>(input.rows());
    if (input.rows() != input.cols()) {
        throw std::invalid_argument("Method is only defined for quadratic matrices.");
    }
    if (size == 0) {
        throw std::invalid_argument(
            "Method is not applicable to matrices of size 0x0. Can not remove rows and columns if there are none.");
    }
    if (size == 1) {
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(0, 0);
    }
    // Shortcut if n = 0 or n = (size -1) to avoid failing Assertion in older Eigen Versions.
    // In older Eigen Versions, the comma-initializer can not handle zero-size vectors or matrices.
    // This is fixed in newer Eigen Versions (>3.2.10 ?).
    if (n == 0) {
        return input.bottomRightCorner(size - 1, size - 1);
    }
    if (n == (size - 1)) {
        return input.topLeftCorner(size - 1, size - 1);
    }
    return (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(size - 1, size - 1) << input.block(0, 0, n, n),
            input.block(0, n + 1, n, size - n - 1),
            input.block(n + 1, 0, size - n - 1, n),
            input.block(n + 1, n + 1, size - n - 1, size - n - 1))
        .finished();
}

/** @brief Checks if a quadratic matrix is symmetric or not.
 * @param[in] MxM quadratic input Matrix
 * @return true if the matrix is symmetric or false otherwise.*/
template <typename T, int Size>
inline bool isSymmetric(const Eigen::Matrix<T, Size, Size>& input) {
    if (input.rows() != input.cols()) {
        throw std::invalid_argument("Method is only defined for quadratic matrices.");
    }
    for (int i = 1; i < input.rows(); i++) {
        if (!(input.diagonal(i).isApprox(input.diagonal(-i)))) {
            return false;
        }
    }
    return true;
}

/** @brief Checks if a quadratic matrix is positive (semi) definite or not.
 * @param[in] MxM quadratic input Matrix
 * @return true if the matrix is positive (semi) definite or false otherwise.*/
template <typename T, int Size>
inline bool isPosSemiDefinit(const Eigen::Matrix<T, Size, Size>& input) {
    if (input.rows() != input.cols()) {
        throw std::invalid_argument("Method is only defined for quadratic matrices.");
    }
    if (input.rows() == 0) {
        return false; // special case for 0D-zero-elements Matrix. Prevents Seg-Fault.
    }
    Eigen::LDLT<Eigen::Matrix<T, Size, Size>> chol(input);
    return (chol.isPositive());
}
} // namespace util_eigen_quadratic_matrices
