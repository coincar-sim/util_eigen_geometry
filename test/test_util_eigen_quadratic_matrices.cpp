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

#include "gtest/gtest.h"

#include "util_eigen_quadratic_matrices.hpp"

using namespace util_eigen_quadratic_matrices;

using eigen_covariance_t = Matrix6d;


class UtilEigenQuadraticMatrices : public ::testing::Test {
protected:
    UtilEigenQuadraticMatrices() {
        covNonSetValuesEigen_ = eigen_covariance_t::Identity(6, 6) * -1;
        covGroundTruthValuesEigen_ = eigen_covariance_t::Zero(6, 6);
        covArbitraryDiagValuesValidEigen_ << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0;
        covArbitraryValuesValidEigen_ << 11, 12, 13, 14, 15, 16, 12, 22, 23, 24, 25, 26, 13, 23, 33, 34, 35, 36, 14, 24,
            34, 44, 45, 46, 15, 25, 35, 45, 55, 56, 16, 26, 36, 46, 56, 66;
        covInvalidUnsymmetricEigen_ << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
            23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36;
        covInvalidBelowZeroEigen_ << -11, 12, 13, 14, 15, 16, 12, 22, 23, 24, 25, 26, 13, 23, 33, 34, 35, 36, 14, 24,
            34, 44, 45, 46, 15, 25, 35, 45, 55, 56, 16, 26, 36, 46, 56, 66;
    }

    eigen_covariance_t covNonSetValuesEigen_;
    eigen_covariance_t covGroundTruthValuesEigen_;
    eigen_covariance_t covArbitraryDiagValuesValidEigen_;
    eigen_covariance_t covArbitraryValuesValidEigen_;
    eigen_covariance_t covInvalidUnsymmetricEigen_;
    eigen_covariance_t covInvalidBelowZeroEigen_;
};

TEST_F(UtilEigenQuadraticMatrices, covarianceMatrixIsSymmetric) {
    EXPECT_TRUE(isSymmetric(covNonSetValuesEigen_));
    EXPECT_TRUE(isSymmetric(covGroundTruthValuesEigen_));
    EXPECT_TRUE(isSymmetric(covArbitraryDiagValuesValidEigen_));
    EXPECT_TRUE(isSymmetric(covArbitraryValuesValidEigen_));
    EXPECT_FALSE(isSymmetric(covInvalidUnsymmetricEigen_));
    EXPECT_TRUE(isSymmetric(covInvalidBelowZeroEigen_));
}

TEST_F(UtilEigenQuadraticMatrices, covarianceMatrixIsPosDefinit) {
    EXPECT_FALSE(isPosSemiDefinit(covNonSetValuesEigen_));
    EXPECT_TRUE(isPosSemiDefinit(covGroundTruthValuesEigen_));
    EXPECT_TRUE(isPosSemiDefinit(covArbitraryDiagValuesValidEigen_));
    EXPECT_TRUE(isPosSemiDefinit(covArbitraryValuesValidEigen_));
    EXPECT_FALSE(isPosSemiDefinit(covInvalidUnsymmetricEigen_));
    EXPECT_FALSE(isPosSemiDefinit(covInvalidBelowZeroEigen_));
}
