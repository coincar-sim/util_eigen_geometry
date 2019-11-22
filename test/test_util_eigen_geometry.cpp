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

#include "util_eigen_geometry.hpp"

using namespace util_eigen_geometry;


TEST(UtilEigenGeometry, positiveFloatModulo) {
    EXPECT_DOUBLE_EQ(0., positiveFloatModulo(0., 0.)) << "float modulo of (0.,0.) is defined to be zero";

    EXPECT_DOUBLE_EQ(0., positiveFloatModulo(10., 1.));
    EXPECT_DOUBLE_EQ(10., positiveFloatModulo(10., 0.));
    EXPECT_DOUBLE_EQ(1., positiveFloatModulo(10., 4.5));

    EXPECT_DOUBLE_EQ(5., positiveFloatModulo(35., 10.));
}

TEST(UtilEigenGeometry, angleDifferenceDouble) {
    EXPECT_NEAR(M_PI / 4., angleDifference(123. + M_PI / 4., 123.), 10.e-9);
    EXPECT_NEAR(-M_PI / 4., angleDifference(123. - M_PI / 4., 123.), 10.e-9);
}


TEST(UtilEigenGeometry, angleDifferenceDegrees) {
    EXPECT_DOUBLE_EQ(40., angleDifferenceDegrees(400., 360.));
    EXPECT_DOUBLE_EQ(-40., angleDifferenceDegrees(400., 440.));
}

TEST(UtilEigenGeometry, cosineSimilarityDouble) {
    EXPECT_DOUBLE_EQ(1., cosineSimilarity(360., 360.));
    EXPECT_DOUBLE_EQ(cos(40.), cosineSimilarity(400., 360.));
    EXPECT_DOUBLE_EQ(cos(-40.), cosineSimilarity(400., 440.));
}

TEST(UtilEigenGeometry, yawFromAffine3d) {
    Eigen::Isometry3d transformation{Eigen::Isometry3d::Identity()};
    transformation.rotate(Eigen::AngleAxisd(M_PI / 4., Eigen::Vector3d::UnitZ()));
    double yaw = yawFromAffine3d(transformation);
    EXPECT_DOUBLE_EQ(M_PI / 4., yaw);

    transformation.rotate(Eigen::AngleAxisd(-M_PI / 2., Eigen::Vector3d::UnitZ()));
    yaw = yawFromAffine3d(transformation);
    EXPECT_DOUBLE_EQ(-M_PI / 4., yaw);
}

TEST(UtilEigenGeometry, affine3dXYFromAffine2d) {
    Eigen::Isometry2d tf(Eigen::Isometry2d::Identity());
    tf.translate(Eigen::Vector2d(1, 1));
    tf.rotate(Eigen::Rotation2Dd(M_PI_2));
    Eigen::Affine3d tf3d = affine3dXYFromAffine2d(tf);
    EXPECT_TRUE(tf3d.translation().isApprox(Eigen::Vector3d(1, 1, 0)));
    Eigen::Affine3d pi2z(Eigen::Affine3d::Identity());
    pi2z.rotate(Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitZ()));
    EXPECT_TRUE(tf3d.rotation().isApprox(pi2z.rotation()));
}

class UtilEigenGeometryPolygons : public ::testing::Test {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
    UtilEigenGeometryPolygons() {
    }

    Eigen::Vector2d p0{0., 0.};
    Eigen::Vector2d p1{1., 1.};
    Eigen::Vector2d p2{2., 2.};

    polygon_t poly01{p0, p1};
    polygon_t poly12{p1, p2};
    polygon_t poly02{p0, p1, p2};
};

TEST_F(UtilEigenGeometryPolygons, angleDifferencePoseToPolygon) {

    Eigen::Affine2d pose{Eigen::Affine2d::Identity()};
    pose.translation().setOnes();

    EXPECT_DOUBLE_EQ(M_PI / 4., angleDifference(pose, poly01));
    EXPECT_DOUBLE_EQ(M_PI / 4., angleDifference(pose, poly12));
    EXPECT_DOUBLE_EQ(M_PI / 4., angleDifference(pose, poly02));
}

TEST_F(UtilEigenGeometryPolygons, cosineSimilarityPoseToPolygon) {

    Eigen::Affine2d pose{Eigen::Affine2d::Identity()};
    pose.translation().setOnes();

    EXPECT_DOUBLE_EQ(cos(M_PI / 4.), cosineSimilarity(pose, poly01));
    EXPECT_DOUBLE_EQ(cos(M_PI / 4.), cosineSimilarity(pose, poly12));
    EXPECT_DOUBLE_EQ(cos(M_PI / 4.), cosineSimilarity(pose, poly02));
}


TEST_F(UtilEigenGeometryPolygons, getClosestId) {
    EXPECT_EQ(0ul, getClosestId(p0, poly02));
    EXPECT_EQ(2ul, getClosestId(Eigen::Vector2d{100., 100.}, poly02));
}

TEST_F(UtilEigenGeometryPolygons, lineStripOrientation) {

    Eigen::Affine2d pose{Eigen::Affine2d::Identity()};
    pose.translation().setOnes();

    EXPECT_DOUBLE_EQ(M_PI / 4., lineStripOrientation(p0, p1));
    EXPECT_DOUBLE_EQ(0., lineStripOrientation(p0, p0));
}

TEST_F(UtilEigenGeometryPolygons, splitPolygonRight) {
    polygon_t poly12new;
    splitPolygonRight(poly02, 0, poly12new);
    EXPECT_DOUBLE_EQ(poly12[0].x(), poly12new[0].x());
    EXPECT_DOUBLE_EQ(poly12[0].y(), poly12new[0].y());
}

TEST_F(UtilEigenGeometryPolygons, canSplitPolygonRight) {
    EXPECT_TRUE(canSplitPolygonRight(poly02, 1));
    EXPECT_FALSE(canSplitPolygonRight(poly02, 2));
}

TEST_F(UtilEigenGeometryPolygons, SamplePolygon) {
    polygon_t sampled = addIntermediateSamplesToPolygon(poly02, 0.1);
    ASSERT_GT(sampled.size(), 1ul);
    for (size_t i = 1; i < sampled.size(); ++i) {
        const auto dist = (sampled.at(i) - sampled.at(i - 1)).norm();
        EXPECT_LE(dist, 0.1);
        EXPECT_GE(dist, 0.05);
    }
}
