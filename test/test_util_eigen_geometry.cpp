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


TEST(UtilEigenGeometry, positiveFloatModulo) { // NOLINT
    EXPECT_DOUBLE_EQ(0., positiveFloatModulo(0., 0.)) << "float modulo of (0.,0.) is defined to be zero";

    EXPECT_DOUBLE_EQ(0., positiveFloatModulo(10., 1.));
    EXPECT_DOUBLE_EQ(10., positiveFloatModulo(10., 0.));
    EXPECT_DOUBLE_EQ(1., positiveFloatModulo(10., 4.5));

    EXPECT_DOUBLE_EQ(5., positiveFloatModulo(35., 10.));
}

TEST(UtilEigenGeometry, angleDifferenceDouble) { // NOLINT
    EXPECT_NEAR(M_PI / 4., angleDifference(123. + M_PI / 4., 123.), 10.e-9);
    EXPECT_NEAR(-M_PI / 4., angleDifference(123. - M_PI / 4., 123.), 10.e-9);
}


TEST(UtilEigenGeometry, angleDifferenceDegrees) { // NOLINT
    EXPECT_DOUBLE_EQ(40., angleDifferenceDegrees(400., 360.));
    EXPECT_DOUBLE_EQ(-40., angleDifferenceDegrees(400., 440.));
}

TEST(UtilEigenGeometry, cosineSimilarityDouble) { // NOLINT
    EXPECT_DOUBLE_EQ(1., cosineSimilarity(360., 360.));
    EXPECT_DOUBLE_EQ(cos(40.), cosineSimilarity(400., 360.));
    EXPECT_DOUBLE_EQ(cos(-40.), cosineSimilarity(400., 440.));
}

TEST(UtilEigenGeometry, yawFromAffine3d) { // NOLINT
    Eigen::Isometry3d transformation{Eigen::Isometry3d::Identity()};
    transformation.rotate(Eigen::AngleAxisd(M_PI / 4., Eigen::Vector3d::UnitZ()));
    double yaw = yawFromAffine3d(transformation);
    EXPECT_DOUBLE_EQ(M_PI / 4., yaw);

    transformation.rotate(Eigen::AngleAxisd(-M_PI / 2., Eigen::Vector3d::UnitZ()));
    yaw = yawFromAffine3d(transformation);
    EXPECT_DOUBLE_EQ(-M_PI / 4., yaw);
}

TEST(UtilEigenGeometry, affine3dXYFromAffine2d) { // NOLINT
    Eigen::Isometry2d tf(Eigen::Isometry2d::Identity());
    tf.translate(Eigen::Vector2d(1, 1));
    tf.rotate(Eigen::Rotation2Dd(M_PI_2));
    Eigen::Affine3d tf3d = affine3dXYFromAffine2d(tf);
    EXPECT_TRUE(tf3d.translation().isApprox(Eigen::Vector3d(1, 1, 0)));
    Eigen::Affine3d pi2z(Eigen::Affine3d::Identity());
    pi2z.rotate(Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitZ()));
    EXPECT_TRUE(tf3d.rotation().isApprox(pi2z.rotation()));
}

TEST(UtilEigenGeometry, isIsometry) { // NOLINT
    // Check 2d.
    {
        Eigen::Affine2d t = Eigen::Affine2d::Identity();
        t.translation().setRandom();

        t.linear() = Eigen::Rotation2Dd(0.34).toRotationMatrix();
        EXPECT_TRUE(isIsometry(t));

        t(0, 0) = std::abs(t(0, 0)) + 5.0;
        EXPECT_FALSE(isIsometry(t));
    }

    // Check 3d.
    {
        Eigen::Affine3d t = Eigen::Affine3d::Identity();
        t.translation().setRandom();

        t.linear() = Eigen::AngleAxisd(0.3, Eigen::Vector3d{1.0, -0.2, 1.0}.normalized()).toRotationMatrix();
        EXPECT_TRUE(isIsometry(t));

        t(0, 0) = std::abs(t(0, 0)) + 5.0;
        EXPECT_FALSE(isIsometry(t));
    }
}

class UtilEigenGeometryOrientationVector : public ::testing::Test {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
    UtilEigenGeometryOrientationVector() {
        isometry_180 = isometry2dFromVector(vector_180);
        isometry_135 = isometry2dFromVector(vector_135);
        isometry_90 = isometry2dFromVector(vector_90);
        isometry_45 = isometry2dFromVector(vector_45);
        isometry0 = isometry2dFromVector(vector0);
        isometry45 = isometry2dFromVector(vector45);
        isometry90 = isometry2dFromVector(vector90);
        isometry135 = isometry2dFromVector(vector135);
        isometry180 = isometry2dFromVector(vector180);
    };

    Eigen::Vector2d vector_180{-1., 0.};  // NOLINT
    Eigen::Vector2d vector_135{-1., -1.}; // NOLINT
    Eigen::Vector2d vector_90{0., -1.};   // NOLINT
    Eigen::Vector2d vector_45{1., -1.};   // NOLINT
    Eigen::Vector2d vector0{1., 0};       // NOLINT
    Eigen::Vector2d vector45{1., 1.};     // NOLINT
    Eigen::Vector2d vector90{0., 1.};     // NOLINT
    Eigen::Vector2d vector135{-1., 1.};   // NOLINT
    Eigen::Vector2d vector180{-1., 0.};   // NOLINT

    Eigen::Isometry2d isometry_180; // NOLINT
    Eigen::Isometry2d isometry_135; // NOLINT
    Eigen::Isometry2d isometry_90;  // NOLINT
    Eigen::Isometry2d isometry_45;  // NOLINT
    Eigen::Isometry2d isometry0;    // NOLINT
    Eigen::Isometry2d isometry45;   // NOLINT
    Eigen::Isometry2d isometry90;   // NOLINT
    Eigen::Isometry2d isometry135;  // NOLINT
    Eigen::Isometry2d isometry180;  // NOLINT
};

TEST_F(UtilEigenGeometryOrientationVector, isometry2dToAndFromVector) { // NOLINT
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, yawFromIsometry2d(isometry_180));
    EXPECT_DOUBLE_EQ(-3. / 4. * M_PI, yawFromIsometry2d(isometry_135));
    EXPECT_DOUBLE_EQ(-2. / 4. * M_PI, yawFromIsometry2d(isometry_90));
    EXPECT_DOUBLE_EQ(-1. / 4. * M_PI, yawFromIsometry2d(isometry_45));
    EXPECT_DOUBLE_EQ(0., yawFromIsometry2d(isometry0));
    EXPECT_DOUBLE_EQ(1. / 4. * M_PI, yawFromIsometry2d(isometry45));
    EXPECT_DOUBLE_EQ(2. / 4. * M_PI, yawFromIsometry2d(isometry90));
    EXPECT_DOUBLE_EQ(3. / 4. * M_PI, yawFromIsometry2d(isometry135));
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, yawFromIsometry2d(isometry180));

    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, yawFromVector(vectorFromIsometry2d(isometry_180)));
    EXPECT_DOUBLE_EQ(-3. / 4. * M_PI, yawFromVector(vectorFromIsometry2d(isometry_135)));
    EXPECT_DOUBLE_EQ(-2. / 4. * M_PI, yawFromVector(vectorFromIsometry2d(isometry_90)));
    EXPECT_DOUBLE_EQ(-1. / 4. * M_PI, yawFromVector(vectorFromIsometry2d(isometry_45)));
    EXPECT_DOUBLE_EQ(0., yawFromVector(vectorFromIsometry2d(isometry0)));
    EXPECT_DOUBLE_EQ(1. / 4. * M_PI, yawFromVector(vectorFromIsometry2d(isometry45)));
    EXPECT_DOUBLE_EQ(2. / 4. * M_PI, yawFromVector(vectorFromIsometry2d(isometry90)));
    EXPECT_DOUBLE_EQ(3. / 4. * M_PI, yawFromVector(vectorFromIsometry2d(isometry135)));
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, yawFromVector(vectorFromIsometry2d(isometry180)));
}

TEST_F(UtilEigenGeometryOrientationVector, yawFromVector) { // NOLINT
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, yawFromVector(vector_180));
    EXPECT_DOUBLE_EQ(-3. / 4. * M_PI, yawFromVector(vector_135));
    EXPECT_DOUBLE_EQ(-2. / 4. * M_PI, yawFromVector(vector_90));
    EXPECT_DOUBLE_EQ(-1. / 4. * M_PI, yawFromVector(vector_45));
    EXPECT_DOUBLE_EQ(0., yawFromVector(vector0));
    EXPECT_DOUBLE_EQ(1. / 4. * M_PI, yawFromVector(vector45));
    EXPECT_DOUBLE_EQ(2. / 4. * M_PI, yawFromVector(vector90));
    EXPECT_DOUBLE_EQ(3. / 4. * M_PI, yawFromVector(vector135));
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, yawFromVector(vector180));
}

TEST_F(UtilEigenGeometryOrientationVector, angleDifferenceIsometryToVector) { // NOLINT
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(isometry_180, vector0));
    EXPECT_DOUBLE_EQ(-3. / 4. * M_PI, angleDifference(isometry_135, vector0));
    EXPECT_DOUBLE_EQ(-2. / 4. * M_PI, angleDifference(isometry_90, vector0));
    EXPECT_DOUBLE_EQ(-1. / 4. * M_PI, angleDifference(isometry_45, vector0));
    EXPECT_DOUBLE_EQ(0., angleDifference(isometry0, vector0));
    EXPECT_DOUBLE_EQ(1. / 4. * M_PI, angleDifference(isometry45, vector0));
    EXPECT_DOUBLE_EQ(2. / 4. * M_PI, angleDifference(isometry90, vector0));
    EXPECT_DOUBLE_EQ(3. / 4. * M_PI, angleDifference(isometry135, vector0));
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(isometry180, vector0));

    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(isometry_135, vector45));
    EXPECT_DOUBLE_EQ(1. / 4. * M_PI, angleDifference(isometry90, vector45));

    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(isometry0, vector_180));
    EXPECT_DOUBLE_EQ(2. / 4. * M_PI, angleDifference(isometry0, vector_90));
    EXPECT_DOUBLE_EQ(-2. / 4. * M_PI, angleDifference(isometry0, vector90));
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(isometry0, vector180));
}

TEST_F(UtilEigenGeometryOrientationVector, angleDifferenceVectorToVector) { // NOLINT
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(vector_180, vector0));
    EXPECT_DOUBLE_EQ(-3. / 4. * M_PI, angleDifference(vector_135, vector0));
    EXPECT_DOUBLE_EQ(-2. / 4. * M_PI, angleDifference(vector_90, vector0));
    EXPECT_DOUBLE_EQ(-1. / 4. * M_PI, angleDifference(vector_45, vector0));
    EXPECT_DOUBLE_EQ(0., angleDifference(vector0, vector0));
    EXPECT_DOUBLE_EQ(1. / 4. * M_PI, angleDifference(vector45, vector0));
    EXPECT_DOUBLE_EQ(2. / 4. * M_PI, angleDifference(vector90, vector0));
    EXPECT_DOUBLE_EQ(3. / 4. * M_PI, angleDifference(vector135, vector0));
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(vector180, vector0));

    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(vector_135, vector45));
    EXPECT_DOUBLE_EQ(1. / 4. * M_PI, angleDifference(vector90, vector45));

    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(vector0, vector_180));
    EXPECT_DOUBLE_EQ(2. / 4. * M_PI, angleDifference(vector0, vector_90));
    EXPECT_DOUBLE_EQ(-2. / 4. * M_PI, angleDifference(vector0, vector90));
    EXPECT_DOUBLE_EQ(-4. / 4. * M_PI, angleDifference(vector0, vector180));
}

class UtilEigenGeometryPolygons : public ::testing::Test {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
    UtilEigenGeometryPolygons() = default;

    Eigen::Vector2d p0{0., 0.}; // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
    Eigen::Vector2d p1{1., 1.}; // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
    Eigen::Vector2d p2{2., 2.}; // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)

    polygon_t poly01{p0, p1};     // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
    polygon_t poly12{p1, p2};     // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
    polygon_t poly02{p0, p1, p2}; // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
};

TEST_F(UtilEigenGeometryPolygons, angleDifferencePoseToPolygon) { // NOLINT

    Eigen::Affine2d pose{Eigen::Affine2d::Identity()};
    pose.translation().setOnes();

    EXPECT_DOUBLE_EQ(M_PI / 4., angleDifference(pose, poly01));
    EXPECT_DOUBLE_EQ(M_PI / 4., angleDifference(pose, poly12));
    EXPECT_DOUBLE_EQ(M_PI / 4., angleDifference(pose, poly02));
}

TEST_F(UtilEigenGeometryPolygons, cosineSimilarityPoseToPolygon) { // NOLINT

    Eigen::Affine2d pose{Eigen::Affine2d::Identity()};
    pose.translation().setOnes();

    EXPECT_DOUBLE_EQ(cos(M_PI / 4.), cosineSimilarity(pose, poly01));
    EXPECT_DOUBLE_EQ(cos(M_PI / 4.), cosineSimilarity(pose, poly12));
    EXPECT_DOUBLE_EQ(cos(M_PI / 4.), cosineSimilarity(pose, poly02));
}


TEST_F(UtilEigenGeometryPolygons, getClosestId) { // NOLINT
    EXPECT_EQ(0ul, getClosestId(p0, poly02));
    EXPECT_EQ(2ul, getClosestId(Eigen::Vector2d{100., 100.}, poly02));
}

TEST_F(UtilEigenGeometryPolygons, lineStripOrientation) { // NOLINT

    Eigen::Affine2d pose{Eigen::Affine2d::Identity()};
    pose.translation().setOnes();

    EXPECT_DOUBLE_EQ(M_PI / 4., lineStripOrientation(p0, p1));
    EXPECT_DOUBLE_EQ(0., lineStripOrientation(p0, p0));
}

TEST_F(UtilEigenGeometryPolygons, splitPolygonRight) { // NOLINT
    polygon_t poly12new;
    splitPolygonRight(poly02, 0, poly12new);
    EXPECT_DOUBLE_EQ(poly12[0].x(), poly12new[0].x());
    EXPECT_DOUBLE_EQ(poly12[0].y(), poly12new[0].y());
}

TEST_F(UtilEigenGeometryPolygons, canSplitPolygonRight) { // NOLINT
    EXPECT_TRUE(canSplitPolygonRight(poly02, 1));
    EXPECT_FALSE(canSplitPolygonRight(poly02, 2));
}

TEST_F(UtilEigenGeometryPolygons, SamplePolygon) { // NOLINT
    polygon_t sampled = addIntermediateSamplesToPolygon(poly02, 0.1);
    ASSERT_GT(sampled.size(), 1ul);
    for (size_t i = 1; i < sampled.size(); ++i) {
        const auto dist = (sampled.at(i) - sampled.at(i - 1)).norm();
        EXPECT_LE(dist, 0.1);
        EXPECT_GE(dist, 0.05);
    }
}
