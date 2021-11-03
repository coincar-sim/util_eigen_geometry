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

#include "util_eigen_geometry.hpp"

#include <sstream>

namespace util_eigen_geometry {

size_t getClosestId(const Eigen::Vector2d& point, const polygon_t& polygon) {

    if (polygon.size() < 2) {
        throw std::invalid_argument("polygon must contain at least two points");
    }

    size_t closestId = 0;
    double smallestDist = (polygon[0] - point).norm();

    for (size_t i = 1; i < polygon.size(); i++) {
        double currentDist = (polygon[i] - point).norm();
        if (currentDist < smallestDist) {
            closestId = i;
            smallestDist = currentDist;
        }
    }
    return closestId;
}


double lineStripOrientation(const Eigen::Vector2d& point1, const Eigen::Vector2d& point2) {
    return normalizeAngleRadians(std::atan2(point2[1] - point1[1], point2[0] - point1[0]));
}

void splitPolygonRight(const polygon_t& inputPolygon, size_t id, polygon_t& outputPolygon) {
    outputPolygon.clear();
    if (id < inputPolygon.size() - 1) {
        outputPolygon.insert(outputPolygon.begin(), inputPolygon.begin() + id + 1, inputPolygon.end());
    } else {
        std::stringstream errorMsg;
        errorMsg << "id must be smaller than polygon.size()-1"
                 << " but inputPolygon.size()=" << inputPolygon.size() << ", id=" << id;
        throw std::invalid_argument(errorMsg.str());
    }
}


double positiveFloatModulo(double x, double y) {
    if (std::abs(y) < 10.e-9) {
        return std::abs(x);
    }
    // Result is always positive; not similar to fmod()
    return x - y * floor(x / y);
}

double normalizeAngleDegrees(double x) {
    return positiveFloatModulo((x + 180.0), 360.0) - 180.0;
}

double normalizeAngleRadians(double x) {
    return positiveFloatModulo((x + M_PI), 2.0 * M_PI) - M_PI;
}

double angleDifference(const double& targetAngle, const double& sourceAngle) {
    double angleDiff = targetAngle - sourceAngle;
    return normalizeAngleRadians(angleDiff);
}

double angleDifference(const Eigen::Vector2d& targetVector, const Eigen::Vector2d& sourceVector) {
    const auto source = sourceVector.normalized();
    const auto target = targetVector.normalized();
    const double angleError = std::atan2(source.x() * target.y() - source.y() * target.x(), source.dot(target));

    return normalizeAngleRadians(angleError);
}

double angleDifference(const Eigen::Affine2d& targetPose, const Eigen::Vector2d& sourceVector) {
    assert(isIsometry(targetPose));
    return angleDifference(Eigen::Isometry2d{targetPose.matrix()}, sourceVector);
}

double angleDifference(const Eigen::Isometry2d& targetPose, const Eigen::Vector2d& sourceVector) {
    return angleDifference(yawFromIsometry2d(targetPose), yawFromVector(sourceVector));
}

double angleDifference(const Eigen::Affine2d& pose, const polygon_t& polygon) {
    assert(isIsometry(pose));
    return angleDifference(Eigen::Isometry2d{pose.matrix()}, polygon);
}

double angleDifference(const Eigen::Isometry2d& pose, const polygon_t& polygon) {
    Eigen::Vector2d point = pose.translation();
    double pointAngle = yawFromIsometry2d(pose);
    size_t i = getClosestId(point, polygon);

    if (i == 0) {
        return angleDifference(lineStripOrientation(polygon[0], polygon[1]), pointAngle);
    }

    if (i == polygon.size() - 1) {
        return angleDifference(lineStripOrientation(polygon[i - 1], polygon[i]), pointAngle);
    }

    double angleDiff1 = angleDifference(lineStripOrientation(polygon[i - 1], polygon[i]), pointAngle);
    double angleDiff2 = angleDifference(lineStripOrientation(polygon[i], polygon[i + 1]), pointAngle);

    if (std::abs(angleDiff1) > std::abs(angleDiff2)) {
        return angleDiff2;
    }
    return angleDiff1;
}

double angleDifferenceDegrees(const double& targetAngle, const double& sourceAngle) {
    double angleDiff = targetAngle - sourceAngle;
    return normalizeAngleDegrees(angleDiff);
}

double cosineSimilarity(const double& angle1, const double& angle2) {
    return cos(angle1 - angle2);
}

double cosineSimilarity(const Eigen::Affine2d& pose, const polygon_t& polygon) {
    assert(isIsometry(pose));
    return cosineSimilarity(Eigen::Isometry2d{pose.matrix()}, polygon);
}

double cosineSimilarity(const Eigen::Isometry2d& pose, const polygon_t& polygon) {
    Eigen::Vector2d point = pose.translation();
    double pointAngle = yawFromIsometry2d(pose);
    size_t i = getClosestId(point, polygon);

    if (i == 0) {
        return cosineSimilarity(lineStripOrientation(polygon[0], polygon[1]), pointAngle);
    }

    if (i == polygon.size() - 1) {
        return cosineSimilarity(lineStripOrientation(polygon[i - 1], polygon[i]), pointAngle);
    }

    double cosSim1 = cosineSimilarity(lineStripOrientation(polygon[i - 1], polygon[i]), pointAngle);
    double cosSim2 = cosineSimilarity(lineStripOrientation(polygon[i], polygon[i + 1]), pointAngle);

    if (cosSim1 > cosSim2) {
        return cosSim1;
    }
    return cosSim2;
}

Eigen::Isometry2d isometry2dFromVector(const Eigen::Vector2d& orientationVector) {
    Eigen::Isometry2d isometry{Eigen::Isometry2d::Identity()};
    return isometry.rotate(yawFromVector(orientationVector));
}

Eigen::Vector2d vectorFromIsometry2d(const Eigen::Isometry2d& pose) {
    Eigen::Vector2d vector{Eigen::Vector2d::UnitX()};
    return pose.rotation() * vector;
}

Eigen::Affine2d affine2dFromXYOfAffine3d(const Eigen::Affine3d& pose) {
    assert(isIsometry(pose));
    return isometry2dFromXYOfIsometry3d(Eigen::Isometry3d{pose.matrix()});
}

Eigen::Isometry2d isometry2dFromXYOfIsometry3d(const Eigen::Isometry3d& pose) {
    Eigen::Isometry2d pose2d = Eigen::Isometry2d::Identity();
    pose2d.linear() = pose.linear().topLeftCorner<2, 2>();
    pose2d.translation() = pose.translation().topRows<2>();
    return pose2d;
}

double yawFromVector(const Eigen::Vector2d& orientationVector) {
    return normalizeAngleRadians(std::atan2(orientationVector.y(), orientationVector.x()));
}

double yawFromAffine2d(const Eigen::Affine2d& pose) {
    assert(isIsometry(pose));
    return yawFromIsometry2d(Eigen::Isometry2d{pose.matrix()});
}

double yawFromIsometry2d(const Eigen::Isometry2d& pose) {
    Eigen::Rotation2D<double> rot;
    rot.fromRotationMatrix(pose.linear());
    return rot.smallestAngle();
}

double yawFromAffine3d(const Eigen::Affine3d& pose) {
    assert(isIsometry(pose));
    return yawFromIsometry3d(Eigen::Isometry3d{pose.matrix()});
}

double yawFromIsometry3d(const Eigen::Isometry3d& pose) {
    return yawFromIsometry2d(isometry2dFromXYOfIsometry3d(pose));
}

bool canSplitPolygonRight(const polygon_t& inputPolygon, const size_t id) {
    return (id < inputPolygon.size() - 1);
}

polygon_t addIntermediateSamplesToPolygon(const polygon_t& polygon, const double maxSampleDist) {
    if (polygon.size() < 2) {
        return polygon;
    }

    polygon_t ret;
    for (size_t i = 1; i < polygon.size(); ++i) {
        const auto& cur = polygon.at(i);
        const auto& prev = polygon.at(i - 1);
        Eigen::Vector2d delta = cur - prev;
        auto length = delta.norm();
        auto nSegments = static_cast<int>(std::ceil(length / maxSampleDist));
        for (int j = 0; j < nSegments; ++j) {
            ret.emplace_back(prev + delta * (static_cast<double>(j) / static_cast<double>(nSegments)));
        }
    }
    ret.emplace_back(polygon.back());
    return ret;
}

Eigen::Affine3d affine3dXYFromAffine2d(const Eigen::Affine2d& pose) {
    assert(isIsometry(pose));
    return isometry3dXYFromIsometry2d(Eigen::Isometry2d{pose.matrix()});
}

Eigen::Isometry3d isometry3dXYFromIsometry2d(const Eigen::Isometry2d& pose) {
    Eigen::Isometry3d pose3d = Eigen::Isometry3d::Identity();
    pose3d.linear().topLeftCorner<2, 2>() = pose.linear();
    pose3d.translation().topRows<2>() = pose.translation();
    return pose3d;
}

} // namespace util_eigen_geometry
