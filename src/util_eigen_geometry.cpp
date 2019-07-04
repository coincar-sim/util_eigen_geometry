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


double angleDifference(const Eigen::Affine2d& pose, const polygon_t& polygon) {
    Eigen::Vector2d point = pose.translation();
    double pointAngle = yawFromAffine2d(pose);
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
    } else {
        return angleDiff1;
    }
}

double angleDifferenceDegrees(const double& targetAngle, const double& sourceAngle) {
    double angleDiff = targetAngle - sourceAngle;
    return normalizeAngleDegrees(angleDiff);
}

double cosineSimilarity(const double& angle1, const double& angle2) {
    return cos(angle1 - angle2);
}

double cosineSimilarity(const Eigen::Affine2d& pose, const polygon_t& polygon) {
    Eigen::Vector2d point = pose.translation();
    double pointAngle = yawFromAffine2d(pose);
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
    } else {
        return cosSim2;
    }
}

Eigen::Affine2d affine2dFromXYOfAffine3d(const Eigen::Affine3d& pose) {
    Eigen::Affine2d pose2d =
        Eigen::Translation2d(pose.translation().topRows<2>()) * pose.linear().topLeftCorner<2, 2>();
    return pose2d;
}

double yawFromAffine2d(const Eigen::Affine2d& pose) {
    Eigen::Rotation2D<double> rot;
    rot.fromRotationMatrix(pose.linear());
    return rot.smallestAngle();
}

double yawFromAffine3d(const Eigen::Affine3d& pose) {
    return yawFromAffine2d(affine2dFromXYOfAffine3d(pose));
}


} // namespace util_eigen_geometry
