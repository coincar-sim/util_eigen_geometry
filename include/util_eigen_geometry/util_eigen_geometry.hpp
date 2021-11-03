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

#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

namespace util_eigen_geometry {

// Cannot change names due to backward compatibility.
// NOLINTNEXTLINE(readability-identifier-naming)
using polygon_t = std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>;
// NOLINTNEXTLINE(readability-identifier-naming)
using polygon_ptr_t = std::shared_ptr<polygon_t>;

size_t getClosestId(const Eigen::Vector2d& point, const polygon_t& polygon);
double lineStripOrientation(const Eigen::Vector2d& point1, const Eigen::Vector2d& point2);
void splitPolygonRight(const polygon_t& inputPolygon, size_t id, polygon_t& outputPolygon);
bool canSplitPolygonRight(const polygon_t& inputPolygon, size_t id);
polygon_t addIntermediateSamplesToPolygon(const polygon_t& polygon, double maxSampleDist);
// void transformPolygon(const Eigen::Affine2d pose, const polygon_t& inputPolygon, polygon_t& outputPolygon);

double positiveFloatModulo(double x, double y);
double normalizeAngleDegrees(double x);
double normalizeAngleRadians(double x);

double angleDifference(const double& targetAngle, const double& sourceAngle);
double angleDifference(const Eigen::Vector2d& targetVector, const Eigen::Vector2d& sourceVector);
double angleDifference(const Eigen::Affine2d& pose, const polygon_t& polygon);
double angleDifference(const Eigen::Isometry2d& pose, const polygon_t& polygon);
double angleDifferenceDegrees(const double& targetAngle, const double& sourceAngle);

double cosineSimilarity(const double& angle1, const double& angle2);
double cosineSimilarity(const Eigen::Affine2d& pose, const polygon_t& polygon);
double cosineSimilarity(const Eigen::Isometry2d& pose, const polygon_t& polygon);

Eigen::Affine2d affine2dFromXYOfAffine3d(const Eigen::Affine3d& pose);
Eigen::Isometry2d isometry2dFromXYOfIsometry3d(const Eigen::Isometry3d& pose);
Eigen::Affine3d affine3dXYFromAffine2d(const Eigen::Affine2d& pose);
Eigen::Isometry3d isometry3dXYFromIsometry2d(const Eigen::Isometry2d& pose);

double yawFromVector(const Eigen::Vector2d& orientationVector);
double yawFromAffine2d(const Eigen::Affine2d& pose);
double yawFromIsometry2d(const Eigen::Isometry2d& pose);
double yawFromAffine3d(const Eigen::Affine3d& pose);
double yawFromIsometry3d(const Eigen::Isometry3d& pose);

template <typename Scalar, int Dim, int Options>
bool isIsometry(const Eigen::Transform<Scalar, Dim, Eigen::Affine, Options>& transformation);

} // namespace util_eigen_geometry

#include "internal/util_eigen_geometry_impl.hpp"
