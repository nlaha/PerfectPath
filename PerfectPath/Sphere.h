#pragma once
#include <stdlib.h>
#include "Object.h"
#include "Ray.h"
#include "glm/glm.hpp"

class Sphere :
    public Object
{
public: 
    glm::vec3 center;
    double radius;

	Sphere(double radius_ = 0, glm::vec3 center_ = glm::vec3(0, 0, 0)) {center = center_; radius = radius_; }
	double intersect(const Ray& ray) const {
		double b = glm::dot(((ray.o - center) * 2.0f), ray.d);
		double c_ = glm::dot((ray.o - center), (ray.o - center)) - (radius * radius);
		double disc = b * b - 4 * c_;
		if (disc < 0) return 0;
		else disc = sqrt(disc);
		double sol1 = -b + disc;
		double sol2 = -b - disc;
		return (sol2 > std::numeric_limits<float>::epsilon()) ? sol2 / 2 : ((sol1 > std::numeric_limits<float>::epsilon()) ? sol1 / 2 : 0);
	}

	glm::vec3 normal(const glm::vec3& p0) const {
		return glm::normalize((p0 - center));
	}
};

