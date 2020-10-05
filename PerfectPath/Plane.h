#pragma once
#include <stdlib.h>
#include "Object.h"
#include "glm/glm.hpp"

class Plane :
    public Object
{
public:
    glm::vec3 n;
    double d;
	Plane(double d_ = 0, glm::vec3 n_ = glm::vec3(0, 0, 0)) {
		d = d_;
		n = n_;
	}
	double intersect(const Ray& ray) const {
		double d0 = glm::dot(n, ray.d);
		if (d0 != 0) {
			double t = -1 * (((glm::dot(n, ray.o)) + d) / d0);
			return (t > std::numeric_limits<float>::epsilon()) ? t : 0;
		}
		else return 0;
	}
	glm::vec3 normal(const glm::vec3& p0) const { return n; }
};

