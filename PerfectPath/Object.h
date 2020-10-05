#pragma once
#include <glm/glm.hpp>
#include <stdlib.h>
#include "Ray.h"

// base class for all 3D objects
class Object
{
public:
	glm::vec3 mat_color;
	double mat_emission;
	int mat_type;
	void setMat(
		glm::vec3 mat_color_ = glm::vec3(0, 0, 0),
		double mat_emission_ = 0,
		int mat_type_ = 0
	) {
		mat_color = mat_color_;
		mat_emission = mat_emission_;
		mat_type = mat_type_;
	}
	virtual double intersect(const Ray&) const = 0;
	virtual glm::vec3 normal(const glm::vec3&) const = 0;
};

