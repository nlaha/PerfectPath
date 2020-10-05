#pragma once
#include <glm/glm.hpp>

struct Ray {
	glm::vec3 o, d;
	Ray(
		glm::vec3 o0 = glm::vec3(0, 0, 0),
		glm::vec3 d0 = glm::vec3(0, 0, 0)
	) {
		o = o0, d = glm::normalize(d0);
	}
};