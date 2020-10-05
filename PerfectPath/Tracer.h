#pragma once

#include <stdlib.h>
#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <ctime>

#include "Object.h"
#include "Plane.h"
#include "Sphere.h"

const double PI = 3.14159265358979323846;
typedef std::unordered_map<std::string, double> pl;

// Helpers for random number generation
std::mt19937 mersenneTwister;
std::uniform_real_distribution<double> uniform;

#define RND (2.0*uniform(mersenneTwister)-1.0)
#define RND2 (uniform(mersenneTwister))

class Intersection {
public:
	Intersection() { t = std::numeric_limits<double>::infinity(); object = nullptr; }
	Intersection(double t_, Object* object_) { t = t_; object = object_; }
	operator bool() { return object != nullptr; }
	double t;
	Object* object;
};

class Scene {
	std::vector<Object*> objects;

public:
	void add(Object* object) {
		objects.push_back(object);
	}

	Intersection intersect(const Ray& ray) const {
		Intersection closestIntersection;
		// intersect all objects, one after the other
		for (auto iter = objects.begin(); iter != objects.end(); ++iter) {
			double t = (*iter)->intersect(ray);
			if (t > std::numeric_limits<double>::epsilon() && t < closestIntersection.t) {
				closestIntersection.t = t;
				closestIntersection.object = *iter;
			}
		}
		return closestIntersection;
	}
};

glm::vec3 modVec(const glm::vec3& a, const glm::vec3& b) {
	return glm::vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}


void ons(const glm::vec3& v1, glm::vec3& v2, glm::vec3& v3) {
	if (std::abs(v1.x) > std::abs(v1.y)) {
		// project to the y = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.f / sqrtf(v1.x * v1.x + v1.z * v1.z);
		v2 = glm::vec3(-v1.z * invLen, 0.0f, v1.x * invLen);
	}
	else {
		// project to the x = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.0f / sqrtf(v1.y * v1.y + v1.z * v1.z);
		v2 = glm::vec3(0.0f, v1.z * invLen, -v1.y * invLen);
	}
	
	v3 = modVec(v1, v2);

}

// Class for generating the Halton low-discrepancy series for Quasi
// Monte Carlo integration.
class Halton {
	double value, inv_base;
public:
	void number(int i, int base) {
		double f = inv_base = 1.0 / base;
		value = 0.0;
		while (i > 0) {
			value += f * (double)(i % base);
			i /= base;
			f *= inv_base;
		}
	}
	void next() {
		double r = 1.0 - value - 0.0000001;
		if (inv_base < r) value += inv_base;
		else {
			double h = inv_base, hh;
			do { hh = h; h *= inv_base; } while (h >= r);
			value += hh + h - 1.0;
		}
	}
	double get() { return value; }
};

// Input is the pixel offset, output is the appropriate coordinate
// on the image plane
glm::vec3 camcr(const double x, const double y, const int width, const int height) {
	double w = width;
	double h = height;
	float fovx = PI / 4;
	float fovy = (h / w) * fovx;
	return glm::vec3(((2 * x - w) / w) * tan(fovx),
		-((2 * y - h) / h) * tan(fovy),
		-1.0);
}

// Uniform sampling on a hemisphere to produce outgoing ray directions.
// courtesy of http://www.rorydriscoll.com/2009/01/07/better-sampling/
glm::vec3 hemisphere(double u1, double u2) {
	const double r = sqrt(1.0 - u1 * u1);
	const double phi = 2 * PI * u2;
	return glm::vec3(cos(phi) * r, sin(phi) * r, u1);
}

void trace(Ray& ray, const Scene& scene, int depth, glm::vec3& clr, pl& params, Halton& hal, Halton& hal2) {

	// Russian roulette: starting at depth 5, each recursive step will stop with a probability of 0.1
	double rrFactor = 1.0;
	if (depth >= 5) {
		const double rrStopProbability = 0.1;
		if (RND2 <= rrStopProbability) {
			return;
		}
		rrFactor = 1.0 / (1.0 - rrStopProbability);
	}

	Intersection intersection = scene.intersect(ray);
	if (!intersection) return;

	// Travel the ray to the hit point where the closest object lies and compute the surface normal there.
	glm::vec3 hp = ray.o + ray.d * float(intersection.t);
	glm::vec3 N = intersection.object->normal(hp);
	ray.o = hp;
	// Add the emission, the L_e(x,w) part of the rendering equation, but scale it with the Russian Roulette
	// probability weight.
	const double emission = intersection.object->mat_emission;
	clr = clr + glm::vec3(emission, emission, emission) * float(rrFactor);

	// Diffuse BRDF - choose an outgoing direction with hemisphere sampling.
	if (intersection.object->mat_type == 1) {
		glm::vec3 rotX = glm::vec3(0, 0, 0);
		glm::vec3 rotY = glm::vec3(0, 0, 0);

		ons(N, rotX, rotY);
		glm::vec3 sampledDir = hemisphere(RND2, RND2);
		glm::vec3 rotatedDir;
		rotatedDir.x = glm::dot(glm::vec3(rotX.x, rotY.x, N.x), sampledDir);
		rotatedDir.y = glm::dot(glm::vec3(rotX.y, rotY.y, N.y), sampledDir);
		rotatedDir.z = glm::dot(glm::vec3(rotX.z, rotY.z, N.z), sampledDir);
		ray.d = rotatedDir;	// already normalized
		float cost = glm::dot(ray.d, N);
		glm::vec3 tmp = glm::vec3(0, 0, 0);
		trace(ray, scene, depth + 1, tmp, params, hal, hal2);
		clr = clr + (tmp * intersection.object->mat_color) * cost * 0.1f * float(rrFactor);
	}

	// Specular BRDF - this is a singularity in the rendering equation that follows
	// delta distribution, therefore we handle this case explicitly - one incoming
	// direction -> one outgoing direction, that is, the perfect reflection direction.
	if (intersection.object->mat_type == 2) {
		double cost = glm::dot(ray.d, N);
		ray.d = glm::normalize(ray.d - N * float(cost * 2));
		glm::vec3 tmp = glm::vec3(0, 0, 0);
		trace(ray, scene, depth + 1, tmp, params, hal, hal2);
		clr = clr + tmp * float(rrFactor);
	}

	// Glass/refractive BRDF - we use the vector version of Snell's law and Fresnel's law
	// to compute the outgoing reflection and refraction directions and probability weights.
	if (intersection.object->mat_type == 3) {
		double n = params["refr_index"];
		double R0 = (1.0 - n) / (1.0 + n);
		R0 = R0 * R0;
		if (glm::dot(N, ray.d) > 0) { // we're inside the medium
			N = N * -1.0f;
			n = 1 / n;
		}
		n = 1 / n;
		double cost1 = (glm::dot(N, ray.d)) * -1.0f; // cosine of theta_1
		double cost2 = 1.0 - n * n * (1.0 - cost1 * cost1); // cosine of theta_2
		double Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0); // Schlick-approximation
		if (cost2 > 0 && RND2 > Rprob) { // refraction direction
			ray.d = glm::normalize((ray.d * float(n)) + (N * float((n * cost1 - sqrt(cost2)))));
		}
		else { // reflection direction
			ray.d = glm::normalize(ray.d + N * float(cost1 * 2));
		}
		glm::vec3 tmp = glm::vec3(0, 0, 0);
		trace(ray, scene, depth + 1, tmp, params, hal, hal2);
		clr = clr + tmp * 1.15f * float(rrFactor);
	}
}

glm::vec3** render(const int width, const int height) {
	srand(time(NULL));
	pl params;
	Scene scene;
	auto add = [&scene](Object* s, glm::vec3 cl, double emission, int type) {
		s->setMat(cl, emission, type);
		scene.add(s);
	};

	// Radius, position, color, emission, type (1=diff, 2=spec, 3=refr) for spheres
	add(new Sphere(1.05, glm::vec3(-0.75, -1.45, -4.4)), glm::vec3(4, 8, 4), 0, 2); // Middle sphere
//	add(new Sphere(0.45,Vec(0.8,-2.05,-3.7)),Vec(10,10,1),0,3); // Right sphere
	add(new Sphere(0.5, glm::vec3(2.0, -2.05, -3.7)), glm::vec3(10, 10, 1), 0, 3); // Right sphere
	add(new Sphere(0.6, glm::vec3(-1.75, -1.95, -3.1)), glm::vec3(4, 4, 12), 0, 1); // Left sphere
	// Position, normal, color, emission, type for planes
	add(new Plane(2.5, glm::vec3(0, 1, 0)), glm::vec3(6, 6, 6), 0, 1); // Bottom plane
	add(new Plane(5.5, glm::vec3(0, 0, 1)), glm::vec3(6, 6, 6), 0, 1); // Back plane
	add(new Plane(2.75, glm::vec3(1, 0, 0)), glm::vec3(10, 2, 2), 0, 1); // Left plane
	add(new Plane(2.75, glm::vec3(-1, 0, 0)), glm::vec3(2, 10, 2), 0, 1); // Right plane
	add(new Plane(3.0, glm::vec3(0, -1, 0)), glm::vec3(6, 6, 6), 0, 1); // Ceiling plane
	add(new Plane(0.5, glm::vec3(0, 0, -1)), glm::vec3(6, 6, 6), 0, 1); // Front plane
	add(new Sphere(0.5, glm::vec3(0, 1.9, -3)), glm::vec3(0, 0, 0), 10000, 1); // Light

	params["refr_index"] = 1.5;
	params["spp"] = 108.0; // samples per pixel

	glm::vec3** pix = new glm::vec3 * [width];
	for (int i = 0; i < width; i++) {
		pix[i] = new glm::vec3[height];
	}

	const float spp = params["spp"];
	// correlated Halton-sequence dimensions
	Halton hal, hal2;
	hal.number(0, 2);
	hal2.number(0, 2);

	#pragma omp parallel for schedule(dynamic) firstprivate(hal,hal2)
	for (int col = 0; col < width; col++) {
		fprintf(stdout, "\rRendering: %1.0fspp %8.2f%%", spp, (double)col / width * 100);
		for (int row = 0; row < height; row++) {
			for (int s = 0; s < spp; s++) {
				glm::vec3 color = glm::vec3(0, 0, 0);
				Ray ray;
				ray.o = (glm::vec3(0, 0, 0)); // rays start out from here
				glm::vec3 cam = camcr(col, row, width, height); // construct image plane coordinates
				cam.x = cam.x + RND / 700; // anti-aliasing for free
				cam.y = cam.y + RND / 700;
				ray.d = glm::normalize(cam - ray.o); // point from the origin to the camera plane
				trace(ray, scene, 0, color, params, hal, hal2);
				//std::cout << color.x;
				pix[col][row] = pix[col][row] + color / spp; // write the contributions
			}
		}
	}

	return pix;
}