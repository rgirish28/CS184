#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {

	// TODO (Part 3): Handle collisions with spheres.
	//

	double dist = (origin - pm.position).norm();

	if(dist <= radius){
		Vector3D point = pm.position + (radius - dist)*(pm.position - origin).unit();
		double distance = (pm.last_position - point).norm();
		Vector3D correction = distance*(point - pm.last_position).unit();
		pm.position = pm.last_position + (correction)*(friction);

	}
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  Misc::draw_sphere(shader, origin, radius * 0.92);
}
