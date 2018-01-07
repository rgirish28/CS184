#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}
double fRand(double fMin, double fMax)
{
    	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

void Cloth::buildGrid() {
  // TODO (Part 1): Build a grid of masses and springs.
  //
  //


	for (int j=0;j<num_height_points;j++){
		for(int i=0;i<num_width_points;i++){
			double x = (double) i * ((double) width/ (double) num_width_points);
			double y = (double) j * ((double) height/ (double) num_height_points);

			Vector3D location;

			if (orientation == 0)
				location = Vector3D(x,1.0,y);
			else
				location = Vector3D(x,y, fRand(-1.0/1000.,1.0/1000.));


			bool pin = false;
			std::vector<int> index = {i,j};
			
			if (std::find(pinned.begin(), pinned.end(),index) != pinned.end())
				pin = true;

			point_masses.push_back(PointMass(location,pin));		
			
		}
	}
	for (int j=0;j<num_height_points;j++){
		for(int i=0;i<num_width_points;i++){			
			
			if(i>0)
				springs.push_back(Spring(&point_masses[j*num_width_points + i - 1], &point_masses[j*num_width_points + i], STRUCTURAL));
			if(j>0)
				springs.push_back(Spring(&point_masses[(j-1)*num_width_points + i], &point_masses[j*num_width_points + i], STRUCTURAL));
		
			if(j>0 && i>0)
				springs.push_back(Spring(&point_masses[(j-1)*num_width_points + i-1], &point_masses[j*num_width_points + i], SHEARING));

			if(j>0&&(i+1)<num_width_points)
				springs.push_back(Spring(&point_masses[(j-1)*num_width_points + i+1], &point_masses[j*num_width_points + i], SHEARING));
	
			if((j-2)>=0)
				springs.push_back(Spring(&point_masses[(j-2)*num_width_points + i], &point_masses[j*num_width_points + i], BENDING));

			if((i+2)<num_width_points)
				springs.push_back(Spring(&point_masses[j*num_width_points + i + 2], &point_masses[j*num_width_points + i], BENDING));
		}
	}
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.

  for (PointMass &pm : point_masses){
 
	  pm.forces = Vector3D(0,0,0);
	  for(Vector3D external : external_accelerations)
		  pm.forces += mass * external;
  }

  for (Spring &s : springs){

	  if(s.spring_type == STRUCTURAL && !cp->enable_structural_constraints)
		  continue;
	  if(s.spring_type == SHEARING && !cp->enable_shearing_constraints)
		  continue;
	  if(s.spring_type == BENDING && !cp->enable_bending_constraints)
		  continue;
	
	  Vector3D pa = s.pm_a->position;
	  Vector3D pb = s.pm_b->position;

	  double F = cp->ks * ( (pa-pb).norm() - s.rest_length);

	  Vector3D F_vec = F *(pa-pb).unit();
	  Vector3D F_opp = F *(pb-pa).unit();

	  s.pm_b->forces+= F_vec;
	  s.pm_a->forces+= F_opp;

  }
	

  // TODO (Part 2): Use Verlet integration to compute new point mass positions

  for (PointMass &pm : point_masses){
 
	  if(pm.pinned)
		  continue;

	  Vector3D acceleration = Vector3D(pm.forces/mass);
	  Vector3D last = Vector3D(pm.last_position);
	  pm.last_position = Vector3D(pm.position);
	  pm.position = pm.last_position + (1.-(cp->damping/100.0)) * (pm.last_position - last) + acceleration * delta_t*delta_t;

}



  // TODO (Part 4): Handle self-collisions.
  // This won't do anything until you complete Part 4.
  build_spatial_map();
  for (PointMass &pm : point_masses) {
    self_collide(pm, simulation_steps);
  }


  // TODO (Part 3): Handle collisions with other primitives.
  // This won't do anything until you complete Part 3.
  for (PointMass &pm : point_masses) {
    for (CollisionObject *co : *collision_objects) {
      co->collide(pm);
    }
  }


  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  

  for (Spring &s : springs){
	  double length = (s.pm_a->position - s.pm_b->position).norm();
	  double last_length = (s.pm_a->last_position - s.pm_b->last_position).norm();
	  
	  double percentage_inc = (fabs((length - last_length)) * 100.0)/last_length;

	  if(percentage_inc > 10.0){
		  
		  double length_dec = (length - 1.1*last_length);
		  if(s.pm_a->pinned)
			  s.pm_b->position = s.pm_b->position + length_dec*(s.pm_b->position - s.pm_a->position).unit();
		  else if(s.pm_b->pinned)
			  s.pm_a->position = s.pm_a->position + length_dec*(s.pm_b->position - s.pm_a->position).unit();
		  else if(!s.pm_a->pinned && s.pm_b->pinned){
			  s.pm_a->position = s.pm_a->position + 0.5*length_dec*(s.pm_b->position - s.pm_a->position).unit();
			  s.pm_b->position = s.pm_b->position + 0.5*length_dec*(s.pm_a->position - s.pm_b->position).unit();
		  }
		  s.rest_length = last_length+length_dec;
	  }
  }

}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.

}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.

}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents
  // membership in some uniquely identified 3D box volume.

  return 0.f;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm, pm + num_width_points, pm + 1));
      triangles.push_back(new Triangle(pm + 1, pm + num_width_points,
                                       pm + num_width_points + 1));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
