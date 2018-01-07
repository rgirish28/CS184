#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{
  void BezierCurve::evaluateStep()
  {
    // TODO Part 1.
    // Perform one step of the Bezier curve's evaluation at t using de Casteljau's algorithm for subdivision.
    // Store all of the intermediate control points into the 2D vector evaluatedLevels.

	  std::vector<Vector2D> current = evaluatedLevels.back();

	  if (current.size()>1){

		  std::vector<Vector2D> newLevel;

		  for( int i = 0; i < current.size() - 1; i++){
			  Vector2D value = (1-t)*current[i] + t*current[i+1];
			  newLevel.push_back(value);
		  }

		  evaluatedLevels.push_back(newLevel);
	  }



  }


  Vector3D BezierPatch::evaluate(double u, double v) const
  {
    // TODO Part 2.
    // Evaluate the Bezier surface at parameters (u, v) through 2D de Casteljau subdivision.
    // (i.e. Unlike Part 1 where we performed one subdivision level per call to evaluateStep, this function
    // should apply de Casteljau's algorithm until it computes the final, evaluated point on the surface)
    
	  std::vector<Vector3D> points;


	  for(int i = 0; i < 4; i++)
		  points.push_back(evaluate1D(controlPoints[i],u));

	  
    return evaluate1D(points,v);
  }

  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> points, double t) const
  {
    // TODO Part 2.
    // Optional helper function that you might find useful to implement as an abstraction when implementing BezierPatch::evaluate.
    // Given an array of 4 points that lie on a single curve, evaluates the Bezier curve at parameter t using 1D de Casteljau subdivision.


	  std::vector<Vector3D> current = points;
	  int step = 0;


	  while (step < 3){

		  std::vector<Vector3D> newLevel;

		  for( int i = 0; i < current.size() - 1; i++){
			  Vector3D value = (1-t)*current[i] + t*current[i+1];
			  newLevel.push_back(value);
		  }

		  current = newLevel;
		  step++;
	  }

    return current.back();
  }



  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // TODO Returns an approximate unit normal at this vertex, computed by
    // TODO taking the area-weighted average of the normals of neighboring
    // TODO triangles, then normalizing.

	  Vector3D n(0,0,0); // initialize a vector to store your normal sum
	  HalfedgeCIter h = halfedge(); // Since we're in a Vertex, this returns a halfedge
	                               // pointing _away_ from that vertex
				       //

	  do {
		  Vector3D a = this->position;
		  Vector3D b = h->twin()->vertex()->position;
		  Vector3D c = h->next()->next()->vertex()->position;
		//  n += cross(b-a,c-a);
		  n+= (cross(b-a,c-a)*cross(b-a,c-a).norm());
		  h = h->next()->next();
		  h = h->twin(); // Bump over to the halfedge pointing _toward_ the vertex.
				       //                // Now h->next() will be another edge on the same face,
				       //                               // sharing the central vertex.
	  }while (h != halfedge());

    return n.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // TODO This method should flip the given edge and return an iterator to the flipped edge.

	  if (e0->halfedge()->isBoundary())
		  return e0;

	  HalfedgeIter h0 = e0->halfedge();
	  HalfedgeIter h1 = h0->next();
	  HalfedgeIter h2 = h1->next();
	  HalfedgeIter h3 = h0->twin();
	  HalfedgeIter h4 = h3->next();
	  HalfedgeIter h5 = h4->next();

	  VertexIter v0 = h0->vertex();
	  VertexIter v1 = h3->vertex();
	  VertexIter v2 = h2->vertex();
	  VertexIter v3 = h5->vertex();

	  FaceIter f0 = h0->face();
	  FaceIter f1 = h3->face();	  
	  
	  h0->setNeighbors( h5 , h3, v2, e0, f0);
          h3->setNeighbors( h2 , h0, v3, e0, f1);
	  h5->setNeighbors( h1 , h5->twin(), v3, h5->edge(),f0);
	  h1->setNeighbors( h0 , h1->twin(), v1, h1->edge(),f0);
	  h2->setNeighbors( h4 , h2->twin(), v2, h2->edge(),f1);
	  h4->setNeighbors( h3 , h4->twin(), v0, h4->edge(),f1);

	  v0->halfedge() = h4;
//	  v1->halfedge() = h1;
	  v2->halfedge() = h0;
//	  v3->halfedge() = h3;

	  f0->halfedge() = h0;
	  f1->halfedge() = h3;


    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // TODO This method should split the given edge and return an iterator to the newly inserted vertex.
    // TODO The halfedge of this vertex should point along the edge that was split, rather than the new edges.

	  if (e0->halfedge()->isBoundary())
		  return e0->halfedge()->vertex();

	  HalfedgeIter h0 = e0->halfedge();
	  HalfedgeIter h1 = h0->next();
	  HalfedgeIter h2 = h1->next();
	  HalfedgeIter h3 = h0->twin();
	  HalfedgeIter h4 = h3->next();
	  HalfedgeIter h5 = h4->next();

	  VertexIter v0 = h0->vertex();
	  VertexIter v1 = h3->vertex();
	  VertexIter v2 = h2->vertex();
	  VertexIter v3 = h5->vertex();

	  FaceIter f0 = h0->face();	  
	 
	  HalfedgeIter hn0 = newHalfedge();
	  HalfedgeIter hn1 = newHalfedge();
	  HalfedgeIter hn2 = newHalfedge();
	  HalfedgeIter hn3 = newHalfedge();
	  HalfedgeIter hn4 = newHalfedge();
	  HalfedgeIter hn5 = newHalfedge();
	  EdgeIter e1 = newEdge();
	  EdgeIter e2 = newEdge();
	  EdgeIter e3 = newEdge();
	  FaceIter f1 = newFace();
	  FaceIter f2 = newFace();
	  FaceIter f3 = newFace();
	  VertexIter midpoint = newVertex();

	  midpoint->position = (v0->position + v1->position)/2.0;
	  midpoint->halfedge() = hn0;

	  h0->setNeighbors( hn2 , hn3, v0, e0, f0);
          h3->setNeighbors( hn4 , hn0, v1, e1, f3);
	  h5->setNeighbors( h3 , h5->twin(), v3, h5->edge(),f3);
	  h1->setNeighbors( hn1 , h1->twin(), v1, h1->edge(),f1);
	  h2->setNeighbors( h0 , h2->twin(), v2, h2->edge(),f0);
	  h4->setNeighbors( hn5 , h4->twin(), v0, h4->edge(),f2);
          hn0->setNeighbors( h1 , h3, midpoint, e1, f1);
          hn1->setNeighbors( hn0 , hn2, v2, e2, f1);
          hn2->setNeighbors( h2 , hn1, midpoint, e2, f0);
          hn3->setNeighbors( h4 , h0, midpoint, e0, f2);
          hn4->setNeighbors( h5 , hn5, midpoint, e3, f3);
          hn5->setNeighbors( hn3 , hn4, v3, e3, f2);

	  f0->halfedge() = h0;
	  f1->halfedge() = hn0;
	  f2->halfedge() = hn3;
	  f3->halfedge() = h3;

	  e1->halfedge() = hn0;
	  e2->halfedge() = hn1;
	  e3->halfedge() = hn4;

	  e0->isNew = false;
	  e1->isNew = false;
	  e2->isNew = true;
	  e3->isNew = true;


    return midpoint;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
    // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse. We will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

	for( VertexIter e = mesh.verticesBegin(); e != mesh.verticesEnd(); e++ )

    	{
		e->isNew = false;

		float n = e->degree();

		float u = 3./(8.*n);

		if (n == 3)
			u = 3./16.;

		Vector3D neighbor_position_sum(0,0,0);

		HalfedgeCIter h = e->halfedge();

		do{
			h = h->twin();
			neighbor_position_sum += h->vertex()->position;
			h = h->next();
		} while (h!= e->halfedge());

		e->newPosition = (1 - n*u) * e->position + u * neighbor_position_sum;


	}



	for( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ ){
	
		Vector3D v1 = e->halfedge()->vertex()->position;
		Vector3D v2 = e->halfedge()->twin()->vertex()->position;
		Vector3D v3 = e->halfedge()->next()->next()->vertex()->position;
		Vector3D v4 = e->halfedge()->twin()->next()->next()->vertex()->position;

		e->newPosition = 3./8. * (v1 + v2) + 1./8. * (v3 + v4);
		e->isNew = false;

	}


	int size = mesh.nEdges();
	int i = 0;

	for( EdgeIter e = mesh.edgesBegin(); i < size;  i++ ){

		if(e->isNew == false){
			VertexIter v = mesh.splitEdge(e);
			v->isNew = true;
			v->position = e->newPosition;
		}
		e++;

	}

	for( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ ){


		int count = 0;

		if (e->isNew){

			HalfedgeIter h = e->halfedge();

			if(h->vertex()->isNew)
		       		count++;
			if(h->twin()->vertex()->isNew)
				count++;
	
			if (count == 1){
				e = mesh.flipEdge(e);
				e->isNew = false;
						
			}
		}

	}

	    for(VertexIter v = mesh.verticesBegin();v!=mesh.verticesEnd();v++){
		    if(v->isNew == false)
			    v->position = v->newPosition;
		    else
			    v->isNew = false;
	    }
  // TODO Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // TODO and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // TODO a vertex of the original mesh.



    // TODO Next, compute the updated vertex positions associated with edges, and store it in Edge::newPosition.


    // TODO Next, we're going to split every edge in the mesh, in any order.  For future
    // TODO reference, we're also going to store some information about which subdivided
    // TODO edges come from splitting an edge in the original mesh, and which edges are new,
    // TODO by setting the flat Edge::isNew.  Note that in this loop, we only want to iterate
    // TODO over edges of the original mesh---otherwise, we'll end up splitting edges that we
    // TODO just split (and the loop will never end!)


    // TODO Now flip any new edge that connects an old and new vertex.


    // TODO Finally, copy the new vertex positions into final Vertex::position.
  }

}
