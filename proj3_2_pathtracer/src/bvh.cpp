    #include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {
    
    BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                       size_t max_leaf_size) {
        
        root = construct_bvh(_primitives, max_leaf_size);
        
    }
    
    BVHAccel::~BVHAccel() {
        if (root) delete root;
    }
    
    BBox BVHAccel::get_bbox() const {
        return root->bb;
    }
    
    void BVHAccel::draw(BVHNode *node, const Color& c) const {
        if (node->isLeaf()) {
            for (Primitive *p : *(node->prims))
                p->draw(c);
        } else {
            draw(node->l, c);
            draw(node->r, c);
        }
    }
    
    void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
        if (node->isLeaf()) {
            for (Primitive *p : *(node->prims))
                p->drawOutline(c);
        } else {
            drawOutline(node->l, c);
            drawOutline(node->r, c);
        }
    }
    
    BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
        
        // TODO Part 2, task 1:
        // Construct a BVH from the given vector of primitives and maximum leaf
        // size configuration. The starter code build a BVH aggregate with a
        // single leaf node (which is also the root) that encloses all the
        // primitives.
        BBox centroid_box, bbox;
        
        for (Primitive *p : prims) {
            BBox bb = p->get_bbox();
            bbox.expand(bb);
            Vector3D c = bb.centroid();
            centroid_box.expand(c);
        }
        
        // You'll want to adjust this code.
        // Right now we just return a single node containing all primitives.
        BVHNode *node = new BVHNode(bbox);
        
    
        
        
        if(prims.size() > max_leaf_size){
            //Recurse left and right
            vector<Primitive *> left;
            vector<Primitive *> right;
            
            
            for (Primitive *p : prims){
                double extent[3] = {bbox.extent.x,bbox.extent.y,bbox.extent.z};
                int choice = std::distance(extent,std::max_element(extent,extent+3));
                
                Vector3D centroid_box_v = centroid_box.centroid();
                Vector3D centroid = p->get_bbox().centroid();
                
                
                if(choice == 0){
                    //Recurse on x
                    if(centroid.x < centroid_box_v.x)
                        left.push_back(p);
                    
                    else
                        right.push_back(p);
                    
                }
                else if (choice ==1){
                    //Recurse on y
                    if(centroid.y < centroid_box_v.y)
                        left.push_back(p);
                    
                    else
                        right.push_back(p);
                    
                }
                else if (choice == 2){
                    //Recurse on z
                    
                    if(centroid.z < centroid_box_v.z)
                        left.push_back(p);
                    else
                        right.push_back(p);
                    
                }
                
            }
            
            
            //If left is empty or right is empty, we are doing another iteration
            //Whic means we should empty out the left and right array
            
            if(left.empty() == true || right.empty() == true){
                node->prims = new vector<Primitive *>(prims);
                return node;
            }
            
            
            
            
            node->l = construct_bvh(left, max_leaf_size);
            node->r = construct_bvh(right, max_leaf_size);
            
        }
        
        
        node->prims = new vector<Primitive *>(prims);
        return node;
        
    }
    
    
    bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
        
        // TODO Part 2, task 3:
        // Implement BVH intersection.
        // Currently, we just naively loop over every primitive.
        
        double t0 = 0;
        double t1 = 0;
        if(node->bb.intersect(ray, t0, t1) == false){
            return false;
        }
        else{
            if(ray.min_t > t1 || ray.max_t < t0){
                return false;
            }
            
            if(node->isLeaf()){
                //
                for (Primitive *p : *(node->prims)) {
                    if (p->intersect(ray)){
                        return true;
                    }
                }
                return false;
            }
            else{
                bool hit1 = intersect(ray, node->l);
                bool hit2 = intersect(ray, node->r);
                return hit1 || hit2;
            }
        }
        
    }
    
    bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {
        
        // TODO Part 2, task 3:
        // Implement BVH intersection.
        // Currently, we just naively loop over every primitive.
        
        double t0 = 0;
        double t1 = 0;
        if(node->bb.intersect(ray, t0, t1) == false){
            return false;
        }
        
        else{
            if(ray.min_t > t1 || ray.max_t < t0){
                return false;
            }
            
            double t = ray.max_t + 1;
            Intersection isect;
            if(node->isLeaf()){
                //
                bool hit = false;
                for (Primitive *p : *(node->prims)) {
                    if (p->intersect(ray, &isect)){
                        if (isect.t < t){
                            t = isect.t;
                            *i = isect;
                        }
                        hit = true;
                    }
                }
                return hit;
            }
            else{
                bool hit1 = intersect(ray, i, node->l);
                bool hit2 = intersect(ray, i, node->r);
                return hit1 || hit2;
            }
        }
    }
    
}  // namespace StaticScene
}  // namespace CGL
