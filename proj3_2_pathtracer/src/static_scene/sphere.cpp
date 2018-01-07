#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {
    
    bool Sphere::test(const Ray& r, double& t1, double& t2) const {
        
        
        double a = dot(r.d, r.d);
        double b = dot(2*(r.o-o), r.d);
        double c = dot(r.o-o, r.o-o) - r2;
        double check = b*b-4*a*c;
        if(check >= 0){
            double tmp1 = (-1*b + sqrt(check))/(2*a);
            double tmp2 = (-1*b - sqrt(check))/(2*a);
            t1 = min(tmp1, tmp2);
            t2 = max(tmp1, tmp2);
            if(t1 < 0){
                t1 = t2; //If t1 is negative, i.e origin of ray is in circle
            }
            if(t1 > r.min_t && t1 < r.max_t){
                r.max_t = t1;
                return true;
            }
            
            else if(t2 > r.min_t && t2 < r.max_t){
                r.max_t = t2;
                return true;
            }
        }
        return false;
        
        
    }
    
    bool Sphere::intersect(const Ray& r) const {
        
        // TODO Part 1, task 4:
        // Implement ray - sphere intersection.
        // Note that you might want to use the the Sphere::test helper here.
        double t1 = 0;
        double t2 = 0;
        bool result = test(r, t1, t2);
        
        return result;
        
    }
    
    bool Sphere::intersect(const Ray& r, Intersection *i) const {
        
        // TODO Part 1m task 4:
        // Implement ray - sphere intersection.
        // Note again that you might want to use the the Sphere::test helper here.
        // When an intersection takes place, the Intersection data should be updated
        // correspondingly.
        double t1 = 0;
        double t2 = 0;
        bool result = test(r, t1, t2);
        
        if(result == true){
            i->t = t1;
            // primitive
            Vector3D normal = r.o + t1*r.d - o;
            
            normal.normalize();
            i->n = normal;
            
            r.max_t = t1;
            i->primitive = this;
            // bsdf
            i->bsdf = this->get_bsdf();
            
        }
        
        return result;
        
        
        
    }
    
    void Sphere::draw(const Color& c) const {
        Misc::draw_sphere_opengl(o, r, c);
    }
    
    void Sphere::drawOutline(const Color& c) const {
        //Misc::draw_sphere_opengl(o, r, c);
    }
    
    
} // namespace StaticScene
} // namespace CGL
