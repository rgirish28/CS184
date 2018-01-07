#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {
    
    bool BBox::intersect(const Ray& r, double& t0, double& t1) const {
        
        // TODO Part 2, task 2:
        // Implement ray - bounding box intersection test
        // If the ray intersected the bounding box within the range given by
        // t0, t1, update t0 and t1 with the new intersection times.
        float t_xmin = std::min((min.x-r.o.x)/r.d.x, (max.x-r.o.x)/r.d.x);
        float t_xmax = std::max((min.x-r.o.x)/r.d.x, (max.x-r.o.x)/r.d.x);
        
        float t_ymin = std::min((min.y-r.o.y)/r.d.y, (max.y-r.o.y)/r.d.y);
        float t_ymax = std::max((min.y-r.o.y)/r.d.y, (max.y-r.o.y)/r.d.y);
        
        float t_zmin = std::min((min.z-r.o.z)/r.d.z, (max.z-r.o.z)/r.d.z);
        float t_zmax = std::max((min.z-r.o.z)/r.d.z, (max.z-r.o.z)/r.d.z);
        
        
        if(t_xmin > t_ymax || t_ymin > t_xmax){ //Line misses completely
            return false;
        }
        if(t_xmin > t_ymin){
            t0 = t_xmin;
        }
        else{
            t0 = t_ymin;
        }
        if(t_xmax < t_ymax){
            t1 = t_xmax;
        }
        else{
            t1 = t_ymax;
        }
        
        if(t0 > t_zmax || t_zmin > t1){
            return false;
        }
        
        if(t0 < t_zmin){
            t0 = t_zmin;
        }
        if(t1 > t_zmax){
            t1 = t_zmax;
        }
        
        return true;
    }
    
    void BBox::draw(Color c) const {
        
        glColor4f(c.r, c.g, c.b, c.a);
        
        // top
        glBegin(GL_LINE_STRIP);
        glVertex3d(max.x, max.y, max.z);
        glVertex3d(max.x, max.y, min.z);
        glVertex3d(min.x, max.y, min.z);
        glVertex3d(min.x, max.y, max.z);
        glVertex3d(max.x, max.y, max.z);
        glEnd();
        
        // bottom
        glBegin(GL_LINE_STRIP);
        glVertex3d(min.x, min.y, min.z);
        glVertex3d(min.x, min.y, max.z);
        glVertex3d(max.x, min.y, max.z);
        glVertex3d(max.x, min.y, min.z);
        glVertex3d(min.x, min.y, min.z);
        glEnd();
        
        // side
        glBegin(GL_LINES);
        glVertex3d(max.x, max.y, max.z);
        glVertex3d(max.x, min.y, max.z);
        glVertex3d(max.x, max.y, min.z);
        glVertex3d(max.x, min.y, min.z);
        glVertex3d(min.x, max.y, min.z);
        glVertex3d(min.x, min.y, min.z);
        glVertex3d(min.x, max.y, max.z);
        glVertex3d(min.x, min.y, max.z);
        glEnd();
        
    }
    
    std::ostream& operator<<(std::ostream& os, const BBox& b) {
        return os << "BBOX(" << b.min << ", " << b.max << ")";
    }
    
} // namespace CGL
