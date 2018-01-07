#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

    void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
        
        Vector3D z = Vector3D(n.x, n.y, n.z);
        Vector3D h = z;
        if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
        else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
        else h.z = 1.0;
        
        z.normalize();
        Vector3D y = cross(h, z);
        y.normalize();
        Vector3D x = cross(z, y);
        x.normalize();
        
        o2w[0] = x;
        o2w[1] = y;
        o2w[2] = z;
    }
    
    // Diffuse BSDF //
    
    Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
        return reflectance * (1.0 / PI);
    }
    
    Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
        // Generate random directions
        
        *wi = sampler.get_sample(pdf);
        
        return reflectance * (1.0 / PI);
    }
    
    // Mirror BSDF //
    
    Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
        //    double eps = 1e-3;
        //
        //    if (fabs(wo[2] - wi[2]) < eps && fabs(wo[0] + wi[0]) < eps && fabs(wo[1] + wi[1]) < eps ) {
        //        return reflectance * (1/std::max(wi[2],1e-8));
        //    }
        
        // Note: Because we already using emission for light source, zero f for direct lighting.
        
        return Spectrum();
    }
    
    Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
        
        // TODO:
        // Implement MirrorBSDF
        
        reflect(wo, wi);
        *pdf = 1;
        
        return reflectance * (1/std::max(wo[2],1e-8));
    }
    
    // Glossy BSDF //
    
    /*
     Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
     return Spectrum();
     }
     Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
     *pdf = 1.0f;
     return reflect(wo, wi, reflectance);
     }
     */
    
    // Refraction BSDF //
    
    Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
        return Spectrum();
    }
    
    Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
        
        // TODO:
        // Implement RefractionBSDF
        
        *pdf = 1;
        
        // Get the initial refract direction.
        bool res = refract(wo, wi, ior);
        if (!res) {
            return Spectrum();
        }
        
        double ni = ior;
        double no = 1;
        if (wo[2] < 0) {
            swap(ni,no);
        }
        
        double ratio = no/ni;
        return transmittance * ratio*ratio * (1/std::max(std::fabs((*wi)[2]),1e-8));
    }
    
    // Glass BSDF //
    
    Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
        
        return Spectrum();
    }
    
    Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
        
        // TODO:
        // Compute Fresnel coefficient and either reflect or refract based on it.
        
        *pdf = 1;
        
        // Get the initial refract direction.
        bool res = refract(wo, wi, ior);
        if (!res) {
            return reflectance * (1/std::max(std::fabs((*wi)[2]),1e-8));
        }
        
        
        // Calculating Fr
        double ni = ior;
        double no = 1;
        double cos_i = std::fabs((*wi)[2]);
        double cos_o = std::fabs(wo[2]);
        if (wo[2] < 0) {
            swap(ni,no);
        }
       
	double R0 = ((1.-ior)*(1.-ior))/((1.+ior)*(1.+ior));

	double R = R0 + (1-R0)*pow((1.-abs_cos_theta(*wi)),5);
       
        if (coin_flip(R)) {    // If we choose reflection
            reflect(wo, wi);
	    *pdf = R;
            // Here we don't need to multiply Fr because we already using randomized strategy to achieve it.
            return R* reflectance * (1/std::max(std::fabs((*wi)[2]),1e-8));
        }
        else{                       // If we choose refraction
        
		float ratio = ior;
		*pdf = 1.-R;
		if (wo[2] > 0) 
			ratio = 1/ratio;

            // Here we don't need to multiply (1-Fr) because we already using randomized strategy to achieve it.
            return (1.-R)*transmittance * ratio * ratio * (1/std::max(std::fabs((*wi)[2]),1e-8));
        }
        
    }
    
    void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
        
        // TODO:
        // Implement reflection of wo about normal (0,0,1) and store result in wi.
        *wi = Vector3D(-wo[0],-wo[1],wo[2]);
    }
    
    bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
        
        // TODO:
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.  
        
        int sign = 1;
        float ratio = ior;
        if (wo[2] > 0) {
            sign = -1;
            ratio = 1/ratio;
        }
        
        float cos2_wi = 1 - ratio*ratio*(1 - wo[2]*wo[2]);
        if (cos2_wi < 0) {
            *wi = Vector3D(-wo[0],-wo[1],wo[2]);
            return false;
        }
        
        *wi = Vector3D(-wo[0]*ratio,-wo[1]*ratio,sign * sqrt(cos2_wi)).unit();
        
        return true;
    }
    
    // Emission BSDF //
    
    Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
        return Spectrum();
    }
    
    Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
        *wi  = sampler.get_sample(pdf);
        return Spectrum();
    }
    
    
    Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
        // TODO: proj3-2, part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Spectrum.        
        
        double cos_i = std::max(wi[2],1e-8);

		
        Spectrum rs = (eta*eta + k*k - 2*eta*cos_i + cos_i*cos_i)/(eta*eta + k*k +2*eta*cos_i + cos_i*cos_i);
        Spectrum rp = ((eta*eta + k*k)*cos_i*cos_i -  2*eta*cos_i + 1) / ((eta*eta + k*k)*cos_i*cos_i +2*eta*cos_i + 1);
        
        Spectrum F = (rs + rp)/2.;

	if(isnan(F.r) || isnan(F.g) || isnan(F.b) || F.r < 0 || F.g < 0 || F.b < 0)
		std::cout<<"Here1"<<std::endl;


	
	F = Spectrum(clamp(F.r,0.,1.),clamp(F.r,0.,1.),clamp(F.r,0.,1.));


	return F;
    }
    
    double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
        
	    if(isnan(Lambda(wi)) || isnan(Lambda(wo)))
			    std::cout<<"Here2"<<std::endl;
	    

	    return clamp(1.0 / (1.0 + Lambda(wi) + Lambda(wo)),0.,1.);


    }
    
    double MicrofacetBSDF::D(const Vector3D& h) {
        // TODO: proj3-2, part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
      

        double theta = getTheta(h);
        
	if(isinf(tan(theta)) || isnan(cos(theta)))
		return 0;

	double tan_theta =std::max(std::fabs(tan(theta)),1e-8);
	double costheta =std::max(std::fabs(cos(theta)),1e-8);
	double expn = (-tan_theta*tan_theta)/(alpha*alpha);
	double den = PI*alpha*alpha*pow(costheta,4.);
	
        
	if(isnan(exp(expn)/den) || exp(expn)/den < 0 || den == 0)	
		std::cout<<"Here4"<<std::endl;

        return clamp(exp(expn)/den,0.,INFINITY);
    }
    
    Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
        // TODO: proj3-2, part 2
        // Implement microfacet model here.
     
	Vector3D n(0,0,1);

	if(dot(wo,n) <= 0 || dot(wi,n) <= 0){
//		std::cout<<"Here5"<<std::endl;
		return Spectrum();
	}
	
	Vector3D h = (wo+wi);
	h.normalize();

	return (F(wi)*G(wo,wi)*D(h))/(4.*std::max(std::fabs(wi[2]),1e-8)*std::max(std::fabs(wo[2]),1e-8));
        
    }
    
    Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
        // TODO: proj3-2, part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF values
	//
	//
/*    
	    Vector2D sample =  sampler.get_sample();
	    double theta = atan(sqrt(-alpha*alpha*log(1-sample.x)));
	    while((-alpha*alpha*log(1-sample.x))<0 || cos(theta) <= 0)
		    sample = sampler.get_sample();
	    theta = atan(sqrt(-alpha*alpha*log(1-sample.x)));
	    double phi = 2*PI*sample.y;

	    Vector3D h = Vector3D(cos(phi)*sin(theta),sin(phi)*sin(theta),std::max(fabs(cos(theta)),1e-8));    
			    
	    *wi = h - wo;
	    wi->normalize();

	    if(isinf(fabs(tan(theta))))
		    return Spectrum();

	    double expn = (-tan(theta)*tan(theta))/(alpha*alpha);
	    double den = alpha*alpha*pow(std::max(fabs(cos(theta)),1e-8),3);

	    double p_theta = (2*std::max(fabs(sin(theta)),1e-8) * fabs(exp(expn)))/den;
	    double p_phi = 1./(2*PI);
	    double p_h = (p_theta*p_phi)/std::max(fabs(sin(theta)),1e-8);
	    
	    
	    *pdf = p_h/(4*std::max(fabs(dot(*wi,h)),1e-8));
*/	    
	    *wi = cosineHemisphereSampler.get_sample(pdf);
	    if(*pdf <= 0)
		    std::cout<<"Here6"<<std::endl;

	    return MicrofacetBSDF::f(wo, *wi);
    }


} // namespace CGL
