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
    *wi = sampler.get_sample(pdf);
    return reflectance * (1.0 / PI);
}


// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 3-2 Part 1 Task 2
  // Implement MirrorBSDF
    
    reflect(wo, wi);
    *pdf = 1;
    
    return reflectance * (1/std::max(fabs(wo[2]),1e-8));
}


// Microfacet BSDF //

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Compute Fresnel term for reflection on dielectric-conductor interface.
    // You will need both eta and etaK, both of which are Spectrum.

	double cos_i = std::max(fabs(wi.z),1e-8);

	Spectrum Rs = ( (eta*eta + k*k) - 2*eta*cos_i + cos_i*cos_i  )/ (   (eta*eta + k*k) + 2*eta*cos_i + cos_i*cos_i  );
	Spectrum Rp = ( (eta*eta + k*k)*cos_i*cos_i - 2*eta*cos_i + 1 )/ (   (eta*eta + k*k)*cos_i*cos_i + 2*eta*cos_i + 1 );

	Spectrum F = (Rs + Rp)/2.;

    return F;
}

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
    // TODO: proj3-2, part 2
    // Compute Beckmann normal distribution function (NDF) here.
    // You will need the roughness alpha.
    
	
	double NdotH = std::max(fabs(h.z),1e-8);

	double r1 = 1./(PI * alpha *alpha *pow(NdotH,4.0));
	double r2 = (NdotH * NdotH - 1.0) / (alpha *alpha * NdotH * NdotH);


	
	return r1 * exp(r2);
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Implement microfacet model here.

	Vector3D n(0,0,1);

	if(wo.z <= 0 || wi.z <= 0)
		return Spectrum();

	Vector3D h = (wi + wo);
	h.normalize();



	return (F(wi)*G(wo,wi)*D(h))/(4.*std::max(fabs(wo.z),1e-8)*std::max(fabs(wi.z),1e-8));


}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    // TODO: proj3-2, part 2
    // *Importance* sample Beckmann normal distribution function (NDF) here.
    // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
    //       and return the sampled BRDF value.

	Vector2D sample;
	double tan_x = -1.;
	double cosx = -1;
	while(tan_x < 0 || cosx <= 0 ){
		sample = sampler.get_sample();
		tan_x = -alpha*alpha*log(1.0-sample.x);
		cosx = cos(atan(sqrt(tan_x)));
	}
	double theta = atan(sqrt(tan_x));
	double phi = 2.0*PI*sample.y;

	Vector3D h = Vector3D(sin(theta)*cos(phi), sin(theta)*sin(phi),cos(theta));

	h.normalize();

	*wi = 2*dot(wo,h)*h.unit() - wo;
	wi->normalize();
	
	double expn = (cos_theta(h) * cos_theta(h) - 1.0) / (alpha *alpha * cos_theta(h)  * cos_theta(h)); 
	double p_theta = (2.0*fabs(sin_theta(h)) * exp(expn))/(alpha*alpha*pow(cos_theta(h) ,3.0));
	double p_phi = 1./(2.0*PI);

	double p = (p_theta*p_phi)/fabs(sin_theta(h));

	if(dot(*wi,h) <= 0){
		*pdf = 0;
		return Spectrum();
	}

	*pdf = p/(4.0*fabs(dot(*wi,h)));


    return MicrofacetBSDF::f(wo, *wi);
}


// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 3-2 Part 1 Task 4
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
        
		double ratio = ior;
		*pdf = 1.-R;
		if (wo[2] > 0) 
			ratio = 1/ratio;

            // Here we don't need to multiply (1-Fr) because we already using randomized strategy to achieve it.
            return (1.-R)*transmittance * ratio * ratio * (1/std::max(std::fabs((*wi)[2]),1e-8));
        }   
      
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 3-2 Part 1 Task 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0],-wo[1],wo[2]);

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 3-2 Part 1 Task 3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When wo.z is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    
    int sign = 1;
    double ratio = ior;
    if (wo[2] > 0) {
        sign = -1;
        ratio = 1/ratio;
    }
    
    double cos2_wi = 1 - ratio*ratio*(1 - wo[2]*wo[2]);
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
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
