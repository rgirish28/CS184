#include "texture.h"
#include "CGL/color.h"

namespace CGL {

// Examines the enum parameters in sp and performs
// the appropriate sampling using the three helper functions below.
Color Texture::sample(const SampleParams &sp) {
  // Part 6: Fill in the functionality for sampling 
  //          nearest or bilinear in mipmap level 0, conditional on sp.psm
  // Part 7: Fill in full sampling (including trilinear), 
  //          conditional on sp.psm and sp.lsm


	if(sp.lsm == L_ZERO) {
		if (sp.psm == P_NEAREST)
			return sample_nearest(sp.uv, 0);
		else if (sp.psm == P_LINEAR)
			return sample_bilinear(sp.uv, 0);
	}
	else if(sp.lsm == L_NEAREST) {
		int level = 0;
	       	if (floor(get_level(sp))>0)
			level = floor(get_level(sp));

		if (sp.psm == P_NEAREST)
			return sample_nearest(sp.uv, level);
		else if (sp.psm == P_LINEAR)
			return sample_bilinear(sp.uv, level);
	}
	else if(sp.lsm == L_LINEAR) {

		int top = ceil(get_level(sp));
		int bottom = floor(get_level(sp));
		float diff = get_level(sp) - bottom;


		if(top < 0 || bottom < 0 ){
			top = 0;
			bottom = 0;
			diff = 0.5;
		}


		if (sp.psm == P_NEAREST)
			return sample_nearest(sp.uv, bottom)*(1.-diff)+sample_nearest(sp.uv, top)*diff;
		else if (sp.psm == P_LINEAR)
			return sample_bilinear(sp.uv, bottom)*(1.-diff)+sample_bilinear(sp.uv, top)*diff;
	}


  return Color();
}


inline void uint8_to_float(float dst[4], unsigned char *src) {
  uint8_t *src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8(unsigned char *dst, float src[4]) {
  uint8_t *dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
  dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
  dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  dst_uint8[3] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[3])));
}


// Given sp.du and sp.dv, returns the appropriate mipmap
// level to use for L_NEAREST or L_LINEAR filtering.
float Texture::get_level(const SampleParams &sp) {
  // Part 7: Fill this in.
  //
  // 
	float du_dx = sp.du.x*mipmap[0].width;
	float du_dy = sp.du.y*mipmap[0].height;

	float dv_dx = sp.dv.x*mipmap[0].width;
	float dv_dy = sp.dv.y*mipmap[0].height;
	
	float L = max(sqrt((du_dx * du_dx)+(dv_dx * dv_dx)),sqrt((du_dy * du_dy)+(dv_dy * dv_dy)));

	return log2f(L);
}

// Indexes into the level'th mipmap
// and returns the nearest pixel to (u,v)
Color Texture::sample_nearest(Vector2D uv, int level) {
  // Part 6: Fill this in.

	int i = (int) floor(uv.x*mipmap[level].width);
	int j = (int) floor(uv.y*mipmap[level].height);

	unsigned char *p = &mipmap[level].texels[0] + 4 * (j*mipmap[level].width + i);
	
	float dst[4];
	
	uint8_to_float(dst,p);

	return Color(dst[0],dst[1],dst[2],dst[3]);

}

// Indexes into the level'th mipmap
// and returns a bilinearly weighted combination of
// the four pixels surrounding (u,v)
Color Texture::sample_bilinear(Vector2D uv, int level) {
  // Part 6: Fill this in.


	int i = floor(uv.x*mipmap[level].width);
	int j = floor(uv.y*mipmap[level].height);


 	float u = (float) uv.x*mipmap[level].width - floor(uv.x*mipmap[level].width);	
 	float v = (float) uv.y*mipmap[level].height - floor(uv.y*mipmap[level].height);

	unsigned char *p1 = &mipmap[level].texels[0] + 4 * (j*mipmap[level].width + i);
 	unsigned char *q1 = &mipmap[level].texels[0] + 4 * (j*mipmap[level].width + i+1);
	unsigned char *r1 = &mipmap[level].texels[0] + 4 * ((j+1)*mipmap[level].width + i);
	unsigned char *s1 = &mipmap[level].texels[0] + 4 * ((j+1)*mipmap[level].width + i+1);

	float p[4],q[4],r[4],s[4];

        uint8_to_float(p,p1);	
        uint8_to_float(q,q1);	
        uint8_to_float(r,r1);	
        uint8_to_float(s,s1);	

	float finalColor[4];

	for(int i=0;i<4;i++)
		finalColor[i] = (1-u)*(1-v)*p[i]+u*(1-v)*q[i]+(1-u)*v*r[i]+u*v*s[i];

	return Color(finalColor[0],finalColor[1],finalColor[2],finalColor[3]);
}



/****************************************************************************/
void Texture::generate_mips(int startLevel) {

  // make sure there's a valid texture
  if (startLevel >= mipmap.size()) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth = mipmap[startLevel].width;
  int baseHeight = mipmap[startLevel].height;
  int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  mipmap.resize(startLevel + numSubLevels + 1);

  int width = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel &level = mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width = max(1, width / 2);
    //assert (width > 0);
    height = max(1, height / 2);
    //assert (height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);
  }

  // create mips
  int subLevels = numSubLevels - (startLevel + 1);
  for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
       mipLevel++) {

    MipLevel &prevLevel = mipmap[mipLevel - 1];
    MipLevel &currLevel = mipmap[mipLevel];

    int prevLevelPitch = prevLevel.width * 4; // 32 bit RGBA
    int currLevelPitch = currLevel.width * 4; // 32 bit RGBA

    unsigned char *prevLevelMem;
    unsigned char *currLevelMem;

    currLevelMem = (unsigned char *)&currLevel.texels[0];
    prevLevelMem = (unsigned char *)&prevLevel.texels[0];

    float wDecimal, wNorm, wWeight[3];
    int wSupport;
    float hDecimal, hNorm, hWeight[3];
    int hSupport;

    float result[4];
    float input[4];

    // conditional differentiates no rounding case from round down case
    if (prevLevel.width & 1) {
      wSupport = 3;
      wDecimal = 1.0f / (float)currLevel.width;
    } else {
      wSupport = 2;
      wDecimal = 0.0f;
    }

    // conditional differentiates no rounding case from round down case
    if (prevLevel.height & 1) {
      hSupport = 3;
      hDecimal = 1.0f / (float)currLevel.height;
    } else {
      hSupport = 2;
      hDecimal = 0.0f;
    }

    wNorm = 1.0f / (2.0f + wDecimal);
    hNorm = 1.0f / (2.0f + hDecimal);

    // case 1: reduction only in horizontal size (vertical size is 1)
    if (currLevel.height == prevLevel.height) {
      //assert (currLevel.height == 1);

      for (int i = 0; i < currLevel.width; i++) {
        wWeight[0] = wNorm * (1.0f - wDecimal * i);
        wWeight[1] = wNorm * 1.0f;
        wWeight[2] = wNorm * wDecimal * (i + 1);

        result[0] = result[1] = result[2] = result[3] = 0.0f;

        for (int ii = 0; ii < wSupport; ii++) {
          uint8_to_float(input, prevLevelMem + 4 * (2 * i + ii));
          result[0] += wWeight[ii] * input[0];
          result[1] += wWeight[ii] * input[1];
          result[2] += wWeight[ii] * input[2];
          result[3] += wWeight[ii] * input[3];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (4 * i), result);
      }

      // case 2: reduction only in vertical size (horizontal size is 1)
    } else if (currLevel.width == prevLevel.width) {
      //assert (currLevel.width == 1);

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        result[0] = result[1] = result[2] = result[3] = 0.0f;
        for (int jj = 0; jj < hSupport; jj++) {
          uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
          result[0] += hWeight[jj] * input[0];
          result[1] += hWeight[jj] * input[1];
          result[2] += hWeight[jj] * input[2];
          result[3] += hWeight[jj] * input[3];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (currLevelPitch * j), result);
      }

      // case 3: reduction in both horizontal and vertical size
    } else {

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = result[3] = 0.0f;

          // convolve source image with a trapezoidal filter.
          // in the case of no rounding this is just a box filter of width 2.
          // in the general case, the support region is 3x3.
          for (int jj = 0; jj < hSupport; jj++)
            for (int ii = 0; ii < wSupport; ii++) {
              float weight = hWeight[jj] * wWeight[ii];
              uint8_to_float(input, prevLevelMem +
                                        prevLevelPitch * (2 * j + jj) +
                                        4 * (2 * i + ii));
              result[0] += weight * input[0];
              result[1] += weight * input[1];
              result[2] += weight * input[2];
              result[3] += weight * input[3];
            }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + currLevelPitch * j + 4 * i, result);
        }
      }
    }
  }
}

}
