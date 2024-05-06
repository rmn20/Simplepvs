#ifndef _pvs_raycast_c
#define _pvs_raycast_c

#include <stdbool.h>

#include "pvs_model.c"

typedef struct {
	float start[3];
	float dir[3];
	
	float length;
	
	float min[3];
	float max[3];
} PVSRay;

typedef struct {
	size_t hitMeshId;
	float hitDistance;
	bool hit;
	bool hitBackFace;
} PVSRayHit;

bool raycastCube(PVSRayHit* hitRes, const PVSRay* ray, const float* min, const float* max);
bool rayAABBtest(const PVSRay* ray, const float* min, const float* max);
bool rayCast(PVSRayHit* hitRes, const PVSMesh* mesh, PVSRay* ray, int skipPolyId);

#endif

#if __INCLUDE_LEVEL__ == 0

#include <float.h>
#include <math.h>

bool rayAABBtest(const PVSRay* ray, const float* min, const float* max) {
	if(ray->dir[0] > 0 && max[0] < ray->start[0]) return false;
	if(ray->dir[0] < 0 && min[0] > ray->start[0]) return false;
	
	if(ray->dir[1] > 0 && max[1] < ray->start[1]) return false;
	if(ray->dir[1] < 0 && min[1] > ray->start[1]) return false;
	
	if(ray->dir[2] > 0 && max[2] < ray->start[2]) return false;
	if(ray->dir[2] < 0 && min[2] > ray->start[2]) return false;
	
	if(ray->length != FLT_MAX) {
		float t = ray->start[0] + ray->dir[0] * ray->length;
		if(ray->dir[0] > 0 && min[0] > t) return false;
		else if(ray->dir[0] < 0 && max[0] < t) return false;
		
		t = ray->start[1] + ray->dir[1] * ray->length;
		if(ray->dir[1] > 0 && min[1] > t) return false;
		else if(ray->dir[1] < 0 && max[1] < t) return false;
		
		t = ray->start[2] + ray->dir[2] * ray->length;
		if(ray->dir[2] > 0 && min[2] > t) return false;
		else if(ray->dir[2] < 0 && max[2] < t) return false;
	}
	
	return true;
}

bool raycastCube(PVSRayHit* hitRes, const PVSRay* ray, const float* min, const float* max) {
    bool insideBox = (ray->start[0] >= min[0]) && (ray->start[0] <= max[0]) &&
                     (ray->start[1] >= min[1]) && (ray->start[1] <= max[1]) &&
                     (ray->start[2] >= min[2]) && (ray->start[2] <= max[2]);
					 
	if(insideBox) {
		hitRes->hit = true;
		hitRes->hitDistance = 0;
		
		return true;
	}

	double invDir[3] = {1.0f / ray->dir[0], 1.0f / ray->dir[1], 1.0f / ray->dir[2]};

	float t[6];
    t[0] = (min[0] - ray->start[0]) * invDir[0];
    t[1] = (max[0] - ray->start[0]) * invDir[0];
	
    t[2] = (min[1] - ray->start[1]) * invDir[1];
    t[3] = (max[1] - ray->start[1]) * invDir[1];
	
    t[4] = (min[2] - ray->start[2]) * invDir[2];
    t[5] = (max[2] - ray->start[2]) * invDir[2];
	
    float tmin = fmaxf(fmaxf(fminf(t[0], t[1]), fminf(t[2], t[3])), fminf(t[4], t[5]));
    float tmax = fminf(fminf(fmaxf(t[0], t[1]), fmaxf(t[2], t[3])), fmaxf(t[4], t[5]));

	if(tmax >= tmin && tmin <= ray->length) {
		if(!hitRes->hit || tmin < hitRes->hitDistance) {
			hitRes->hit = true;
			hitRes->hitDistance = tmin;
		}
		
		return true;
	}

    return false;
}

#include <x86intrin.h>

inline bool _isPointOnPolygon2D(float px, float py, 
		float x1, float y1, float x2, float y2, float x3, float y3) {
	
	//m128
	__m128 x123 = _mm_set_ps(x1, x2, x3, 0);
	//__m128 x231 = _mm_set_ps(x2, x3, x1, 0);
    __m128 x231 = _mm_shuffle_ps(x123, x123, 0b10011100); //less instructions, but is it actually faster?
	
	x231 = _mm_sub_ps(x231, x123);
	__m128 pxxx = _mm_set_ps1(px);
	pxxx = _mm_sub_ps(pxxx, x123);
	
	__m128 y123 = _mm_set_ps(y1, y2, y3, 0);
	//__m128 y231 = _mm_set_ps(y2, y3, y1, 0);
    __m128 y231 = _mm_shuffle_ps(y123, y123, 0b10011100);
	
	y231 = _mm_sub_ps(y231, y123);
	__m128 pyyy = _mm_set_ps1(py);
	pyyy = _mm_sub_ps(pyyy, y123);
	
	x231 = _mm_mul_ps(x231, pyyy);
	y231 = _mm_mul_ps(y231, pxxx);
	
	return _mm_movemask_epi8(_mm_castps_si128(_mm_cmple_ps(x231, y231))) == 0xffff;
	
	//m256
	/*__m256 xy123 = _mm256_set_ps(x1, x2, x3, 0, y1, y2, y3, 0);
    //__m256 xy231 = _mm256_set_ps(x2, x3, x1, 0, y2, y3, y1, 0);
    __m256 xy231 = _mm256_shuffle_ps(xy123, xy123, 0b10011100);
	
	xy231 = _mm256_sub_ps(xy231, xy123);
	
	__m256 pxy = _mm256_set_m128(_mm_set_ps1(px), _mm_set_ps1(py));
	pxy = _mm256_sub_ps(pxy, xy123);
	
	xy231 = _mm256_mul_ps(xy231, _mm256_permute2f128_ps(pxy, pxy, 1));
	
	return _mm_movemask_epi8(_mm_castps_si128(_mm_cmple_ps(_mm256_extractf128_ps(xy231, 1), _mm256_extractf128_ps(xy231, 0)))) == 0xffff;*/
	
	//original
	/*return  (x2-x1)*(py-y1) <= (px-x1)*(y2-y1) &&
			(x3-x2)*(py-y2) <= (px-x2)*(y3-y2) && 
			(x1-x3)*(py-y3) <= (px-x3)*(y1-y3);*/
			
	/*
	return  !((x2-x1)*(py-y1) > (px-x1)*(y2-y1) ||
			(x3-x2)*(py-y2) > (px-x2)*(y3-y2) || 
			(x1-x3)*(py-y3) > (px-x3)*(y1-y3));
	*/
}

inline bool _isPointOnPolygon(const float* start, const float* dir, const float distance, const float* verts, const float* normal) {
	
	float nx = fabsf(normal[0]);
	float ny = fabsf(normal[1]);
	float nz = fabsf(normal[2]);

	if(nx >= ny && nx >= nz) {
		float p0 = start[2] + dir[2] * distance;
		float p1 = start[1] + dir[1] * distance;
		
		if(normal[0] >= 0) {
			return _isPointOnPolygon2D(p0, p1, verts[0 + 2], verts[0 + 1], verts[3 + 2], verts[3 + 1], verts[6 + 2], verts[6 + 1]);
		} else {
			return _isPointOnPolygon2D(p0, p1, verts[6 + 2], verts[6 + 1], verts[3 + 2], verts[3 + 1], verts[0 + 2], verts[0 + 1]);
		}
	}

	if(ny >= nx && ny >= nz) {
		float p0 = start[0] + dir[0] * distance;
		float p1 = start[2] + dir[2] * distance;
		
		if(normal[1] >= 0) {
			return _isPointOnPolygon2D(p0, p1, verts[0 + 0], verts[0 + 2], verts[3 + 0], verts[3 + 2], verts[6 + 0], verts[6 + 2]);
		} else {
			return _isPointOnPolygon2D(p0, p1, verts[6 + 0], verts[6 + 2], verts[3 + 0], verts[3 + 2], verts[0 + 0], verts[0 + 2]);
		}
	}

	if(nz >= nx && nz >= ny) {
		float p0 = start[0] + dir[0] * distance;
		float p1 = start[1] + dir[1] * distance;
		
		if(normal[2] <= 0) {
			return _isPointOnPolygon2D(p0, p1, verts[0 + 0], verts[0 + 1], verts[3 + 0], verts[3 + 1], verts[6 + 0], verts[6 + 1]);
		} else {
			return _isPointOnPolygon2D(p0, p1, verts[6 + 0], verts[6 + 1], verts[3 + 0], verts[3 + 1], verts[0 + 0], verts[0 + 1]);
		}
	}
	
	return true;
}

inline float _raycastFace(
	const float* verts, const float* norm, 
	const float* start, const float* dir, float maxDist, 
	bool* backfaceRet
	) {
	
	float dot = -(dir[0] * norm[0] + dir[1] * norm[1] + dir[2] * norm[2]);
	float dot2 = (start[0] - verts[0]) * norm[0] + (start[1] - verts[1]) * norm[1] + (start[2] - verts[2]) * norm[2];
	
	if(fabsf(dot) < 0.000001f && fabsf(dot2) > 0.000001f) return FLT_MAX; //ray is parallel to plane and isnt inside polygon
	
	float distance = dot == 0 ? 0 : dot2 / dot;
	bool backface = dot < 0;
	//if(distance == 0 && backface) return FLT_MAX;
	if(distance < 0 || distance > maxDist) return FLT_MAX;
	
	//float pos[3];
	
	//pos[0] = start[0] + dir[0] * distance;
	
	/*int axis = 0;
	if(pos[axis] < verts[axis] && pos[axis] < verts[axis + 3] && pos[axis] < verts[axis + 6]) return FLT_MAX;
	else if(pos[axis] > verts[axis] && pos[axis] > verts[axis + 3] && pos[axis] > verts[axis + 6]) return FLT_MAX;*/
	
	//pos[1] = start[1] + dir[1] * distance;
	
	/*axis = 1;
	if(pos[axis] < verts[axis] && pos[axis] < verts[axis + 3] && pos[axis] < verts[axis + 6]) return FLT_MAX;
	else if(pos[axis] > verts[axis] && pos[axis] > verts[axis + 3] && pos[axis] > verts[axis + 6]) return FLT_MAX;*/
	
	//pos[2] = start[2] + dir[2] * distance;
	
	/*axis = 2;
	if(pos[axis] < verts[axis] && pos[axis] < verts[axis + 3] && pos[axis] < verts[axis + 6]) return FLT_MAX;
	else if(pos[axis] > verts[axis] && pos[axis] > verts[axis + 3] && pos[axis] > verts[axis + 6]) return FLT_MAX;*/
	
	if(_isPointOnPolygon(start, dir, distance, verts, norm)) {
		*backfaceRet = backface;
		return distance;
	}
	
	return FLT_MAX;
}

bool rayCast(PVSRayHit* hitRes, const PVSMesh* mesh, PVSRay* ray, int skipPolyId) {
	
	bool collision = false;
	float maxDistance = ray->length;
	
	float rayMin[3];
	float rayMax[3];
	
	for(int axis=0; axis<3; axis++) {
		
		if(ray->dir[axis] > 0) {
			rayMin[axis] = ray->start[axis];
			rayMax[axis] = ray->start[axis] + ray->dir[axis] * maxDistance;
		} else if(ray->dir[axis] < 0) {
			rayMin[axis] = ray->start[axis] + ray->dir[axis] * maxDistance;
			rayMax[axis] = ray->start[axis];
		}
	}
	
	for(int i=0; i<mesh->facesCount; i++) {
		
		//if(i == skipPolyId) continue;
		
		float* verts = mesh->vertsNorms + i * 12;
		
		int axis = 0;
		if(rayMax[axis] < verts[axis] && rayMax[axis] < verts[axis + 3] && rayMax[axis] < verts[axis + 6]) continue;
		if(rayMin[axis] > verts[axis] && rayMin[axis] > verts[axis + 3] && rayMin[axis] > verts[axis + 6]) continue;
		
		axis = 1;
		if(rayMax[axis] < verts[axis] && rayMax[axis] < verts[axis + 3] && rayMax[axis] < verts[axis + 6]) continue;
		if(rayMin[axis] > verts[axis] && rayMin[axis] > verts[axis + 3] && rayMin[axis] > verts[axis + 6]) continue;
		
		axis = 2;
		if(rayMax[axis] < verts[axis] && rayMax[axis] < verts[axis + 3] && rayMax[axis] < verts[axis + 6]) continue;
		if(rayMin[axis] > verts[axis] && rayMin[axis] > verts[axis + 3] && rayMin[axis] > verts[axis + 6]) continue;
		
		bool backface;
		float dist = _raycastFace(verts, verts + 9, ray->start, ray->dir, maxDistance, &backface);

		if(dist != FLT_MAX) {
			
			if(dist < hitRes->hitDistance || !hitRes->hit || (hitRes->hit && hitRes->hitBackFace && !backface && dist == hitRes->hitDistance)) {
				hitRes->hit = true;
				hitRes->hitDistance = dist;
				hitRes->hitMeshId = mesh->id;
				hitRes->hitBackFace = backface;
				
				//_calcRayBounds(ray, maxDistance, rayMin, rayMax);
				for(int axis=0; axis<3; axis++) {
					
					if(ray->dir[axis] > 0) {
						rayMax[axis] += ray->dir[axis] * (maxDistance - dist);
					} else if(ray->dir[axis] < 0) {
						rayMin[axis] += ray->dir[axis] * (maxDistance - dist);
					}
				}
				
				collision = true;
				maxDistance = dist;
			}
		}
		
	}
	
	return collision;
}

#endif