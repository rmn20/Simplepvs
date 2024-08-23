#ifndef _pvs_c
#define _pvs_c

#include <stdlib.h>

#include "pvs_model.c"

typedef struct {
	size_t meshesCount;
	size_t visMeshesCount;
	
	char* visible;
} PVSResult;

typedef struct {
	size_t min[3];
	size_t max[3];
	
	char* visMesh;
} PVSCell;

typedef struct {
	size_t cellsX, cellsY, cellsZ;
	size_t meshes;
	
	PVSCell* cellGrid;
	
	float min[3];
	float max[3];
} PVSdb;

void pvsInitDB(PVSdb* db, PVSModel mdl, float cellSize);
void pvsUnloadDB(PVSdb* db);

size_t pvsCompute(PVSdb* db, PVSModel mdl, float raysDensity, float adaptRaysCountDistance, float* dpos, float* dhit);
PVSResult pvsGetGridVisibility(PVSdb* db, float x, float y, float z);

cvector_vector_type(PVSCell) pvsJoinCells(PVSdb* db, PVSModel mdl, size_t cellsCount, float minError,bool useFacesArea);

#endif

#if __INCLUDE_LEVEL__ == 0

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <limits.h>

#include "cvector.h"

#include "pvs_raycast.c"

void pvsInitDB(PVSdb* db, PVSModel mdl, float cellSize) {
	//Create cells
	db->cellsX = (size_t) roundf((mdl.max[0] - mdl.min[0]) / cellSize);
	db->cellsY = (size_t) roundf((mdl.max[1] - mdl.min[1]) / cellSize);
	db->cellsZ = (size_t) roundf((mdl.max[2] - mdl.min[2]) / cellSize);
	
	if(db->cellsX == 0) db->cellsX = 1;
	if(db->cellsY == 0) db->cellsY = 1;
	if(db->cellsZ == 0) db->cellsZ = 1;
	
	//db->cells = NULL;
	db->cellGrid = malloc(sizeof(PVSCell) * db->cellsX * db->cellsY * db->cellsZ);
	assert(db->cellGrid);
	
	//Fill cells
	db->meshes = cvector_size(mdl.meshes);
	size_t bytesCount = (cvector_size(mdl.meshes) / 8) + ((cvector_size(mdl.meshes) % 8) > 0 ? 1 : 0);
	
	for(size_t z=0; z<db->cellsZ; z++) {
		for(size_t y=0; y<db->cellsY; y++) {
			
			for(size_t x=0; x<db->cellsX; x++) {
				
				PVSCell* cell = db->cellGrid + (x + y * db->cellsX + z * db->cellsX * db->cellsY);
				
				cell->visMesh = calloc(bytesCount, sizeof(char));
				assert(cell->visMesh);
			}
		}
	}
	
	//Copy min and max to db
	memcpy(db->min, mdl.min, sizeof(float) * 3);
	memcpy(db->max, mdl.max, sizeof(float) * 3);
	
	db->max[1] -= 0.4f;
}

void pvsUnloadDB(PVSdb* db) {
	size_t cells = db->cellsX * db->cellsY * db->cellsZ ;
	
	for(size_t i = 0; i<cells; i++) {
		PVSCell* cell = db->cellGrid + i;
		
		free(cell->visMesh);
	}
	
	free(db->cellGrid);
}

inline void _getCellPos(PVSdb* db, const float* pos, size_t* out) {
	out[0] = (size_t) fminf(db->cellsX - 1, fmaxf(pos[0] - db->min[0], 0) * db->cellsX / (db->max[0] - db->min[0]));
	out[1] = (size_t) fminf(db->cellsY - 1, fmaxf(pos[1] - db->min[1], 0) * db->cellsY / (db->max[1] - db->min[1]));
	out[2] = (size_t) fminf(db->cellsZ - 1, fmaxf(pos[2] - db->min[2], 0) * db->cellsZ / (db->max[2] - db->min[2]));
}

inline void _getCellPosF(PVSdb* db, const float* pos, float* out) {
	out[0] = fminf(db->cellsX, fmaxf(pos[0] - db->min[0], 0) * db->cellsX / (db->max[0] - db->min[0]));
	out[1] = fminf(db->cellsY, fmaxf(pos[1] - db->min[1], 0) * db->cellsY / (db->max[1] - db->min[1]));
	out[2] = fminf(db->cellsZ, fmaxf(pos[2] - db->min[2], 0) * db->cellsZ / (db->max[2] - db->min[2]));
}

inline void _getCellBounds(PVSdb* db, const size_t* pos, float* min, float* max) {
	min[0] = pos[0] * (db->max[0] - db->min[0]) / db->cellsX + db->min[0];
	min[1] = pos[1] * (db->max[1] - db->min[1]) / db->cellsY + db->min[1];
	min[2] = pos[2] * (db->max[2] - db->min[2]) / db->cellsZ + db->min[2];
	
	max[0] = (pos[0] + 1) * (db->max[0] - db->min[0]) / db->cellsX + db->min[0];
	max[1] = (pos[1] + 1) * (db->max[1] - db->min[1]) / db->cellsY + db->min[1];
	max[2] = (pos[2] + 1) * (db->max[2] - db->min[2]) / db->cellsZ + db->min[2];
}

inline size_t mins(size_t a, size_t b) {
	if(a < b) return a;
	return b;
} 

inline size_t maxs(size_t a, size_t b) {
	if(a > b) return a;
	return b;
}

inline float randf() {
	return (float)rand() / RAND_MAX;
}

inline float pow2f(float x) {
	return x*x;
}

/*void _pvsProcessRaycastResult(PVSdb* db, const int startMeshId, const PVSRayHit* rayHit, const PVSRay* ray) {
	//Start and end cells position
	size_t startPos[3];
	_getCellPos(db, ray->start, startPos);
	
	size_t hitPos[3];
	float tmp[3] = {
		ray->start[0] + ray->dir[0] * rayHit->hitDistance,
		ray->start[1] + ray->dir[1] * rayHit->hitDistance,
		ray->start[2] + ray->dir[2] * rayHit->hitDistance
	};
	_getCellPos(db, tmp, hitPos);
	
	//Cell-space raytrace
	PVSRayHit tmpHit;
	float cellMin[3];
	float cellMax[3];
	
	size_t cellPos[3] = {startPos[0], startPos[1], startPos[2]};
	
	while(true) {
		
		PVSCell* cell = db->cellGrid + (cellPos[0] + cellPos[1] * db->cellsX + cellPos[2] * db->cellsX * db->cellsY);

		if(rayHit->hit) cell->visMesh[rayHit->hitMeshId / 8] |= 1 << (rayHit->hitMeshId % 8);
		cell->visMesh[startMeshId / 8] |= 1 << (startMeshId % 8);
		
		bool inHitPos = cellPos[0] == hitPos[0] && cellPos[1] == hitPos[1] && cellPos[2] == hitPos[2];
		if(rayHit->hit && inHitPos) break;
		
		//X plane test
		if(ray->dir[0] < 0 && cellPos[0] > 0) {
			cellPos[0] -= 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[0] += 1;
		} else if(ray->dir[0] > 0 && cellPos[0] < db->cellsX - 1) {
			cellPos[0] += 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[0] -= 1;
		}
		
		//Y plane test
		if(ray->dir[1] < 0 && cellPos[1] > 0) {
			cellPos[1] -= 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[1] += 1;
		} else if(ray->dir[1] > 0 && cellPos[1] < db->cellsY - 1) {
			cellPos[1] += 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[1] -= 1;
		}
		
		//Z plane test
		if(ray->dir[2] < 0 && cellPos[2] > 0) {
			cellPos[2] -= 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[2] += 1;
		} else if(ray->dir[2] > 0 && cellPos[2] < db->cellsZ - 1) {
			cellPos[2] += 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[2] -= 1;
		}
		
		//not hits?
		break;
	}
}*/

void _pvsProcessRaycastResult(PVSdb* db, const int startMeshId, const PVSRayHit* rayHit, const PVSRay* ray) {
	//Start and end cells position
	size_t startPos[3];
	float startPosF[3];
	_getCellPos(db, ray->start, startPos);
	_getCellPosF(db, ray->start, startPosF);
	
	size_t hitPos[3];
	float hitPosF[3];
	float tmp[3] = {
		ray->start[0] + ray->dir[0] * rayHit->hitDistance,
		ray->start[1] + ray->dir[1] * rayHit->hitDistance,
		ray->start[2] + ray->dir[2] * rayHit->hitDistance
	};
	_getCellPos(db, tmp, hitPos);
	_getCellPosF(db, tmp, hitPosF);
	
	//Prepare stuff for voxel traversal
	float maxDist = sqrtf(pow2f(hitPosF[0] - startPosF[0]) + pow2f(hitPosF[1] - startPosF[1]) + pow2f(hitPosF[2] - startPosF[2]));
	if(!rayHit->hit) maxDist = FLT_MAX;
	
	float len = sqrtf(pow2f(ray->dir[0] * db->cellsX / (db->max[0] - db->min[0])) + pow2f(ray->dir[1] * db->cellsY / (db->max[1] - db->min[1])) + pow2f(ray->dir[2] * db->cellsZ / (db->max[2] - db->min[2])));
	float dx = (ray->dir[0] * db->cellsX / (db->max[0] - db->min[0])) / len;
	float dy = (ray->dir[1] * db->cellsY / (db->max[1] - db->min[1])) / len;
	float dz = (ray->dir[2] * db->cellsZ / (db->max[2] - db->min[2])) / len;
	
	size_t ix = startPos[0];
	size_t iy = startPos[1];
	size_t iz = startPos[2];
	
	int stepx = (dx > 0) ? 1 : -1;
	int stepy = (dy > 0) ? 1 : -1;
	int stepz = (dz > 0) ? 1 : -1;
	
	float txDelta = fabs(1.0 / dx);
	float tyDelta = fabs(1.0 / dy);
	float tzDelta = fabs(1.0 / dz);
	
	float xdist = (stepx > 0) ? (ix + 1 - startPosF[0]) : (startPosF[0] - ix);
	float ydist = (stepy > 0) ? (iy + 1 - startPosF[1]) : (startPosF[1] - iy);
	float zdist = (stepz > 0) ? (iz + 1 - startPosF[2]) : (startPosF[2] - iz);
	
	float txMax = (dx != 0) ? txDelta * xdist : FLT_MAX;
	float tyMax = (dy != 0) ? tyDelta * ydist : FLT_MAX;
	float tzMax = (dz != 0) ? tzDelta * zdist : FLT_MAX;
	
	float t = 0;
	while(t <= maxDist) {
		
		PVSCell* cell = db->cellGrid + (ix + iy * db->cellsX + iz * db->cellsX * db->cellsY);

		if(rayHit->hit) cell->visMesh[rayHit->hitMeshId / 8] |= 1 << (rayHit->hitMeshId % 8);
		cell->visMesh[startMeshId / 8] |= 1 << (startMeshId % 8);
		
		bool inHitPos = ix == hitPos[0] && iy == hitPos[1] && iz == hitPos[2];
		if(rayHit->hit && inHitPos) break;
		
		// advance t to next nearest voxel boundary
		if(txMax < tyMax) {
			if(txMax < tzMax) {
				if(stepx == -1 && ix == 0) break;
				else if(stepx == 1 && ix == (db->cellsX - 1)) break;
				
				ix += stepx;
				t = txMax;
				txMax += txDelta;
			} else {
				if(stepz == -1 && iz == 0) break;
				else if(stepz == 1 && iz == (db->cellsZ - 1)) break;
				
				iz += stepz;
				t = tzMax;
				tzMax += tzDelta;
			}
		} else {
			if(tyMax < tzMax) {
				if(stepy == -1 && iy == 0) break;
				else if(stepy == 1 && iy == (db->cellsY - 1)) break;
				
				iy += stepy;
				t = tyMax;
				tyMax += tyDelta;
			} else {
				if(stepz == -1 && iz == 0) break;
				else if(stepz == 1 && iz == (db->cellsZ - 1)) break;
				
				iz += stepz;
				t = tzMax;
				tzMax += tzDelta;
			}
		}

	}
}

float _sign(float a) {
	if(a > 0) return 1;
	else if(a < 0) return -1;
	else return 0;
}


typedef struct {
	size_t mdlIndex;
	float hitDist;
} MdlToRaycast;

int mdlHitDistanceCompare(const void *a, const void *b) {
	float mdl1 = ((MdlToRaycast*) a)->hitDist;
	float mdl2 = ((MdlToRaycast*) b)->hitDist;
	
	if(mdl1 < mdl2) return -1;
	else if(mdl1 > mdl2) return 1;
	else return 0;
}

size_t pvsCompute(PVSdb* db, PVSModel mdl, float raysDensity, float adaptRaysCountDistance, float* dpos, float* dhit) {
	srand(time(NULL));
	
	size_t raysCount = 0;
	
	for(size_t meshId=0; meshId<cvector_size(mdl.meshes); meshId++) {
		PVSMesh mesh = mdl.meshes[meshId];
			
		size_t meshRaysCount = 0;
		
		//lazy multithreading..
		#pragma omp parallel
		{
			size_t startFace = mesh.facesCount * omp_get_thread_num() / omp_get_max_threads();
			size_t endFace = mesh.facesCount * (omp_get_thread_num() + 1) / omp_get_max_threads();
				
			//for lazy multithreading..
			const size_t maxResults = 1024; //processing all results shouldnt take too long!
			PVSRay hitRays[maxResults];
			PVSRayHit hitResults[maxResults];
			size_t hitResMeshId[maxResults];
			size_t hitResultsCount = 0, missCounts = 0;
			
			cvector_vector_type(MdlToRaycast) mdlsToRaycast = NULL;
			cvector_reserve(mdlsToRaycast, cvector_size(mdl.meshes));
			
			size_t localMeshRaysCount = 0;
			
			for(size_t faceId=startFace; faceId<endFace; faceId++) {
				
				float* verts = mesh.vertsNorms + faceId * 4 * 3;
				float* faceNorm = verts + 9;
				
				//Calculate rays count for current polygon size
				float distAB = sqrtf(pow2f(verts[0] - verts[3]) + pow2f(verts[1] - verts[3 + 1]) + pow2f(verts[2] - verts[3 + 2]));
				float distBC = sqrtf(pow2f(verts[3] - verts[6]) + pow2f(verts[3 + 1] - verts[6 + 1]) + pow2f(verts[3 + 2] - verts[6 + 2]));
				float distAC = sqrtf(pow2f(verts[0] - verts[6]) + pow2f(verts[1] - verts[6 + 1]) + pow2f(verts[2] - verts[6 + 2]));
				
				float S = (distAB + distBC + distAC) / 2;
				float faceArea = sqrtf(S * (S-distAB) * (S-distBC) * (S-distAC) / 2);
				//float faceArea = distAB + distBC + distAC;
				
				size_t faceRays = (size_t) round(raysDensity * faceArea);
				if(faceRays == 0) faceRays = 1;
				
				localMeshRaysCount += faceRays;
				
				//Used for adaptive rays count
				float maxHitDist = 1;
				size_t adaptFaceRays = faceRays;
				
				//Cast ray on a per mesh basis
				for(size_t i=0; i<adaptFaceRays; i++) {
					
					PVSRay ray;
				
					//Starting position inside triangle
					float a_r = randf();
					float b_r = randf();
					
					for(int axis=0; axis<3; axis++) {
						ray.start[axis] = (1 - sqrt(a_r)) * verts[axis] + (sqrt(a_r) * (1 - b_r)) * verts[3 + axis] + (b_r * sqrt(a_r)) * verts[6 + axis];
						//ray.start[axis] += faceNorm[axis] * 0.01;
						ray.start[axis] += _sign(faceNorm[axis]) * fabs(ray.start[axis] * 0.0000002);
					}
					
					bool skip = false;
					for(int axis=0; axis<3; axis++) {
						skip |= ray.start[axis] > db->max[axis];
						skip |= ray.start[axis] < db->min[axis];
						
						if(skip) break;
					}
					
					if(skip) continue;
					
					//Starting direction
					//Generating uniformly distributed random directions 
					while(true) {
						ray.dir[0] = (randf() + randf()) * 2 - 1;
						ray.dir[1] = (randf() + randf()) * 2 - 1;
						ray.dir[2] = (randf() + randf()) * 2 - 1;
						
						float len = sqrtf(pow2f(ray.dir[0]) + pow2f(ray.dir[1]) + pow2f(ray.dir[2]));
						if(len == 0) continue;
						
						ray.dir[0] /= len;
						ray.dir[1] /= len;
						ray.dir[2] /= len;
						break;
					}
					
					//Rotate direction toward face normal
					float dot = ray.dir[0] * faceNorm[0] + ray.dir[1] * faceNorm[1] + ray.dir[2] * faceNorm[2];
					
					if(dot < 0) {
						ray.dir[0] *= -1;
						ray.dir[1] *= -1;
						ray.dir[2] *= -1;
					}
					
					//Limit ray length to scene AABB for faster raycasting
					ray.length = FLT_MAX;
					
					if(ray.dir[0] > 0) ray.length = (db->max[0] - ray.start[0]) / ray.dir[0];
					else if(ray.dir[0] < 0) ray.length = (db->min[0] - ray.start[0]) / ray.dir[0];
					
					if(ray.dir[1] > 0) ray.length = fminf(ray.length, (db->max[1] - ray.start[1]) / ray.dir[1]);
					else if(ray.dir[1] < 0) ray.length = fminf(ray.length, (db->min[1] - ray.start[1]) / ray.dir[1]);
					
					if(ray.dir[2] > 0) ray.length = fminf(ray.length, (db->max[2] - ray.start[2]) / ray.dir[2]);
					else if(ray.dir[2] < 0) ray.length = fminf(ray.length, (db->min[2] - ray.start[2]) / ray.dir[2]);
					
					//Form list of meshes to raycast
					cvector_clear(mdlsToRaycast);
					
					for(size_t t=0; t<cvector_size(mdl.meshes); t++) {
						PVSMesh* meshCast = mdl.meshes + t;
						PVSRayHit tmpHit = {0};
						
						if(!rayAABBtest(&ray, meshCast->min, meshCast->max)) continue;
						if(!raycastCube(&tmpHit, &ray, meshCast->min, meshCast->max)) continue;
						
						MdlToRaycast mdlToRaycast = {t, tmpHit.hitDistance};
						cvector_push_back(mdlsToRaycast, mdlToRaycast);
					}
					
					qsort(mdlsToRaycast, cvector_size(mdlsToRaycast), sizeof(MdlToRaycast), &mdlHitDistanceCompare);
					
					//Raycasting
					PVSRayHit rayHit = {0};
					
					for(int t=0; t<cvector_size(mdlsToRaycast); t++) {
						size_t mdlIndex = mdlsToRaycast[t].mdlIndex;
						PVSMesh* meshCast = mdl.meshes + mdlIndex;
						
						if(rayHit.hit && ray.length < mdlsToRaycast[t].hitDist) break;
					
						bool hit = rayCast(&rayHit, meshCast, &ray, mdlIndex==meshId?faceId:-1);
						if(hit) ray.length = rayHit.hitDistance;
					}
				
					if(!rayHit.hitBackFace) {
					//if(!rayHit.hitBackFace && (rayHit.hit || ray.dir[1] > 0)) {
					//if(!rayHit.hitBackFace && rayHit.hit) {//(rayHit.hit || ray.start[1] + ray.dir[1] * ray.length >= db->max[1])) {
						hitResults[hitResultsCount] = rayHit;
						hitRays[hitResultsCount] = ray;
						hitResMeshId[hitResultsCount] = mesh.id;
						hitResultsCount++;
					} else missCounts++;
					
					//Update required rays count
					if(rayHit.hit && !rayHit.hitBackFace && rayHit.hitDistance / adaptRaysCountDistance > maxHitDist) {
						maxHitDist = rayHit.hitDistance / adaptRaysCountDistance;
						
						adaptFaceRays = (size_t) round(faceRays * pow2f(maxHitDist));
					}
					
					if(!rayHit.hit && ray.length / adaptRaysCountDistance > maxHitDist) {
						maxHitDist = ray.length / adaptRaysCountDistance;
						
						adaptFaceRays = (size_t) round(faceRays * pow2f(maxHitDist));
					}

					//lazy multithreading..
					if(hitResultsCount == maxResults) {
						#pragma omp critical
						{
							for(size_t resId = 0; resId < hitResultsCount; resId++) {
								_pvsProcessRaycastResult(db, hitResMeshId[resId], hitResults + resId, hitRays + resId);
								
								//if(hitResMeshId[resId] == 11 && (hitRays[resId].start[2] + hitRays[resId].dir[2] * hitRays[resId].length) < -40 && hitRays[resId].start[1] < 5) {
								/*if((hitResMeshId[resId] == 0 && hitResMeshHitId[resId] == 3) || (hitResMeshId[resId] == 3 && hitResMeshHitId[resId] == 0)) {
									dpos[0] = hitRays[resId].start[0];
									dpos[1] = hitRays[resId].start[1];
									dpos[2] = hitRays[resId].start[2];
									
									dhit[0] = dpos[0] + hitRays[resId].dir[0] * hitRays[resId].length;
									dhit[1] = dpos[1] + hitRays[resId].dir[1] * hitRays[resId].length;
									dhit[2] = dpos[2] + hitRays[resId].dir[2] * hitRays[resId].length;
								}*/
							}
							
							raysCount += hitResultsCount + missCounts;
						}
						
						hitResultsCount = 0;
						missCounts = 0;
					}
				}
			}
			
			#pragma omp barrier
			#pragma omp critical
			{
				for(size_t resId = 0; resId < hitResultsCount; resId++) {
					if(!hitResults[resId].hitBackFace) _pvsProcessRaycastResult(db, hitResMeshId[resId], hitResults + resId, hitRays + resId);
								
					/*if(hitResMeshId[resId] == 11 && (hitRays[resId].start[2] + hitRays[resId].dir[2] * hitRays[resId].length) < -40 && hitRays[resId].start[1] < 5) {
						dpos[0] = hitRays[resId].start[0];
						dpos[1] = hitRays[resId].start[1];
						dpos[2] = hitRays[resId].start[2];
									
						dhit[0] = dpos[0] + hitRays[resId].dir[0] * hitRays[resId].length;
						dhit[1] = dpos[1] + hitRays[resId].dir[1] * hitRays[resId].length;
						dhit[2] = dpos[2] + hitRays[resId].dir[2] * hitRays[resId].length;
					}*/
				}
				
				raysCount += hitResultsCount;
				meshRaysCount += localMeshRaysCount;
			}
			
			cvector_free(mdlsToRaycast);
			
			printf("expected rays: %llu, real rays count: %llu\n", (size_t) round(raysDensity * mesh.facesArea), meshRaysCount);
		}
		
		int printEvery = cvector_size(mdl.meshes) / 10;
		if(printEvery == 0) printEvery = 1;
		
		if((meshId % printEvery) == 0) {
			printf("%d %%\n", (int) (meshId * 100 / cvector_size(mdl.meshes)));
		}
	}
	
	return raysCount;
}

PVSResult pvsGetGridVisibility(PVSdb* db, float x, float y, float z) {
	PVSResult res = {};
	
	float tmp[3] = {x, y, z};
	size_t cellPos[3] = {};
	_getCellPos(db, tmp, cellPos);
	
	res.meshesCount = db->meshes;
	res.visMeshesCount = 0;
	res.visible = db->cellGrid[cellPos[0] + cellPos[1] * db->cellsX + cellPos[2] * db->cellsX * db->cellsY].visMesh;
	
	for(size_t i=0; i<res.meshesCount; i++) {
		res.visMeshesCount += (res.visible[i / 8] & (1 << (i % 8))) != 0;
	}
	
	return res;
}

size_t _cellsVisMeshesCount(PVSdb* db, char* tmpVisMeshes, size_t visMeshesSize, size_t* min, size_t* max) {
	
	memset(tmpVisMeshes, 0, visMeshesSize * sizeof(char));
	
	for(size_t x=min[0]; x<max[0]; x++) {
		for(size_t y=min[1]; y<max[1]; y++) {
			for(size_t z=min[2]; z<max[2]; z++) {
				
				PVSCell* cell = db->cellGrid + (x + y * db->cellsX + z * db->cellsX * db->cellsY);
				
				for(size_t i=0; i<visMeshesSize; i++) {
					tmpVisMeshes[i] |= cell->visMesh[i];
				}
				
			}
		}
	}
	
	size_t visCount = 0;
	
	for(size_t i=0; i<visMeshesSize; i++) {
		int cellVisible = tmpVisMeshes[i];
		
		for(int t=0; t<8; t++) visCount += (cellVisible >> t) & 1;
	}
	
	return visCount;
}

double _cellsError(PVSdb* db, PVSModel* mdl, char* tmpVisMeshes, size_t visMeshesSize, size_t* min, size_t* max, bool useFacesArea) {
	//Calculate error in cell min-max
	
	memset(tmpVisMeshes, 0, visMeshesSize * sizeof(char));
	
	for(size_t x=min[0]; x<max[0]; x++) {
		for(size_t y=min[1]; y<max[1]; y++) {
			for(size_t z=min[2]; z<max[2]; z++) {
				
				PVSCell* cell = db->cellGrid + (x + y * db->cellsX + z * db->cellsX * db->cellsY);
				
				for(size_t i=0; i<visMeshesSize; i++) {
					tmpVisMeshes[i] |= cell->visMesh[i];
				}
				
			}
		}
	}
	
	double error = 0;
	
	for(size_t x=min[0]; x<max[0]; x++) {
		for(size_t y=min[1]; y<max[1]; y++) {
			for(size_t z=min[2]; z<max[2]; z++) {
				
				PVSCell* cell = db->cellGrid + (x + y * db->cellsX + z * db->cellsX * db->cellsY);
				
				for(size_t i=0; i<db->meshes; i++) {
					int meshError = (tmpVisMeshes[i / 8] ^ cell->visMesh[i / 8]) >> (i&7);
					meshError &= 1;
					
					error += meshError * mdl->meshes[i].facesArea * mdl->meshes[i].facesCount;//(useFacesArea ? mdl->meshes[i].facesArea : mdl->meshes[i].facesCount);
				}
				/*for(size_t i=0; i<visMeshesSize; i++) {
					int cellError = tmpVisMeshes[i] ^ cell->visMesh[i];
					
					for(int t=0; t<8; t++) error += (cellError >> t) & 1;
				}*/
				
			}
		}
	}
	
	return error;
}

inline void _cellFindSplit(PVSdb* db, PVSModel* mdl, char* tmpVisMeshes, size_t visMeshesSize, PVSCell* cell, int* splitAxis, double* minSplitError, double* minMaxSplitError, size_t* splitPos, bool useFacesArea) {
	//Find a axis aligned split plane
	*splitAxis = -1;
	*minSplitError = 0;
	*minMaxSplitError = 0;
	*splitPos = 0;
	
	for(int axis=0; axis<3; axis++) {
		
		size_t minPos = cell->min[axis];
		size_t maxPos = cell->max[axis];
		
		for(size_t pos=minPos+1; pos<maxPos; pos++) {
			
			size_t splitMin[3] = {cell->min[0], cell->min[1], cell->min[2]};
			size_t splitMax[3] = {cell->max[0], cell->max[1], cell->max[2]};
			
			splitMax[axis] = pos;
			double err1 = _cellsError(db, mdl, tmpVisMeshes, visMeshesSize, splitMin, splitMax, useFacesArea);
			
			splitMin[axis] = pos;
			splitMax[axis] = cell->max[axis];
			double err2 = _cellsError(db, mdl, tmpVisMeshes, visMeshesSize, splitMin, splitMax, useFacesArea);
			
			if(*splitAxis == -1 || (err1 + err2) < *minSplitError) {
				*splitAxis = axis;
				*minSplitError = err1 + err2;
				*minMaxSplitError = err1 > err2 ? err1 : err2;
				*splitPos = pos;
				
				if(*minSplitError == 0) return;
			}
		}
	}
}

cvector_vector_type(PVSCell) pvsJoinCells(PVSdb* db, PVSModel mdl, size_t cellsCount, float minError, bool splitClever) {
	
	cvector_vector_type(PVSCell) cells = NULL;
	
	size_t visMeshesSize = (db->meshes / 8) + ((db->meshes % 8) > 0 ? 1 : 0);
	char* tmpVisMeshes = malloc(visMeshesSize);
	assert(tmpVisMeshes);
	
	PVSCell firstCell = {0};
	firstCell.max[0] = db->cellsX;
	firstCell.max[1] = db->cellsY;
	firstCell.max[2] = db->cellsZ;
	
	cvector_push_back(cells, firstCell);
	
	bool useFacesArea = true;
	
	while(cvector_size(cells) < cellsCount) {
		
		//Find a cell with maximum total error
		bool splitCellExist = false;
		size_t splitCellId;
		float splitCellError = 0;
		
		for(size_t i=0; i<cvector_size(cells); i++) {
			PVSCell* cell = cells + i;
			
			//cell can't be split
			//if(cell->min[0] == cell->max[0] + 1 && cell->min[1] == cell->max[1] + 1 && cell->min[2] + 1 == cell->max[2]) continue;
			//error will be 0 anyway
			
			double err = _cellsError(db, &mdl, tmpVisMeshes, visMeshesSize, cell->min, cell->max, useFacesArea);
			if(err < splitCellError) continue;
			if(err < minError) continue;
			
			//Find an axis aligned split plane
			if(splitClever) {
				int splitAxis = -1;
				double minSplitError;
				double minMaxSplitError;
				size_t splitPos;
				_cellFindSplit(db, &mdl, tmpVisMeshes, visMeshesSize, cell, &splitAxis, &minSplitError, &minMaxSplitError, &splitPos, useFacesArea);
				
				if(minSplitError > err) printf("wtf\n");
				
				err = err - minSplitError;
			}
			
			if(err > splitCellError) {
				splitCellExist = true;
				splitCellId = i;
				splitCellError = err;
			}
		}
		
		//There is no point in splitting more cells
		if(!splitCellExist) break;
		PVSCell* splitCell = cells + splitCellId;
		
		//Find a axis aligned split plane
		int splitAxis = -1;
		double minSplitError;
		double minMaxSplitError;
		size_t splitPos;
		_cellFindSplit(db, &mdl, tmpVisMeshes, visMeshesSize, splitCell, &splitAxis, &minSplitError, &minMaxSplitError, &splitPos, useFacesArea);
		
		//There is nothing to split (this shouldnt be possible?)
		//if(splitAxis == -1) break;
		
		PVSCell cell2 = {0};
		for(int axis=0; axis<3; axis++) {
			cell2.min[axis] = splitCell->min[axis];
			cell2.max[axis] = splitCell->max[axis];
		}
		//memcpy(cell2.min, splitCell->min, sizeof(float) * 3);
		//memcpy(cell2.max, splitCell->max, sizeof(float) * 3);
		
		splitCell->max[splitAxis] = splitPos;
		cell2.min[splitAxis] = splitPos;
		
		cvector_push_back(cells, cell2);
		
		if(_cellsVisMeshesCount(db, tmpVisMeshes, visMeshesSize, cell2.min, cell2.max) == 0) cvector_pop_back(cells);
		if(_cellsVisMeshesCount(db, tmpVisMeshes, visMeshesSize, splitCell->min, splitCell->max) == 0) cvector_erase(cells, splitCellId);
	}
	
	//Compute visibility for every cell
	for(size_t i=0; i<cvector_size(cells); i++) {
		PVSCell* cell = cells + i;
		
		cell->visMesh = calloc(visMeshesSize, sizeof(char));
		assert(cell->visMesh);
		
		size_t* min = cell->min;
		size_t* max = cell->max;
		
		for(size_t x=min[0]; x<max[0]; x++) {
			for(size_t y=min[1]; y<max[1]; y++) {
				for(size_t z=min[2]; z<max[2]; z++) {
					
					PVSCell* cell2 = db->cellGrid + (x + y * db->cellsX + z * db->cellsX * db->cellsY);
					
					for(size_t t=0; t<visMeshesSize; t++) {
						cell->visMesh[t] |= cell2->visMesh[t];
					}
					
				}
			}
		}
	}
	
	return cells;
}

#endif