#ifndef _pvs_model_c
#define _pvs_model_c

#include <stdbool.h>
#include <stdlib.h>

#include "cvector.h"

//todo mesh, model unload!

typedef struct {
	size_t id;
	float min[3];
	float max[3];
	double facesArea;
	
	int facesCount;
	float* vertsNorms;
} PVSMesh;

typedef struct {
	cvector_vector_type(PVSMesh) meshes;
	
	float min[3];
	float max[3];
} PVSModel;

void calcNorm(float* normal, float* face);
bool pvsLoadObj(PVSModel* mdl, char* inputPath);

#endif

#if __INCLUDE_LEVEL__ == 0

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

//Calculate face normal
void calcNorm(float* normal, float* face) {
	float ax = face[0], ay = face[1], az = face[2];
	float bx = face[3], by = face[4], bz = face[5];
	float cx = face[6], cy = face[7], cz = face[8];

	float x = (by - ay) * (cz - az) - (bz - az) * (cy - ay);
	float y = (bz - az) * (cx - ax) - (bx - ax) * (cz - az);
	float z = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);

	double len = sqrt(x * x + y * y + z * z);
	
	if(len == 0) {
		//found a degenerate..
		normal[0] = 0;
		normal[1] = 1;
		normal[2] = 0;
	} else {
		normal[0] = (float) (x / len);
		normal[1] = (float) (y / len);
		normal[2] = (float) (z / len);
	}
}

//Load vertex data
cvector_vector_type(float) _objLoadVertex(char* line, cvector_vector_type(float) vertices, size_t attribs) {
	
	char* token = strtok(line, " ");
	size_t tokens = 0;
	
    while(token) {
        if(tokens > 0 && tokens-1 < attribs) cvector_push_back(vertices, (float) atof(token));
        token = strtok(NULL, " ");
		tokens++;
    }
	
	return vertices;
}

//Load face data
cvector_vector_type(float) _objLoadFace(char* line, cvector_vector_type(float) vertices, cvector_vector_type(float) meshVerts) {
	
	char* token = strtok(line, " ");
	size_t tokens = 0;
	size_t vert = 0;
	
	//parse vertex by vertex
    while(token) {
		
        if(tokens > 0) {
			int v = 0;//, uv = 0, n = 0;
			
			char* tmp = token;
			size_t attrib = 0;
			
			//Parse every attribute
			while(tmp) {
				if(attrib == 0) v = atoi(tmp);
				/*else if(attrib == 1) uv = atoi(tmp);
				else if(attrib == 2) n = atoi(tmp);*/
				
				attrib++;
				tmp = strchr(tmp, '/');
				if(tmp) tmp++;
			}
			
			//polygon triangulation!
			if(vert >= 3) {
				for(int i=0; i<3; i++) cvector_push_back(meshVerts, meshVerts[cvector_size(meshVerts) - 3*vert]);
				for(int i=0; i<3; i++) cvector_push_back(meshVerts, meshVerts[cvector_size(meshVerts) - 6]);

				vert += 2;
			}
			
			//adding vertices to a mesh
			if(v < 0) {
				for(int i=0; i<3; i++) cvector_push_back(meshVerts, vertices[cvector_size(vertices) + v * 3 + i]);
			} else if(v > 0) {
				for(int i=0; i<3; i++) cvector_push_back(meshVerts, vertices[(v-1) * 3 + i]);
			} else {
				//???
			}
			
			vert++;
		}
		
        token = strtok(NULL, " ");
		tokens++;
    }
	
	return meshVerts;
}

float inline pow2f(float a) {
	return a*a;
}

double inline pow2d(double a) {
	return a*a;
}

//Add mesh to a model
cvector_vector_type(float) _objAddMesh(PVSModel* mdl, PVSMesh* meshData, cvector_vector_type(float) meshVerts) {
	
	//Calculate AABB
	float mins[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
	float maxs[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
	
	for(int i=0; i<cvector_size(meshVerts); i+=3) {
		
		for(int axis=0; axis<3; axis++) {
			mins[axis] = fminf(meshVerts[i + axis], mins[axis]);
			maxs[axis] = fmaxf(meshVerts[i + axis], maxs[axis]);
		}
	}
	
	memcpy(&(meshData->min), mins, sizeof(float) * 3);
	memcpy(&(meshData->max), maxs, sizeof(float) * 3);
	
	//Copy vertex data
	meshData->facesCount = cvector_size(meshVerts) / 9;
	meshData->vertsNorms = malloc(sizeof(float) * (cvector_size(meshVerts) * 4 / 3) );
	assert(meshData->vertsNorms);
	
	meshData->facesArea = 0;
	
	for(size_t i=0; i<cvector_size(meshVerts) / 9; i++) {
		
		float* verts = meshData->vertsNorms + i * 12;
		
		memcpy(verts, meshVerts + i * 9, 9 * sizeof(float));
		calcNorm(verts + 9, verts);
		
		{
			float distAB = sqrtf(pow2f(verts[0] - verts[3]) + pow2f(verts[1] - verts[3 + 1]) + pow2f(verts[2] - verts[3 + 2]));
			float distBC = sqrtf(pow2f(verts[3] - verts[6]) + pow2f(verts[3 + 1] - verts[6 + 1]) + pow2f(verts[3 + 2] - verts[6 + 2]));
			float distAC = sqrtf(pow2f(verts[0] - verts[6]) + pow2f(verts[1] - verts[6 + 1]) + pow2f(verts[2] - verts[6 + 2]));
			
			float S = (distAB + distBC + distAC) / 2;
			float faceArea = sqrtf(S * (S-distAB) * (S-distBC) * (S-distAC) / 2);
			
			meshData->facesArea += faceArea;
		}
	}
	
	cvector_clear(meshVerts);
	
	//Add mesh to model and update mesh data
	cvector_push_back(mdl->meshes, *meshData);
	
	meshData->id++;
	meshData->vertsNorms = NULL;
	
	return meshVerts;
}

//load obj wavefront model
bool pvsLoadObj(PVSModel* mdl, char* inputPath) {
	
	FILE* f = fopen(inputPath, "rt");
	if(!f) {
		printf("Can't open file\n");
		
		return false;
	}
	
	cvector_vector_type(float) vertices = NULL;
	cvector_vector_type(float) meshVerts = NULL;
	
	PVSMesh meshData;
	meshData.id = 0;
	
	//Parse line by line
	char line[1024] = {'\0'};
	
	while(fgets(line, 1024, f)) {
		size_t size = strlen(line) - 1; //skip new line character
		if(!size) continue;
		
		if(strstr(line, "o ") == line || strstr(line, "g ") == line) {
			
			if(cvector_size(meshVerts) > 0) {
				meshVerts = _objAddMesh(mdl, &meshData, meshVerts);
			}
			
		} else if(strstr(line, "v ") == line) {
			vertices = _objLoadVertex(line, vertices, 3);
		} else if(strstr(line, "f ") == line) {
			meshVerts = _objLoadFace(line, vertices, meshVerts);
		}
		
	}
	
	//Add last vertices to a mesh
	if(cvector_size(meshVerts) > 0) {
		meshVerts = _objAddMesh(mdl, &meshData, meshVerts);
	}
	
	fclose(f);
	
	//AABB
	float mins[3] = {FLT_MAX};
	float maxs[3] = {-FLT_MAX};
	
	for(size_t i=0; i<cvector_size(mdl->meshes); i++) {
		
		for(int axis=0; axis<3; axis++) {
			mins[axis] = fminf(mdl->meshes[i].min[axis], mins[axis]);
			maxs[axis] = fmaxf(mdl->meshes[i].max[axis], maxs[axis]);
		}
	}
	
	memcpy(mdl->min, mins, sizeof(float) * 3);
	memcpy(mdl->max, maxs, sizeof(float) * 3);
	
	printf("Meshes: %d\n", (int) cvector_size(mdl->meshes));
	
	cvector_free(vertices);
	cvector_free(meshVerts);

	return true;
}

#endif