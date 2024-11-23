#ifndef _pvs_c
#define _pvs_c

#include <stdlib.h>
#include <stdint.h>

#include "cvector.h"

#include "raylib.h"

typedef struct {
	size_t min[3];
	size_t max[3];
	double error;
	
	uint32_t* visMesh;
} PVSCell;

typedef struct {
	float min[3];
	float max[3];
	
	float facesArea;
	size_t faces;
} PVSMeshData;

typedef struct {
	size_t cellsX, cellsY, cellsZ;
	PVSCell* cellGrid;
	
	size_t meshes;
	PVSMeshData* meshesData;
	
	float min[3];
	float max[3];
} PVSdb;

typedef struct {
	size_t meshesCount;
	size_t visMeshesCount;
	
	uint32_t* visible;
} PVSResult;

void pvsInitDB(PVSdb* db, float cellSize, Model* mdl);
void pvsUnloadDB(PVSdb* db);

size_t pvsCompute(PVSdb* db, Model* mdl, Vector3 camPos);
PVSResult pvsGetGridVisibility(PVSdb* db, float x, float y, float z);

cvector_vector_type(PVSCell) pvsJoinCells(PVSdb* db, size_t cellsCount, float minError, bool useFacesArea);

#endif

#if __INCLUDE_LEVEL__ == 0

#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <limits.h>

#include "cvector.h"

#include "rlgl.h"
#include "raymath.h"

inline size_t _getCellPtr(PVSdb* db, size_t x, size_t y, size_t z) {
	return x + y * db->cellsX + z * db->cellsX * db->cellsY;
}

void pvsInitDB(PVSdb* db, float cellSize, Model* mdl) {
	//Compute meshes data
	db->meshes = mdl->meshCount;
	db->meshesData = malloc(sizeof(PVSMeshData) * db->meshes);
	assert(db->meshesData);
	
	for(int axis = 0; axis < 3; axis++) {
		db->min[axis] = FLT_MAX;
		db->max[axis] = -FLT_MAX;
	}
	
	for(int i = 0; i < mdl->meshCount; i++) {
		PVSMeshData* meshData = db->meshesData + i;
		const Mesh mesh = mdl->meshes[i];
		
		BoundingBox aabb = GetMeshBoundingBox(mesh);
		
		memcpy(meshData->min, &aabb.min, sizeof(Vector3));
		memcpy(meshData->max, &aabb.max, sizeof(Vector3));
		
		for(int axis = 0; axis < 3; axis++) {
			db->min[axis] = fminf(db->min[axis], meshData->min[axis]);
			db->max[axis] = fmaxf(db->max[axis], meshData->max[axis]);
		}
		
		meshData->faces = mesh.triangleCount;
		meshData->facesArea = 0;
		
		Vector3* verts = (Vector3*) mesh.vertices;
		
		for(int i = 0; i < mesh.triangleCount; i++) {
			Vector3 a, b, c;
			
			if(mesh.indices != NULL) {
				a = verts[mesh.indices[i * 3]];
				b = verts[mesh.indices[i * 3 + 1]];
				c = verts[mesh.indices[i * 3 + 2]];
			} else {
				a = verts[i * 3];
				b = verts[i * 3 + 1];
				c = verts[i * 3 + 2];
			}
			
			float distAB = Vector3Distance(a, b);
			float distBC = Vector3Distance(b, c);
			float distAC = Vector3Distance(a, c);
			
			float S = (distAB + distBC + distAC) / 2;
			meshData->facesArea += sqrtf(S * (S-distAB) * (S-distBC) * (S-distAC) / 2);
		}
	}
	
	//db->max[1] -= 0.4f;
	
	//Create cells
	db->cellsX = (size_t) roundf((db->max[0] - db->min[0]) / cellSize);
	db->cellsY = (size_t) roundf((db->max[1] - db->min[1]) / cellSize);
	db->cellsZ = (size_t) roundf((db->max[2] - db->min[2]) / cellSize);
	
	if(db->cellsX == 0) db->cellsX = 1;
	if(db->cellsY == 0) db->cellsY = 1;
	if(db->cellsZ == 0) db->cellsZ = 1;
	
	//db->cells = NULL;
	db->cellGrid = malloc(sizeof(PVSCell) * db->cellsX * db->cellsY * db->cellsZ);
	assert(db->cellGrid);
	
	//Fill cells
	size_t dataPerCell = (db->meshes / 32) + ((db->meshes % 32) > 0 ? 1 : 0);
	
	for(size_t z=0; z<db->cellsZ; z++) {
		for(size_t y=0; y<db->cellsY; y++) {
			for(size_t x=0; x<db->cellsX; x++) {
				PVSCell* cell = db->cellGrid + _getCellPtr(db, x, y, z);
				
				cell->visMesh = calloc(dataPerCell, sizeof(uint32_t));
				assert(cell->visMesh);
			}
		}
	}
}

void pvsUnloadDB(PVSdb* db) {
	size_t cells = db->cellsX * db->cellsY * db->cellsZ ;
	
	for(size_t i = 0; i<cells; i++) {
		PVSCell* cell = db->cellGrid + i;
		
		free(cell->visMesh);
	}
	
	free(db->cellGrid);
	free(db->meshesData);
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

inline float randf() {
	return (float)rand() / RAND_MAX;
}

inline float pow2f(float x) {
	return x*x;
}

typedef struct GpuData {
    float matVP[16 * 6];
    float matVPInv[16 * 6];
	
	uint32_t gridSize[4];
	float aabbMin[4];
	float aabbMax[4];
	/*uint32_t outAABB[4 * 6];
	uint32_t camPos[4];*/
	uint32_t dataPerCell;
} GpuData;

size_t pvsCompute(PVSdb* db, Model* mdl, Vector3 camPos) {
	srand(time(NULL));
	
	//Create render target
	int size = 128; //should be multiple of 8
	RenderTexture2D target = {0};

	target.id = rlLoadFramebuffer(size * 3, size * 2); // Load an empty framebuffer
	assert(target.id > 0);
	rlEnableFramebuffer(target.id);
	
	//Create render textures
	target.texture.id = rlLoadTexture(NULL, size * 3, size * 2, RL_PIXELFORMAT_UNCOMPRESSED_R16, 1);
	target.texture.width = size * 3;
	target.texture.height = size * 2;
	target.texture.format = RL_PIXELFORMAT_UNCOMPRESSED_R16;
	target.texture.mipmaps = 1;
	
	rlTextureParameters(target.texture.id, RL_TEXTURE_MIN_FILTER, RL_TEXTURE_FILTER_NEAREST);
	rlTextureParameters(target.texture.id, RL_TEXTURE_MAG_FILTER, RL_TEXTURE_FILTER_NEAREST);
	
	target.depth.id = rlLoadTextureDepth(size * 3, size * 2, false);
	target.depth.width = size * 3;
	target.depth.height = size * 2;
	target.depth.format = 19; //DEPTH_COMPONENT_24BIT?
	target.depth.mipmaps = 1;
	
	rlTextureParameters(target.depth.id, RL_TEXTURE_MIN_FILTER, RL_TEXTURE_FILTER_NEAREST);
	rlTextureParameters(target.depth.id, RL_TEXTURE_MAG_FILTER, RL_TEXTURE_FILTER_NEAREST);
	
	//Attach depth texture to FBO
	rlFramebufferAttach(target.id, target.texture.id, RL_ATTACHMENT_COLOR_CHANNEL0, RL_ATTACHMENT_TEXTURE2D, 0);
	rlFramebufferAttach(target.id, target.depth.id, RL_ATTACHMENT_DEPTH, RL_ATTACHMENT_TEXTURE2D, 0);

	//Check if fbo is complete with attachments (valid)
	assert(rlFramebufferComplete(target.id));
	rlDisableFramebuffer();
	
	//Load material and shader
	Shader renderShader = LoadShader("shaders/cubemap-render.vs", "shaders/cubemap-render.fs");
	int meshIdUniform = GetShaderLocation(renderShader, "meshId");
	
	Material mat = LoadMaterialDefault();
	mat.shader = renderShader;
	
	//Load cell visibility fill shader
	char *cellVisShaderText = LoadFileText("shaders/cubemap2cell.glsl");
    unsigned int cellVisShaderId = rlCompileShader(cellVisShaderText, RL_COMPUTE_SHADER);
    unsigned int cellVisProgram = rlLoadComputeShaderProgram(cellVisShaderId);
    UnloadFileText(cellVisShaderText);
	
	//Create data struct for use in compute shader
	GpuData gpuData = {0};
	
	gpuData.gridSize[0] = db->cellsX;
	gpuData.gridSize[1] = db->cellsY;
	gpuData.gridSize[2] = db->cellsZ;
	
	gpuData.aabbMin[0] = db->min[0];
	gpuData.aabbMin[1] = db->min[1];
	gpuData.aabbMin[2] = db->min[2];
	
	gpuData.aabbMax[0] = db->max[0];
	gpuData.aabbMax[1] = db->max[1];
	gpuData.aabbMax[2] = db->max[2];
	
	/*gpuData.camPos[0] = (camPos.x - db->min[0]) * db->cellsX / (db->max[0] - db->min[0]);
	gpuData.camPos[1] = (camPos.y - db->min[1]) * db->cellsY / (db->max[1] - db->min[1]);
	gpuData.camPos[2] = (camPos.z - db->min[2]) * db->cellsZ / (db->max[2] - db->min[2]);*/
	
	gpuData.dataPerCell = (db->meshes / 32) + ((db->meshes % 32) > 0 ? 1 : 0);
	
	size_t visBufferSize = db->cellsX * db->cellsY * db->cellsZ * gpuData.dataPerCell;
	uint32_t* tmpVisBuffer = malloc(visBufferSize * sizeof(uint32_t));
	assert(tmpVisBuffer);
	
	for(size_t z=0; z<db->cellsZ; z++) {
		for(size_t y=0; y<db->cellsY; y++) {
			for(size_t x=0; x<db->cellsX; x++) {
				
				PVSCell* cell = db->cellGrid + _getCellPtr(db, x, y, z);
				
				memcpy(
					tmpVisBuffer + _getCellPtr(db, x, y, z) * gpuData.dataPerCell, 
					cell->visMesh, 
					gpuData.dataPerCell * sizeof(uint32_t)
				);
			}
		}
	}
	
	unsigned int sbGrid = rlLoadShaderBuffer(visBufferSize * sizeof(uint32_t), tmpVisBuffer, RL_DYNAMIC_COPY);
	unsigned int sbData = rlLoadShaderBuffer(sizeof(GpuData), &gpuData, RL_DYNAMIC_COPY);
	
	//Setup compute shader data
	rlEnableShader(cellVisProgram);
	
	rlSetUniformSampler(rlGetLocationUniform(cellVisProgram, "cubeIdxTex"), 0);
	rlSetUniformSampler(rlGetLocationUniform(cellVisProgram, "cubeDepthTex"), 1);
					
	rlActiveTextureSlot(0);
	rlEnableTexture(target.texture.id);
	rlActiveTextureSlot(1);
	rlEnableTexture(target.depth.id);
	
	rlBindShaderBuffer(sbGrid, 1);
	rlBindShaderBuffer(sbData, 2);
	
	rlEnableShader(renderShader.id);
	
	rlBindShaderBuffer(sbData, 2);
					
	rlSetVertexAttribute(renderShader.locs[SHADER_LOC_VERTEX_POSITION], 3, RL_FLOAT, 0, 0, 0);
	rlEnableVertexAttribute(renderShader.locs[SHADER_LOC_VERTEX_POSITION]);
					
	rlDisableShader();
	
	//Setup render
	size_t cubemapsCount = 0;
	size_t cubemapsPerCell = 8;
	
	rlDisableColorBlend();
	rlDisableBackfaceCulling();
	rlEnableDepthTest();// Enable DEPTH_TEST for 3D
	rlClearColor(0, 0, 0, 255);
	//glEnable(0x9346); //CONSERVATIVE_RASTERIZATION_NV
	
	Vector3 cellSize = {
		(db->max[0] - db->min[0]) / db->cellsX,
		(db->max[1] - db->min[1]) / db->cellsY,
		(db->max[2] - db->min[2]) / db->cellsZ
	};
	
	int prevPercentage = -1;
	clock_t start = clock();
	
	rlMatrixMode(RL_PROJECTION);
    rlPushMatrix();
	rlLoadIdentity();               // Reset current matrix (projection)

	//Setup perspective projection
	const float top = RL_CULL_DISTANCE_NEAR * tan(90*0.5*DEG2RAD);
	const float right = top;

	rlFrustum(-right, right, -top, top, RL_CULL_DISTANCE_NEAR, RL_CULL_DISTANCE_FAR);
	Matrix matProj = rlGetMatrixProjection();
	
    rlPopMatrix();
	
	rlMatrixMode(RL_MODELVIEW);     // Switch back to modelview matrix
	
	rlEnableFramebuffer(target.id);
	rlViewport(0, 0, size * 3, size * 2);
	
	glEnable(12288); //GL_CLIP_DISTANCE0
	glEnable(12288 + 1); //GL_CLIP_DISTANCE1
	glEnable(12288 + 2); //GL_CLIP_DISTANCE2
	glEnable(12288 + 3); //GL_CLIP_DISTANCE3
	
	for(size_t z=0; z<db->cellsZ; z++) {
		for(size_t y=0; y<db->cellsY; y++) {
			for(size_t x=0; x<db->cellsX; x++) {
				
				for(size_t cellCubemap = 0; cellCubemap < cubemapsPerCell; cellCubemap++) {
					Vector3 camPos = {
						db->min[0] + cellSize.x * (x + randf()),
						db->min[1] + cellSize.y * (y + randf()),
						db->min[2] + cellSize.z * (z + randf()),
					};
					
					rlEnableShader(renderShader.id);
					rlClearScreenBuffers();
					
					for(int viewport = 0; viewport < 6; viewport++) {
						Vector3 target = (Vector3) {(viewport >> 1) == 0, (viewport >> 1) == 1, (viewport >> 1) == 2};
						target = Vector3Add(camPos, Vector3Scale(target, (viewport & 1) ? -1 : 1));
						
						Vector3 upVec = (Vector3) {0, -1, 0};
						if((viewport >> 1) == 1) upVec = (Vector3) {0, 0, (viewport & 1) ? -1 : 1};

						//Setup Camera view
						Matrix matView = MatrixLookAt(camPos, target, upVec);
						
						Matrix* vpMat = (Matrix*) (gpuData.matVP + viewport * 16);
						Matrix* vpMatInv = (Matrix*) (gpuData.matVPInv + viewport * 16);
						
						*vpMat = MatrixMultiply(matView, matProj);
						*vpMatInv = MatrixInvert(*vpMat);
					}
					
					rlUpdateShaderBuffer(sbData, &gpuData, sizeof(Matrix) * 6 * 2, 0);
					
					for(int i=0; i<mdl->meshCount; i++) {
						int meshId = i + 1;
						rlSetUniform(meshIdUniform, &meshId, RL_SHADER_UNIFORM_INT, 1);
						
						Mesh* mesh = &mdl->meshes[i];
						
						if(!rlEnableVertexArray(mesh->vaoId)) {
							// Bind mesh VBO data: vertex position (shader-location = 0)
							rlEnableVertexBuffer(mesh->vboId[0]);

							if(mesh->indices != NULL) rlEnableVertexBufferElement(mesh->vboId[6]);
						}
						
						if(mesh->indices != NULL) rlDrawVertexArrayElementsInstanced(0, mesh->triangleCount * 3, 0, 6);
						else rlDrawVertexArrayInstanced(0, mesh->vertexCount, 6);
					}
					
					//Process cubemap results
					rlEnableShader(cellVisProgram);
					rlComputeShaderDispatch(size / 8, size / 8, 6);
					
					cubemapsCount++;
				}
				
				int percentage = cubemapsCount * 100 / (db->cellsX * db->cellsY * db->cellsZ * cubemapsPerCell);
				
				if(prevPercentage != percentage) {
					prevPercentage = percentage;
					printf("%d %%\n", percentage);
					SwapScreenBuffer();
				}
			}
		}
	}
	
	//Copy result to cpu
	rlReadShaderBuffer(sbGrid, tmpVisBuffer, visBufferSize * sizeof(uint32_t), 0);
	rlReadShaderBuffer(sbData, &gpuData, sizeof(GpuData), 0);
	
	float raysComputeTime = (float) (clock() - start) / CLOCKS_PER_SEC;
	printf("CUBEMAPS COMPUTE TIME %.2f sec\n", raysComputeTime);
	
	for(size_t z=0; z<db->cellsZ; z++) {
		for(size_t y=0; y<db->cellsY; y++) {
			for(size_t x=0; x<db->cellsX; x++) {
				
				PVSCell* cell = db->cellGrid + _getCellPtr(db, x, y, z);
				uint32_t* tmpVis = tmpVisBuffer + _getCellPtr(db, x, y, z) * gpuData.dataPerCell;
				
				for(int i = 0; i < gpuData.dataPerCell; i++) {
					cell->visMesh[i] |= tmpVis[i];
				}
			}
		}
	}
	
	//just kiddin ;)
	rlDisableFramebuffer();
	rlDisableShader();
	
	rlDisableVertexArray();
	rlDisableVertexBuffer();
	rlDisableVertexBufferElement();
					
	rlUnloadShaderProgram(cellVisProgram);
	rlUnloadShaderBuffer(sbGrid);
	rlUnloadShaderBuffer(sbData);
	free(tmpVisBuffer);
	
	rlViewport(0, 0, GetRenderWidth(), GetRenderHeight());
	rlEnableColorBlend();
	rlEnableBackfaceCulling();
	rlDisableDepthTest();
	
	glDisable(12288); //GL_CLIP_DISTANCE0
	glDisable(12288 + 1); //GL_CLIP_DISTANCE0
	glDisable(12288 + 2); //GL_CLIP_DISTANCE0
	glDisable(12288 + 3); //GL_CLIP_DISTANCE0
	
	while(!IsKeyPressed(KEY_ENTER)) {
		BeginDrawing();
		ClearBackground(RAYWHITE);
		DrawTextureEx(target.texture, (Vector2) {0, 0}, 0, 1, WHITE);
		
		/*for(int view = 0; view < 6; view++) {
			printf("%d %d %d %d %d\n", view, (int) gpuData.outAABB[view * 4], (int) gpuData.outAABB[view * 4 + 1], (int) gpuData.outAABB[view * 4 + 2], (int) gpuData.outAABB[view * 4 + 3]);
			DrawRectangleLines(
				gpuData.outAABB[view * 4] + (view % 3) * size, 
				gpuData.outAABB[view * 4 + 1] + (view / 3) * size, 
				gpuData.outAABB[view * 4 + 2] - gpuData.outAABB[view * 4], 
				gpuData.outAABB[view * 4 + 3] - gpuData.outAABB[view * 4 + 1], 
				GREEN
			);
		}*/
		
		EndDrawing();
	}
	
	rlUnloadFramebuffer(target.id);
	UnloadShader(renderShader);
	mat.shader = (Shader) {0};
	UnloadMaterial(mat);
	
	return cubemapsCount;
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
		res.visMeshesCount += (res.visible[i / 32] & (1 << (i % 32))) != 0;
	}
	
	return res;
}

size_t _cellsVisMeshesCount(PVSdb* db, uint32_t* tmpVisMeshes, size_t visMeshesSize, size_t* min, size_t* max) {
	
	memset(tmpVisMeshes, 0, visMeshesSize * sizeof(uint32_t));
	
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
		
		for(int t=0; t<32; t++) visCount += (cellVisible >> t) & 1;
	}
	
	return visCount;
}

double _cellsError(PVSdb* db, uint32_t* tmpVisMeshes, size_t visMeshesSize, size_t* min, size_t* max, bool useFacesArea) {
	//Calculate error in cell min-max
	
	memset(tmpVisMeshes, 0, visMeshesSize * sizeof(uint32_t));
	
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
					int meshError = (tmpVisMeshes[i / 32] ^ cell->visMesh[i / 32]) >> (i&31);
					meshError &= 1;
					
					error += meshError * db->meshesData[i].facesArea * db->meshesData[i].faces;//(useFacesArea ? mdl->meshes[i].facesArea : mdl->meshes[i].facesCount);
				}
				//for(size_t i=0; i<visMeshesSize; i++) {
				//	int cellError = tmpVisMeshes[i] ^ cell->visMesh[i];
				//	
				//	for(int t=0; t<8; t++) error += (cellError >> t) & 1;
				//}
				
			}
		}
	}
	
	return error;
}

inline void _cellFindSplit(PVSdb* db, uint32_t* tmpVisMeshes, size_t visMeshesSize, PVSCell* cell, int* splitAxis, double* minSplitError, double* minMaxSplitError, size_t* splitPos, bool useFacesArea) {
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
			double err1 = _cellsError(db, tmpVisMeshes, visMeshesSize, splitMin, splitMax, useFacesArea);
			
			splitMin[axis] = pos;
			splitMax[axis] = cell->max[axis];
			double err2 = _cellsError(db, tmpVisMeshes, visMeshesSize, splitMin, splitMax, useFacesArea);
			
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

cvector_vector_type(PVSCell) pvsJoinCells(PVSdb* db, size_t cellsCount, float minError, bool splitClever) {
	
	cvector_vector_type(PVSCell) cells = NULL;
	
	size_t visMeshesSize = (db->meshes / 32) + ((db->meshes % 32) > 0 ? 1 : 0);
	uint32_t* tmpVisMeshes = malloc(visMeshesSize * sizeof(uint32_t));
	assert(tmpVisMeshes);
	
	PVSCell firstCell = {0};
	firstCell.max[0] = db->cellsX;
	firstCell.max[1] = db->cellsY;
	firstCell.max[2] = db->cellsZ;
	firstCell.error = -1;
	
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
			
			double err = cell->error;
			if(err < 0) {
				err = _cellsError(db, tmpVisMeshes, visMeshesSize, cell->min, cell->max, useFacesArea);
				cell->error = err;
			}
			
			if(err < splitCellError) continue;
			if(err < minError) continue;
			
			//Find an axis aligned split plane
			if(splitClever) {
				int splitAxis = -1;
				double minSplitError;
				double minMaxSplitError;
				size_t splitPos;
				_cellFindSplit(db, tmpVisMeshes, visMeshesSize, cell, &splitAxis, &minSplitError, &minMaxSplitError, &splitPos, useFacesArea);
				
				if(minSplitError > err) printf("wtf\n");
				
				err = err - minMaxSplitError;
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
		_cellFindSplit(db, tmpVisMeshes, visMeshesSize, splitCell, &splitAxis, &minSplitError, &minMaxSplitError, &splitPos, useFacesArea);
		
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
		
		splitCell->error = -1;
		cell2.error = -1;
		
		cvector_push_back(cells, cell2);
		
		if(_cellsVisMeshesCount(db, tmpVisMeshes, visMeshesSize, cell2.min, cell2.max) == 0) cvector_pop_back(cells);
		if(_cellsVisMeshesCount(db, tmpVisMeshes, visMeshesSize, splitCell->min, splitCell->max) == 0) cvector_erase(cells, splitCellId);
	}
	
	//Compute visibility for every cell
	for(size_t i=0; i<cvector_size(cells); i++) {
		PVSCell* cell = cells + i;
		
		cell->visMesh = calloc(visMeshesSize, sizeof(uint32_t));
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