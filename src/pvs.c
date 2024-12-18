#ifndef _pvs_c
#define _pvs_c

#include <stdlib.h>
#include <stdint.h>

#include "cvector.h"

#include "raylib.h"
#include "raymath.h"

typedef struct {
	Vector3 min, max;
	
	float facesArea;
	int faces;
} PVSMeshData;

typedef struct {
	Vector3 min, max;
	
	int meshCount;
	PVSMeshData* meshData;
} PVSModelData;

typedef struct {
	int gridSize[3];
	Vector3 min, max;
	
	int meshCount;
	
	int intsPerCell;
	uint32_t* cells;
} PVSdb;

typedef struct {
	int meshCount;
	int visMeshCount;
	
	uint32_t* visible;
} PVSResult;

void pvsInitModelData(PVSModelData* mdlData, Model* mdl);
void pvsUnloadModelData(PVSModelData* mdlData);

void pvsInitDB(PVSdb* db, PVSModelData* mdlData, float cellSize);
void pvsUnloadDB(PVSdb* db);

size_t pvsCompute(PVSdb* db, Model* mdl, Vector3 camPos);
PVSResult pvsGetVisData(PVSdb* db, Vector3 camPos);

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
	return x + y * db->gridSize[0] + z * db->gridSize[0] * db->gridSize[1];
}

void pvsInitModelData(PVSModelData* mdlData, Model* mdl) {
	//Allocate memory for meshes data
	mdlData->meshCount = mdl->meshCount;
	mdlData->meshData = malloc(sizeof(PVSMeshData) * mdl->meshCount);
	assert(mdlData->meshData);
	
	//Preparations for model AABB calculation
	mdlData->min = (Vector3) {FLT_MAX, FLT_MAX, FLT_MAX};
	mdlData->max = (Vector3) {-FLT_MAX, -FLT_MAX, -FLT_MAX};
	
	//Compute data for each mesh
	for(int i = 0; i < mdl->meshCount; i++) {
		PVSMeshData* meshData = mdlData->meshData + i;
		const Mesh mesh = mdl->meshes[i];
		
		BoundingBox meshAABB = GetMeshBoundingBox(mesh);
		
		//Copy AABB to mesh data
		meshData->min = meshAABB.min;
		meshData->max = meshAABB.max;
		
		//Update model AABB
		mdlData->min = Vector3Min(mdlData->min, meshAABB.min);
		mdlData->max = Vector3Max(mdlData->max, meshAABB.max);
		
		//Calculate faces count and summ of faces area
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
			float faceArea = sqrtf(S * (S-distAB) * (S-distBC) * (S-distAC) / 2);
			
			meshData->facesArea += faceArea;
		}
	}
}

void pvsUnloadModelData(PVSModelData* mdlData) {
	free(mdlData->meshData);
}

void pvsInitDB(PVSdb* db, PVSModelData* mdlData, float cellSize) {
	//Copy AABB from model data
	db->min = mdlData->min;
	db->max = mdlData->max;
	
	//Calculate grid size
	Vector3 gridSize = Vector3Subtract(db->max, db->min);
	gridSize = Vector3Scale(gridSize, 1.0f / cellSize);
	
	db->gridSize[0] = (int) roundf(gridSize.x);
	db->gridSize[1] = (int) roundf(gridSize.y);
	db->gridSize[2] = (int) roundf(gridSize.z);
	
	if(db->gridSize[0] < 1) db->gridSize[0] = 1;
	if(db->gridSize[1] < 1) db->gridSize[1] = 1;
	if(db->gridSize[2] < 1) db->gridSize[2] = 1;
	
	//Allocate visibility cells
	db->meshCount = mdlData->meshCount;
	db->intsPerCell = (db->meshCount / 32) + ((db->meshCount % 32) > 0 ? 1 : 0);
	
	db->cells = calloc(db->gridSize[0] * db->gridSize[1] * db->gridSize[2] * db->intsPerCell, sizeof(uint32_t));
	assert(db->cells);
}

void pvsUnloadDB(PVSdb* db) {
	free(db->cells);
}

inline float randf() {
	return (float) rand() / RAND_MAX;
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
	
	gpuData.gridSize[0] = db->gridSize[0];
	gpuData.gridSize[1] = db->gridSize[1];
	gpuData.gridSize[2] = db->gridSize[2];
	
	gpuData.aabbMin[0] = db->min.x;
	gpuData.aabbMin[1] = db->min.y;
	gpuData.aabbMin[2] = db->min.z;
	
	gpuData.aabbMax[0] = db->max.x;
	gpuData.aabbMax[1] = db->max.y;
	gpuData.aabbMax[2] = db->max.z;
	
	/*gpuData.camPos[0] = (camPos.x - db->min[0]) * db->cellsX / (db->max[0] - db->min[0]);
	gpuData.camPos[1] = (camPos.y - db->min[1]) * db->cellsY / (db->max[1] - db->min[1]);
	gpuData.camPos[2] = (camPos.z - db->min[2]) * db->cellsZ / (db->max[2] - db->min[2]);*/
	
	gpuData.dataPerCell = db->intsPerCell;
	
	size_t visBufferSize = db->gridSize[0] * db->gridSize[1] * db->gridSize[2] * gpuData.dataPerCell;
	uint32_t* tmpVisBuffer = malloc(visBufferSize * sizeof(uint32_t));
	assert(tmpVisBuffer);
	
	memcpy(
		tmpVisBuffer, 
		db->cells, 
		gpuData.dataPerCell * sizeof(uint32_t) * db->gridSize[0] * db->gridSize[1] * db->gridSize[2]
	);
	
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
	
	Vector3 cellSize = Vector3Divide(
		Vector3Subtract(db->max, db->min), 
		(Vector3) {db->gridSize[0], db->gridSize[1], db->gridSize[2]}
	);
	
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
	
	for(size_t z=0; z<db->gridSize[2]; z++) {
		for(size_t y=0; y<db->gridSize[1]; y++) {
			for(size_t x=0; x<db->gridSize[0]; x++) {
				
				for(size_t cellCubemap = 0; cellCubemap < cubemapsPerCell; cellCubemap++) {
					Vector3 camPos = Vector3Add(
						db->min, 
						Vector3Multiply(cellSize, (Vector3) {x + randf(), y + randf(), z + randf()})
					);
					
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
				
				int percentage = cubemapsCount * 100 / (db->gridSize[0] * db->gridSize[1] * db->gridSize[2] * cubemapsPerCell);
				
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
	
	for(size_t i = 0; i < db->gridSize[0] * db->gridSize[1] * db->gridSize[2] * db->intsPerCell; i++) {
		db->cells[i] = tmpVisBuffer[i];
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

inline int _clamp(int x, int min, int max) {
	if(x < min) return min;
	else if(x > max) return max;
	else return x;
}

PVSResult pvsGetVisData(PVSdb* db, Vector3 camPos) {
	PVSResult res = {};
	
	Vector3 cellPos = Vector3Subtract(camPos, db->min);
	cellPos = Vector3Multiply(cellPos, (Vector3) {db->gridSize[0], db->gridSize[1], db->gridSize[2]});
	cellPos = Vector3Divide(cellPos, Vector3Subtract(db->max, db->min));
	
	int cellPosX = _clamp((int) cellPos.x, 0, db->gridSize[0]);
	int cellPosY = _clamp((int) cellPos.y, 0, db->gridSize[1]);
	int cellPosZ = _clamp((int) cellPos.z, 0, db->gridSize[2]);
	
	res.meshCount = db->meshCount;
	res.visMeshCount = 0;
	res.visible = db->cells + (cellPosX + cellPosY * db->gridSize[0] + cellPosZ * db->gridSize[0] * db->gridSize[1]);
	
	for(size_t i=0; i<res.meshCount; i++) {
		res.visMeshCount += (res.visible[i / 32] & (1 << (i % 32))) != 0;
	}
	
	return res;
}

#endif