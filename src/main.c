#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <time.h>

#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"
#include "cvector.h"

#include "pvs_model.c"
#include "pvs_raycast.c"
#include "pvs.c"

int main(void) {
	InitWindow(1280, 720, "pvs test");
	SetWindowState(FLAG_WINDOW_RESIZABLE);
	
    SetTargetFPS(60);
	
	Model mdl = LoadModel("test2.glb");
	Texture2D tex = LoadTexture("tex2.png");
	GenTextureMipmaps(&tex);
	//SetTextureFilter(tex, TEXTURE_FILTER_POINT);
	rlTextureParameters(tex.id, RL_TEXTURE_MIN_FILTER, RL_TEXTURE_FILTER_NEAREST_MIP_LINEAR);
	rlTextureParameters(tex.id, RL_TEXTURE_MAG_FILTER, RL_TEXTURE_FILTER_NEAREST_MIP_LINEAR);
	//rlTextureParameters(tex.id, RL_TEXTURE_FILTER_ANISOTROPIC, 16);
	mdl.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = tex;
	
	/*int size;
	char* data = LoadFileData("out.lrb", &size);*/
	
	Model cube = LoadModelFromMesh(GenMeshCube(1, 1, 1));
	
	#define GLSL_VERSION 330
	
	Shader shader = LoadShader(TextFormat("shaders/shader-glsl%i.vs", GLSL_VERSION),
                               TextFormat("shaders/shader-glsl%i.fs", GLSL_VERSION));
	shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(shader, "viewPos");
	
	mdl.materials[0].shader = shader;
	cube.materials[0].shader = shader;
	
	Vector3 camPos = (Vector3) {0};
	Vector3 pvsCamPos = (Vector3) {0};
	float rotX = 4.100000, rotY = 58.299999;
	
	int t = 0;
	
	PVSRay ray = {0};
	PVSRayHit rayHit = {0};
	
	//todo unload!!!
	bool usePVS = true;
	bool lockPVS = false;
	bool useCellData = false;
	bool hideWorld = false;
	
	PVSdb pvsDB = {0};
	PVSModel pvsMdl = {NULL};
	pvsLoadObj(&pvsMdl, "test2.obj");
	
	pvsInitDB(&pvsDB, pvsMdl, 3);
	size_t totalRays = 0;
	double raysComputeTime = 0;
	
	cvector_vector_type(PVSCell) joinedCells = NULL;
	
	PVSResult pvsRes = {0};
	int cellsAverageSize = 10;
	int cellsDebug = -1;
	
	//cvector_vector_type(float) tracedRays = NULL;
	
	float posAnim[6*3] = {
		4.430673, 2.525841, 11.464261,
		4.401867, 2.180269, -21.533922,
		25.260057, 3.054857, -21.514534,
		25.945866, 2.275901, 1.857101,
		47.035301, 2.374130, 0.804444,
		47.749432, 2.810467, -35.547024
	};
	
	int posAnimFrames = 6;
	int posAnimId = -1;
	clock_t posAnimStart = 0;
	
	printf("THREADS: %d\n", omp_get_max_threads());
	
	float dpos[3] = {0};
	float dhit[3] = {0};
	
	while(!WindowShouldClose()) {
		
		if(IsWindowFocused()) {
			Vector2 delta = GetMouseDelta();
			rotY -= delta.x / GetRenderHeight() * 90;
			rotX += delta.y / GetRenderHeight() * 90;
			if(rotX < -89) rotX = -89;
			if(rotX > 89) rotX = 89;
			
			SetMousePosition(GetRenderWidth()/2, GetRenderHeight()/2);
		}
		
		//Movement
		float speed = 0.12 * (IsKeyDown(KEY_LEFT_SHIFT) ? 3 : 1) * (IsKeyDown(KEY_LEFT_CONTROL) ? 0.333 : 1);
		
		if(IsKeyDown(KEY_W)) {
			Vector3 vec = (Vector3) {0, 0, speed};
			
			vec = Vector3RotateByAxisAngle(vec, (Vector3) {1, 0, 0}, rotX * M_PI / 180);
			vec = Vector3RotateByAxisAngle(vec, (Vector3) {0, 1, 0}, rotY * M_PI / 180);
			
			camPos = Vector3Add(camPos, vec);
		}
		if(IsKeyDown(KEY_S)) {
			Vector3 vec = (Vector3) {0, 0, -speed};
			
			vec = Vector3RotateByAxisAngle(vec, (Vector3) {1, 0, 0}, rotX * M_PI / 180);
			vec = Vector3RotateByAxisAngle(vec, (Vector3) {0, 1, 0}, rotY * M_PI / 180);
			
			camPos = Vector3Add(camPos, vec);
		}
		if(IsKeyDown(KEY_A)) {
			Vector3 vec = (Vector3) {speed, 0, 0};
			
			vec = Vector3RotateByAxisAngle(vec, (Vector3) {0, 1, 0}, rotY * M_PI / 180);
			
			camPos = Vector3Add(camPos, vec);
		}
		if(IsKeyDown(KEY_D)) {
			Vector3 vec = (Vector3) {-speed, 0, 0};
			
			vec = Vector3RotateByAxisAngle(vec, (Vector3) {0, 1, 0}, rotY * M_PI / 180);
			
			camPos = Vector3Add(camPos, vec);
		}
		
		//Animation
		if(IsKeyPressed(KEY_E)) {
			if(IsKeyDown(KEY_LEFT_SHIFT)) {
				posAnimId = 0;
				posAnimStart = clock();
			} else {
				printf("%f, %f, %f,\n", camPos.x, camPos.y, camPos.z);
			}
		}
		
		//PVS Contols
		if(IsKeyPressed(KEY_Q)) {
			Vector3 vec = (Vector3) {0, 0, speed};
			
			vec = Vector3RotateByAxisAngle(vec, (Vector3) {1, 0, 0}, rotX * M_PI / 180);
			vec = Vector3RotateByAxisAngle(vec, (Vector3) {0, 1, 0}, rotY * M_PI / 180);
			
			ray.start[0] = camPos.x;
			ray.start[1] = camPos.y;
			ray.start[2] = camPos.z;
			
			ray.dir[0] = vec.x;
			ray.dir[1] = vec.y;
			ray.dir[2] = vec.z;
			
			ray.length = FLT_MAX;
			
			rayHit.hit = false;
			
			for(int i=0; i<cvector_size(pvsMdl.meshes); i++) {
				rayCast(&rayHit, pvsMdl.meshes + i, &ray, -1);
			}
			
			if(rayHit.hit) printf("hit! ");
			else printf("miss! ");
			
			printf("%f %f %f %f\n", ray.start[0], ray.start[1], ray.start[2], rayHit.hitDistance);
		}
		
		
		if(IsKeyPressed(KEY_L)) lockPVS = !lockPVS;
		if(IsKeyPressed(KEY_P)) usePVS = !usePVS;
		if(IsKeyPressed(KEY_H)) hideWorld = !hideWorld;
		
		if(IsKeyPressed(KEY_X)) {
			clock_t start = clock();
			totalRays += pvsCompute(&pvsDB, pvsMdl, 32, 8, dpos, dhit);
			//totalRays += pvsCompute(&pvsDB, pvsMdl, 6.54f, 8, dpos, dhit);
			//totalRays += pvsCompute(&pvsDB, pvsMdl, 48, 99999, dpos, dhit);
			//totalRays += pvsCompute(&pvsDB, pvsMdl, 48, 8, dpos, dhit);
			//totalRays += pvsCompute(&pvsDB, pvsMdl, 48/5.78, 99999, dpos, dhit);
			raysComputeTime += (double) (clock() - start) / CLOCKS_PER_SEC;
		}
		
		if((IsKeyPressed(KEY_RIGHT) || IsKeyPressedRepeat(KEY_RIGHT)) && !IsKeyDown(KEY_LEFT_SHIFT)) cellsAverageSize += cellsAverageSize>=5?1:1;
		if((IsKeyPressed(KEY_LEFT) || IsKeyPressedRepeat(KEY_LEFT)) && !IsKeyDown(KEY_LEFT_SHIFT) && cellsAverageSize > 1) cellsAverageSize -= cellsAverageSize>5?1:1;
		
		if((IsKeyPressed(KEY_RIGHT) || IsKeyPressedRepeat(KEY_RIGHT)) && IsKeyDown(KEY_LEFT_SHIFT)) cellsDebug += (cellsDebug == mdl.meshCount-1) ? 0 : 1;
		if((IsKeyPressed(KEY_LEFT) || IsKeyPressedRepeat(KEY_LEFT)) && IsKeyDown(KEY_LEFT_SHIFT)) cellsDebug -= (cellsDebug == -1) ? 0 : 1;
		
		if(IsKeyPressed(KEY_C)) useCellData = !useCellData;
		
		if(IsKeyPressed(KEY_J)) {
			
			size_t cellsCount = (pvsDB.max[0] - pvsDB.min[0]) * (pvsDB.max[1] - pvsDB.min[1]) * (pvsDB.max[2] - pvsDB.min[2]);
			//cellsCount *= pvsDB.cellsX * pvsDB.cellsY * pvsDB.cellsZ;
			//cellsCount /= pvsDB.cellsX * pvsDB.cellsY * pvsDB.cellsZ;
			cellsCount /= cellsAverageSize * cellsAverageSize * cellsAverageSize;
			
			cvector_free(joinedCells);
			joinedCells = pvsJoinCells(&pvsDB, pvsMdl, cellsCount, 0, IsKeyDown(KEY_LEFT_SHIFT));
		}
		
		//Update animation
		if(posAnimId != -1) {
			
			while(1) {
				
				if(posAnimId >= posAnimFrames-1) {
					posAnimId = -1;
					break;
					
				} else {
					clock_t curTime = clock();
					clock_t spent = curTime - posAnimStart;
					
					Vector3 startPos = (Vector3) {posAnim[posAnimId * 3], posAnim[posAnimId * 3 + 1], posAnim[posAnimId * 3 + 2]};
					Vector3 endPos = (Vector3) {posAnim[posAnimId * 3 + 3], posAnim[posAnimId * 3 + 1 + 3], posAnim[posAnimId * 3 + 2 + 3]};
					
					clock_t frameLength = (clock_t) (Vector3Length(Vector3Subtract(endPos, startPos)) * CLOCKS_PER_SEC / 8);
					
					if(spent > frameLength) {
						posAnimStart += frameLength;
						posAnimId++;
						pvsCamPos = endPos;
						
					} else {
						pvsCamPos = Vector3Lerp(startPos, endPos, (double)spent / frameLength);
						break;
					}
				}
			}
			
		}
		
		//Draw
		SetShaderValue(shader, shader.locs[SHADER_LOC_VECTOR_VIEW], &camPos, SHADER_UNIFORM_VEC3);
		
		if(!lockPVS) pvsCamPos = camPos;
		pvsRes = pvsGetGridVisibility(&pvsDB, pvsCamPos.x, pvsCamPos.y, pvsCamPos.z);
		
		BeginDrawing();
		ClearBackground(RAYWHITE);
		
		//Setup camera
		Vector3 target = (Vector3) {0, 0, 1};
			
		target = Vector3RotateByAxisAngle(target, (Vector3) {1, 0, 0}, rotX * M_PI / 180);
		target = Vector3RotateByAxisAngle(target, (Vector3) {0, 1, 0}, rotY * M_PI / 180);
		target = Vector3Add(camPos, target);
		
		Camera3D cam = (Camera3D) {
			camPos, 
			target,
			(Vector3) {0, 1, 0},             // Camera up vector (rotation over its axis)
			70,             // Camera field-of-view aperture in Y (degrees) in perspective, used as near plane width in orthographic
			CAMERA_PERSPECTIVE         // Camera projection: CAMERA_PERSPECTIVE or CAMERA_ORTHOGRAPHIC
		};
		
		BeginMode3D(cam);
		
		//Camera
		DrawModelEx(
			cube, 
			(Vector3) pvsCamPos, 
			(Vector3) {0, 1, 0}, 
			0,
			(Vector3) {0.5, 0.5, 0.5}, 
			(Color) {255, 255, 0, 255}
		);
		
		//printf("%d\n", janua_res.model_ids_count);
		//printf("%d\n", mdl.meshCount);
		
		float avgCellSize = 0;
		
		if(useCellData) {
			pvsRes.visMeshesCount = 0;
			
			for(size_t i=0; i<cvector_size(joinedCells); i++) {
				Vector3 min = (Vector3) {joinedCells[i].min[0], joinedCells[i].min[1], joinedCells[i].min[2]};
				Vector3 max = (Vector3) {joinedCells[i].max[0], joinedCells[i].max[1], joinedCells[i].max[2]};
				
				Vector3 gridSize = (Vector3) {pvsDB.cellsX, pvsDB.cellsY, pvsDB.cellsZ};
				Vector3 dbSize = (Vector3) {pvsDB.max[0] - pvsDB.min[0], pvsDB.max[1] - pvsDB.min[1], pvsDB.max[2] - pvsDB.min[2]};
				Vector3 dbPos = (Vector3) {pvsDB.min[0], pvsDB.min[1], pvsDB.min[2]};
			
				min = Vector3Multiply(min, dbSize);
				min = Vector3Divide(min, gridSize);
				min = Vector3Add(min, dbPos);
				
				max = Vector3Multiply(max, dbSize);
				max = Vector3Divide(max, gridSize);
				max = Vector3Add(max, dbPos);
				
				Vector3 size = Vector3Subtract(max, min);
				
				avgCellSize += size.x * size.y * size.z;
				
				if(
					pvsCamPos.x >= min.x && pvsCamPos.x <= max.x && 
					pvsCamPos.y >= min.y && pvsCamPos.y <= max.y &&
					pvsCamPos.z >= min.z && pvsCamPos.z <= max.z) {
					
					pvsRes.meshesCount = pvsDB.meshes;
					pvsRes.visible = joinedCells[i].visMesh;
					
					for(size_t t=0; t<pvsRes.meshesCount; t++) {
						pvsRes.visMeshesCount += (pvsRes.visible[t / 8] & (1 << (t % 8))) != 0;
					}
				}
			}
			
			avgCellSize /= cvector_size(joinedCells);
		}
		
		rlSetLineWidth(5);
		if(!hideWorld) for(int i=0; i<mdl.meshCount; i++) {
			
			bool visible = !usePVS || (pvsRes.visMeshesCount > 0 && (pvsRes.visible[i / 8] & (1 << (i % 8))) != 0);
			int matId = 0;//
			matId = mdl.meshMaterial[i];
			
			rlSetCullFace(RL_CULL_FACE_FRONT);
			Color colBack = visible?((Color){64,64,64,255}):((Color){64,0,0,255});
			if(cellsDebug == i) colBack = (Color){0,64,0,255};
			
			mdl.materials[matId].maps[MATERIAL_MAP_ALBEDO].color = colBack;//visible?RED:((Color){64,0,0,255});
			DrawMesh(mdl.meshes[i], mdl.materials[matId], mdl.transform);
			
			rlSetCullFace(RL_CULL_FACE_BACK);
			Color col = visible?WHITE:RED;
			if(cellsDebug == i) col = (Color){0,255,0,255};
			
			mdl.materials[matId].maps[MATERIAL_MAP_ALBEDO].color = col;
			DrawMesh(mdl.meshes[i], mdl.materials[matId], mdl.transform);
			
			Matrix transTmp = (Matrix) {
				1.05f, 0, 0, -camPos.x * 0.05f,
				0, 1.05f, 0, -camPos.y * 0.05f,
				0, 0, 1.05f, -camPos.z * 0.05f,
				0, 0, 0, 1
			};
			
			rlEnableWireMode();
			mdl.materials[matId].maps[MATERIAL_MAP_ALBEDO].color = colBack;
			DrawMesh(mdl.meshes[i], mdl.materials[matId], MatrixMultiply(mdl.transform, transTmp));
			rlDisableWireMode();
			
			if(IsKeyDown(KEY_LEFT_CONTROL)) {
				float* min = pvsMdl.meshes[i].min;
				float* max = pvsMdl.meshes[i].max;
				
				DrawCubeWiresV((Vector3) {(min[0] + max[0]) * 0.5f, (min[1] + max[1]) * 0.5f, (min[2] + max[2]) * 0.5f}, (Vector3) {max[0] - min[0], max[1] - min[1], max[2] - min[2]}, RED);
			}
		}
		
		rlSetLineWidth(1);
		
		DrawLine3D(
			(Vector3) {dpos[0], dpos[1], dpos[2]},
			(Vector3) {dhit[0], dhit[1], dhit[2]},
			BLUE
		);
		
		DrawModelEx(
			cube, 
			(Vector3) {dpos[0], dpos[1], dpos[2]},
			(Vector3) {0, 1, 0}, 
			0,
			(Vector3) {0.5, 0.5, 0.5}, 
			(Color) {255, 255, 0, 255}
		);
	
		if(rayHit.hit) {
			Color col = (Color) {rayHit.hitBackFace * 255, 255 - rayHit.hitBackFace * 255, 0, 255};
			
			DrawLine3D(
				(Vector3) {ray.start[0], ray.start[1], ray.start[2]},
				(Vector3) {
					ray.start[0] + ray.dir[0] * rayHit.hitDistance, 
					ray.start[1] + ray.dir[1] * rayHit.hitDistance, 
					ray.start[2] + ray.dir[2] * rayHit.hitDistance
				},
				col
			);
			
			DrawModelEx(
				cube, 
				(Vector3) {
					ray.start[0] + ray.dir[0] * rayHit.hitDistance, 
					ray.start[1] + ray.dir[1] * rayHit.hitDistance, 
					ray.start[2] + ray.dir[2] * rayHit.hitDistance
				},
				(Vector3) {0, 1, 0}, 
				0,
				(Vector3) {0.1, 0.1, 0.1}, 
				col
			);
		}
		
		//rlDisableDepthMask();
		for(size_t i=0; i<cvector_size(joinedCells); i++) {
			
			Vector3 min = (Vector3) {joinedCells[i].min[0], joinedCells[i].min[1], joinedCells[i].min[2]};
			Vector3 max = (Vector3) {joinedCells[i].max[0], joinedCells[i].max[1], joinedCells[i].max[2]};
			
			Vector3 gridSize = (Vector3) {pvsDB.cellsX, pvsDB.cellsY, pvsDB.cellsZ};
			Vector3 dbSize = (Vector3) {pvsDB.max[0] - pvsDB.min[0], pvsDB.max[1] - pvsDB.min[1], pvsDB.max[2] - pvsDB.min[2]};
			Vector3 dbPos = (Vector3) {pvsDB.min[0], pvsDB.min[1], pvsDB.min[2]};
			
			Color col = (Color) {min.x * 255 / (gridSize.x - 1), 0, min.z * 255 / (gridSize.z - 1), 255};
			//col = (Color) {255, 0, 0, 255};
			
			min = Vector3Multiply(min, dbSize);
			min = Vector3Divide(min, gridSize);
			min = Vector3Add(min, dbPos);
			
			max = Vector3Multiply(max, dbSize);
			max = Vector3Divide(max, gridSize);
			max = Vector3Add(max, dbPos);
			
			DrawCubeWires(
				Vector3Scale(Vector3Add(min, max), 0.5),
				max.x - min.x - 0.5, max.y - min.y - 0.5, max.z - min.z - 0.5,
				col
			);
				
			if(cellsDebug != -1 && ((joinedCells[i].visMesh[cellsDebug/8] >> (cellsDebug&7)) & 1) == 0) {
				DrawModelEx(
					cube, 
					Vector3Scale(Vector3Add(min, max), 0.5),
					(Vector3) {0, 1, 0}, 
					0,
					(Vector3) {max.x - min.x - 0.5, max.y - min.y - 0.5, max.z - min.z - 0.5}, 
					RED
				);
			}
		}
		/*rlDrawRenderBatchActive();
		rlEnableDepthMask();*/
		
		/*for(size_t i=0; i<cvector_size(tracedRays); i+=3) {
			DrawCubeNormals(
				(Vector3) {
					tracedRays[i], 
					tracedRays[i + 1], 
					tracedRays[i + 2]
				},
				0.02, 0.02, 0.02,
				(Color) {255, 0, 0, 128}
			);
		}*/
		
		rlDrawRenderBatchActive();
		
		EndMode3D();
		
		{
			Font fnt = GetFontDefault();
			int line = 0;
			
			DrawText(TextFormat("(P)VS on: %s", usePVS ? "on" : "off"), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("(L)ock PVS position: %s", lockPVS ? "on" : "off"), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("Visible meshes: %ld", pvsRes.visMeshesCount), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			line++;
			
			DrawText(TextFormat("PVS calc (X) time: %.2f sec", raysComputeTime), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("Total rays: %ld", totalRays), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("Rays per sec: %d", (int) (totalRays / raysComputeTime)), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			line++;
			
			DrawText(TextFormat("Cells av. size (arrows) after (J)oin: %d m", cellsAverageSize), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("Avg size: %.2f m", (float) pow(avgCellSize, 1 / 3.0)), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("Cells count: %d", cvector_size(joinedCells)), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("Use (C)ell data: %s", useCellData ? "on" : "off"), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("(H)ide world: %s", hideWorld ? "on" : "off"), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
		}

        EndDrawing();
		t++;
	}
	
	pvsUnloadDB(&pvsDB);

	UnloadTexture(tex);
	UnloadModel(mdl);
	UnloadModel(cube);
	UnloadShader(shader);
	
    CloseWindow();
	
	return 0;
}
