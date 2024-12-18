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

#include "pvs.c"

typedef struct {
	float hitDist;
	int mdlIndex;
} MdlToRaycast;

int mdlHitDistanceCompare(const void *a, const void *b) {
	float mdl1 = *(float*)a;
	float mdl2 = *(float*)b;
	
	if(mdl1 < mdl2) return -1;
	else if(mdl1 > mdl2) return 1;
	else return 0;
}

int main(void) {
	SetConfigFlags(FLAG_MSAA_4X_HINT);
	InitWindow(1280, 720, "pvs test");
	SetWindowState(FLAG_WINDOW_RESIZABLE);
	
    SetTargetFPS(60);
	
	Model mdl = LoadModel("de_dust.glb");
	Texture2D tex = LoadTexture("44.png");
	GenTextureMipmaps(&tex);
	//SetTextureFilter(tex, TEXTURE_FILTER_POINT);
	rlTextureParameters(tex.id, RL_TEXTURE_MIN_FILTER, RL_TEXTURE_FILTER_NEAREST_MIP_LINEAR);
	rlTextureParameters(tex.id, RL_TEXTURE_MAG_FILTER, RL_TEXTURE_FILTER_NEAREST_MIP_LINEAR);
	//rlTextureParameters(tex.id, RL_TEXTURE_FILTER_ANISOTROPIC, 16);
	mdl.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = tex;
	
	int aniso = 0;
	
	Model cube = LoadModelFromMesh(GenMeshCube(1, 1, 1));
	
	Shader shader = LoadShader("shaders/fresnel.vs", "shaders/fresnel.fs");
	shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(shader, "viewPos");
	
	mdl.materials[0].shader = shader;
	cube.materials[0].shader = shader;
	
	Vector3 camPos = (Vector3) {0};
	Vector3 pvsCamPos = (Vector3) {0};
	float rotX = 4.100000, rotY = 58.299999;
	
	int t = 0;
	
	//todo unload!!!
	bool usePVS = true;
	bool lockPVS = false;
	bool useCellData = false;
	bool hideWorld = false;
	
	PVSModelData pvsMdlData = {0};
	PVSdb pvsDB = {0};
	
	pvsInitModelData(&pvsMdlData, &mdl);
	pvsInitDB(&pvsDB, &pvsMdlData, 3);
	
	size_t totalRays = 0;
	double raysComputeTime = 0;
	
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
	
	float dpos[3] = {0};
	float dhit[3] = {0};
	
	while(!WindowShouldClose()) {
		//Camera rotation
		if(IsWindowFocused()) {
			Vector2 delta = GetMouseDelta();
			rotY -= delta.x / GetRenderHeight() * 90;
			rotX += delta.y / GetRenderHeight() * 90;
			
			if(rotX < -89) rotX = -89;
			if(rotX > 89) rotX = 89;
			
			SetMousePosition(GetRenderWidth() / 2, GetRenderHeight() / 2);
		}
		
		//Camera movement
		Vector3 moveVec = (Vector3) {0};
		
		if(IsKeyDown(KEY_W)) moveVec.z += 1;
		if(IsKeyDown(KEY_S)) moveVec.z -= 1;
		if(IsKeyDown(KEY_A)) moveVec.x += 1;
		if(IsKeyDown(KEY_D)) moveVec.x -= 1;
		
		if(moveVec.x != 0 || moveVec.y != 0 || moveVec.z != 0) {
			float speed = 0.12f * (IsKeyDown(KEY_LEFT_SHIFT) ? 3 : 1) * (IsKeyDown(KEY_LEFT_CONTROL) ? 0.333f : 1);
			speed *= GetFrameTime() / (1 / 60.0f);
			
			moveVec = Vector3Scale(Vector3Normalize(moveVec), speed);
			
			moveVec = Vector3RotateByAxisAngle(moveVec, (Vector3) {1, 0, 0}, rotX * M_PI / 180);
			moveVec = Vector3RotateByAxisAngle(moveVec, (Vector3) {0, 1, 0}, rotY * M_PI / 180);
			
			camPos = Vector3Add(camPos, moveVec);
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
		
		
		if(IsKeyPressed(KEY_L)) lockPVS = !lockPVS;
		if(IsKeyPressed(KEY_P)) usePVS = !usePVS;
		if(IsKeyPressed(KEY_H)) hideWorld = !hideWorld;
		
		if(IsKeyPressed(KEY_X)) {
			clock_t start = clock();
			totalRays += pvsCompute(&pvsDB, &mdl, camPos);
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
		
		if(IsKeyPressed(KEY_F)) {
			aniso ^= 1;
			rlTextureParameters(tex.id, RL_TEXTURE_FILTER_ANISOTROPIC, aniso * 16);
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
		pvsRes = pvsGetVisData(&pvsDB, pvsCamPos);
		
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
		
		rlSetLineWidth(5);
		if(!hideWorld) for(int i=0; i<mdl.meshCount; i++) {
			
			bool visible = !usePVS || (pvsRes.visMeshCount > 0 && (pvsRes.visible[i / 32] & (1 << (i % 32))) != 0);
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
		
		rlDrawRenderBatchActive();
		
		EndMode3D();
		
		{
			Font fnt = GetFontDefault();
			int line = 0;
			
			DrawText(TextFormat("(P)VS on: %s", usePVS ? "on" : "off"), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("(L)ock PVS position: %s", lockPVS ? "on" : "off"), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("Visible meshes: %ld", pvsRes.visMeshCount), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			line++;
			
			DrawText(TextFormat("PVS calc (X) time: %.2f sec", raysComputeTime), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			DrawText(TextFormat("Total cubemaps: %ld (%d per sec)", totalRays, (int) (totalRays / raysComputeTime)), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
			line++;
			
			DrawText(TextFormat("(H)ide world: %s", hideWorld ? "on" : "off"), 0, fnt.baseSize * line * 2, fnt.baseSize * 2, BLACK); line++;
		}

        EndDrawing();
		t++;
	}
	
	pvsUnloadDB(&pvsDB);
	pvsUnloadModelData(&pvsMdlData);

	UnloadTexture(tex);
	UnloadModel(mdl);
	UnloadModel(cube);
	UnloadShader(shader);
	
    CloseWindow();
	
	return 0;
}
