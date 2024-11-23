#version 430

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

layout(std430, binding = 1) restrict buffer gridData {
    uint cellGrid[];
};

layout(std430, binding = 2) restrict buffer commonData {
    mat4 matVP[6];
    mat4 matVPInv[6];
	
	uvec4 gridSize;
	vec4 aabbMin;
	vec4 aabbMax;
	/*uvec4 outAABB[6];
	uvec4 camPos;*/
	uint dataPerCell;
};

uniform sampler2D cubeIdxTex;
uniform sampler2D cubeDepthTex;

void main() {
	uint viewId = gl_GlobalInvocationID.z;
	int texSize = textureSize(cubeIdxTex, 0).y / 2;
	
	ivec2 texelPos = ivec2(gl_GlobalInvocationID.xy) + ivec2((viewId % 3) * texSize, viewId / 3 * texSize);
	
	uint meshId = uint(texelFetch(cubeIdxTex, texelPos, 0).r);
	if(meshId == 0) return;
	meshId--;
	
	//Calculate ray data
	mat4 invMat = matVPInv[viewId];
	
	vec4 viewPos = vec4(float(gl_GlobalInvocationID.x), float(gl_GlobalInvocationID.y), 0.0, 1.0);
	viewPos.xy = viewPos.xy / float(texSize) * 2.0 - 1.0;
	
	vec4 rayStart = viewPos * invMat;
	rayStart.xyz /= rayStart.w;
	
	viewPos.z = texelFetch(cubeDepthTex, texelPos, 0).r * 2.0f - 1.0f;
	vec4 rayEnd = viewPos * invMat;
	rayEnd.xyz /= rayEnd.w;
	
	vec3 rayDir = (rayEnd - rayStart).xyz;
	float rayLen = length(rayDir);
	float rayScaledLen = length(rayDir * vec3(gridSize.xyz) / (aabbMax - aabbMin).xyz);
	rayDir /= rayLen;
	
	//Ray trace using branchless DDA (https://www.shadertoy.com/view/4dX3zl)
	vec3 mapRayPos = (rayStart - aabbMin).xyz * vec3(gridSize.xyz) / (aabbMax - aabbMin).xyz;
	ivec3 mapPos = ivec3(floor(mapRayPos - vec3(lessThan(ceil(mapRayPos), vec3(0.0))))); //hopefully proper rounding
	
	vec3 mapRayDir = normalize(rayDir * vec3(gridSize.xyz) / (aabbMax - aabbMin).xyz);
	
	ivec3 stepIsPositive = ivec3(greaterThan(mapRayDir, vec3(0.0)));
	ivec3 step = ivec3(greaterThan(mapRayDir, vec3(0.0))) * 2 - 1;
	
	vec3 tDelta = vec3(1.0) / abs(mapRayDir);
	
	vec3 tMax = ((mapPos - mapRayPos) * step + step * 0.5 + 0.5) * tDelta;
	
	#define FLT_MAX 3.402823466e+38
	if(mapRayDir.x == 0.0) tMax.x = FLT_MAX;
	if(mapRayDir.y == 0.0) tMax.y = FLT_MAX;
	if(mapRayDir.z == 0.0) tMax.z = FLT_MAX;
	
	float t = 0;
	while(t <= rayScaledLen) {
		//Update cell visibility
		if(greaterThanEqual(mapPos, ivec3(0)) == bvec3(true) && lessThan(mapPos, ivec3(gridSize.xyz)) == bvec3(true)) {
			uint cellPtr = (mapPos.x + mapPos.y * gridSize.x + mapPos.z * gridSize.x * gridSize.y) * dataPerCell;
			atomicOr(cellGrid[cellPtr + (meshId >> 5)], 1 << (meshId & 31));
		}
		
		if(tMax.x < tMax.y) {
			if(tMax.x < tMax.z) {
				if(step.x == -1 && mapPos.x <= 0) break;
				else if(step.x == 1 && mapPos.x >= (gridSize.x - 1)) break;
				
				mapPos.x += step.x;
				t = tMax.x;
				tMax.x += tDelta.x;
			} else {
				if(step.z == -1 && mapPos.z <= 0) break;
				else if(step.z == 1 && mapPos.z >= (gridSize.z - 1)) break;
				
				mapPos.z += step.z;
				t = tMax.z;
				tMax.z += tDelta.z;
			}
		} else {
			if(tMax.y < tMax.z) {
				if(step.y == -1 && mapPos.y <= 0) break;
				else if(step.y == 1 && mapPos.y >= (gridSize.y - 1)) break;
				
				mapPos.y += step.y;
				t = tMax.y;
				tMax.y += tDelta.y;
			} else {
				if(step.z == -1 && mapPos.z <= 0) break;
				else if(step.z == 1 && mapPos.z >= (gridSize.z - 1)) break;
				
				mapPos.z += step.z;
				t = tMax.z;
				tMax.z += tDelta.z;
			}
		}
	}
}