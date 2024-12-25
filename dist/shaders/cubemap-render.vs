#version 430

// Input vertex attributes
in vec3 vertexPosition;

struct Viewport {
	mat4 mat;
	mat4 invMat;
};

layout(std430, binding = 2) readonly restrict buffer commonData {
	Viewport views[64 * 6];
	
    uvec4 gridSize;
	vec4 aabbMin;
	vec4 aabbMax;
	uvec2 dataPerCell; //datapercell, cubemapsinimage
};

layout(std430, binding = 3) readonly restrict buffer instanceData {
    uint instanceViewId[];
};

void main() {
	uint viewId = instanceViewId[gl_InstanceID];
	
    // Calculate final vertex position
	vec4 pos = vec4(vertexPosition, 1.0) * views[viewId].mat;
	
	gl_ClipDistance[0] = dot(pos, vec4(0.0, 1.0, 0.0, 1.0));
	gl_ClipDistance[1] = dot(pos, vec4(0.0, -1.0, 0.0, 1.0));
	gl_ClipDistance[2] = dot(pos, vec4(1.0, 0.0, 0.0, 1.0));
	gl_ClipDistance[3] = dot(pos, vec4(-1.0, 0.0, 0.0, 1.0));
	
	vec2 offset = vec2((float(viewId % 6) * 2.0 + 1.0) / 6.0 - 1.0, (float(viewId / 6) * 2.0 + 1.0) / dataPerCell.y - 1.0);
	offset *= pos.w;
	
    gl_Position = pos * vec4(1.0 / 6.0, 1.0 / dataPerCell.y, 1.0, 1.0) + vec4(offset.x, offset.y, 0.0, 0.0);
}
