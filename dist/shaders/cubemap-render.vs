#version 430
// Input uniform values
uniform int viewId;

// Input vertex attributes
in vec3 vertexPosition;

layout(std430, binding = 2) readonly restrict buffer commonData {
    mat4 matVP[6];
    mat4 matVPInv[6];
	
	uvec4 gridSize;
	vec4 aabbMin;
	vec4 aabbMax;
	uint dataPerCell;
};

void main() {
    // Calculate final vertex position
	vec4 pos = vec4(vertexPosition, 1.0) * matVP[viewId];
	
	gl_ClipDistance[0] = dot(pos, vec4(0.0, 1.0, 0.0, 1.0));
	gl_ClipDistance[1] = dot(pos, vec4(0.0, -1.0, 0.0, 1.0));
	gl_ClipDistance[2] = dot(pos, vec4(1.0, 0.0, 0.0, 1.0));
	gl_ClipDistance[3] = dot(pos, vec4(-1.0, 0.0, 0.0, 1.0));
	
	vec2 offset = vec2(float(viewId % 3) / 1.5 - 2.0 / 3.0, float(viewId / 3) - 0.5);
	offset *= pos.w;
	
    gl_Position = pos * vec4(1.0 / 3.0, 0.5, 1.0, 1.0) + vec4(offset.x, offset.y, 0.0, 0.0);
}
