#version 430
// Input uniform values
uniform int meshId;

// Output fragment color
out vec4 outputId;

void main() {
	outputId = vec4(float(gl_FrontFacing ? meshId : 0), 0.0, 0.0, 1.0);
}