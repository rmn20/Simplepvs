#version 430

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;
in vec4 vertexColor;

// Input uniform values
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
uniform mat4 matView;

// Output vertex attributes (to fragment shader)
out vec3 fragPosition;
out vec2 fragTexCoord;
out vec4 fragColor;
out vec3 fragNormal;

// NOTE: Add here your custom variables

void main() {
	const float pi = 3.14159265359;
    // Send vertex attributes to fragment shader
    fragPosition = vec3(matModel*vec4(vertexPosition, 1.0));
    fragTexCoord = vertexTexCoord;
    fragColor = vertexColor;
    fragNormal = normalize(vec3(matNormal * vec4(vertexNormal, 1.0)));

    // Calculate final vertex position
	vec4 pos = mvp * vec4(vertexPosition, 1.0);
    gl_Position = pos;// * vec4(0.5, 0.5, 1.0, 1.0);
	
	/*gl_ClipDistance[0] = dot(pos, vec4(0.0, 1.0, 0.0, 1.0));
	gl_ClipDistance[1] = dot(pos, vec4(0.0, -1.0, 0.0, 1.0));
	gl_ClipDistance[2] = dot(pos, vec4(1.0, 0.0, 0.0, 1.0));
	gl_ClipDistance[3] = dot(pos, vec4(-1.0, 0.0, 0.0, 1.0));*/
}
