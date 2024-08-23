#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;
in vec3 viewDir;

// Input uniform values
uniform sampler2D texture0;
uniform vec4 colDiffuse;

// Output fragment color
out vec4 finalColor;

uniform vec3 viewPos;

void main()
{
    // Texel color fetching from texture sampler
    vec4 texelColor = texture(texture0, fragTexCoord);
	
	vec4 col = texelColor * colDiffuse * fragColor;
	
	/*if(length(fragNormal) > 0.000001) */col.rgb *= pow(abs(dot(fragNormal, normalize(viewPos - fragPosition))), 0.5) * 0.5 + 0.5;

	finalColor = col;
    //finalColor = pow(finalColor, vec4(1.0/2.2));
}
