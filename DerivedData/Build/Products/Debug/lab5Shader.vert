
attribute vec4 position;
attribute vec4 normal1;
attribute vec4 color1;
attribute vec2 textureCoord;
attribute vec4 tangent;

uniform mat4 local2clip;
uniform mat4 modelviewMatrix;
uniform mat3 normalMatrix;
uniform mat4 shadowMatrix;

varying vec4 pcolor;
varying vec3 N;
varying vec4 p;
varying vec2 fragTextureCoord;
varying vec4 shadowCoord;
varying vec3 lightDir;
varying vec3 viewDir;

uniform float Time; // The animation time
// Wave parameters
uniform float K;         // Wavenumber
uniform float Velocity; // Wave's velocity
uniform float Amp;         // Wave's amplitude

vec3 Gravity = vec3(0.0,-0.1,0.0);
float ParticleLifetime = 15.0;

uniform bool flagAnimation;
varying float flagAnimation1;
varying float Transp;

void main(){









    vec4 lightPosition = vec4(15.0, 15.0, 15.0, 1.0);

    N = normalize(normalMatrix * normal1.xyz);  //normal

    vec3 N1 = normalize(normalMatrix * vec3(0.0, 1.0, 0.0));
    vec3 T = normalize(normalMatrix * vec3(0.0, 0.0, 1.0));
    vec3 B = normalize(normalMatrix * vec3(1.0, 0.0, 0.0));

    mat3 toObjectLocal = mat3(T.x, B.x, N1.x,
                             T.y, B.y, N1.y,
                             T.z, B.z, N1.z);
    
    p = modelviewMatrix * position; //vertex position
    
    lightDir = normalize(toObjectLocal * (lightPosition - p).xyz);
    viewDir = normalize(toObjectLocal * normalize(-p.xyz));
    
    pcolor = color1;    
    gl_Position = local2clip * position;
    fragTextureCoord = textureCoord;
    
    shadowCoord = shadowMatrix * position;
    

    
    //Animation Particles
    if (flagAnimation == true) {
        vec3 pos = vec3(0.0);
        Transp = 0.0;

        // Particle doesn't exist until the start time
        if( Time >tangent.w ) {
            float t = Time - tangent.w;
        
            if( t <ParticleLifetime ) {
                pos = tangent.xyz * t + Gravity * t * t;
                Transp = 1.0 - t / ParticleLifetime;
            }
        }
        gl_Position = local2clip * vec4(pos, 1.0);
        flagAnimation1 = 1.0;
    }
    else {
        flagAnimation1 = 0.0;
    }
}


