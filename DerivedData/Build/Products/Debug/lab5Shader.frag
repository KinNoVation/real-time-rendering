uniform mat3 inverseModelview;

varying vec4 pcolor;
varying vec3 N;
varying vec4 p;
varying vec2 fragTextureCoord;
varying vec3 lightDir;
varying vec3 viewDir;

varying vec4 shadowCoord;

uniform sampler2D tex1;
uniform bool flagTexture;
uniform bool flagCubmap;
uniform bool flagShadowMap;
uniform samplerCube cubeMap;
uniform bool flagShadowMatrix;

uniform sampler2DShadow shadowMap;
//uniform sampler2D shadowMap;

uniform sampler2D normalMap;
uniform bool flagNormalMap;

uniform bool flagCartoon;
float EdgeThreshold = 0.0001;
uniform bool flagEdge;

int Width = 1024;
int Height = 1024;

varying float flagAnimation1;
varying float Transp;

float luma( vec3 color ) {
    return 0.2126 * color.r + 0.7152 * color.g +
            0.0722 * color.b;
}

vec4 pass2()
{
    float dx = 1.0 / float(Width);
    float dy = 1.0 / float(Height);
    
    float s00 = luma(texture2D( tex1,
                     fragTextureCoord + vec2(-dx,dy) ).rgb);
    float s10 = luma(texture2D( tex1,
                     fragTextureCoord + vec2(-dx,0.0) ).rgb);
    float s20 = luma(texture2D( tex1,
                     fragTextureCoord + vec2(-dx,-dy) ).rgb);
    float s01 = luma(texture2D( tex1,
                     fragTextureCoord + vec2(0.0,dy) ).rgb);
    float s21 = luma(texture2D( tex1,
                     fragTextureCoord + vec2(0.0,-dy) ).rgb);
    float s02 = luma(texture2D( tex1,
                     fragTextureCoord + vec2(dx, dy) ).rgb);
    float s12 = luma(texture2D( tex1,
                     fragTextureCoord + vec2(dx, 0.0) ).rgb);
    float s22 = luma(texture2D( tex1,
                     fragTextureCoord + vec2(dx, -dy) ).rgb);

    float sx = s00 + 2.0 * s10 + s20 - (s02 + 2.0 * s12 + s22);
    float sy = s00 + 2.0 * s01 + s02 - (s20 + 2.0 * s21 + s22);

    float dist = sx * sx + sy * sy;

    if( dist>EdgeThreshold )
    return vec4(1.0);
    else
        return vec4(0.0,0.0,0.0,1.0);
}

void main() {
    
    vec3 texColor = vec3(0.0);
    vec4 cubeMapTex = vec4(0.0);
    vec3 totalColor = vec3(0.0);
    
    vec3 N1 = normalize(texture2D(normalMap, fragTextureCoord).rgb);
    

    
    vec4 lightPosition = vec4(15.0, 15.0, 15.0, 1.0);
    
    vec3 L = normalize((lightPosition - p).xyz); //light vector, points to light source
    vec3 Camera = normalize((vec4(0.0) - p).xyz);
    vec3 V = normalize(-p.xyz);
    //    vec3 R = reflect(-L, N);
    vec3 R = reflect(-Camera, N);
    vec3 R1 = inverseModelview * R;
    vec3 ambient = vec3(0.0);
    float LDotN = abs(dot(L, N));   //here we should use max(dot(L, N), 0.0);

    vec3 diffuse = vec3(pcolor) * LDotN;
    vec3 specular = vec3(0.0);
    
    if (flagCartoon == true) {
        float cartoon = floor(LDotN * 4.0) / 4.0;
        diffuse = vec3(pcolor) * cartoon;
    }

    if (LDotN > 0.0) {
        specular = vec3(1.0) * pow(max(dot(R,V), 0.0), 20.0);
    }

    
    if(fragTextureCoord.x > 0.0 && fragTextureCoord.y > 0.0 && (flagTexture == true)) {
        if (flagEdge == true) {
            texColor = pass2().rgb;
        }
        else {
            texColor = texture2D(tex1, fragTextureCoord).rgb;
        }
//        texColor = texture2D(normalMap, fragTextureCoord).rgb;
//        texColor = texture2D(shadowMap, fragTextureCoord).rgb;
        totalColor = ambient + texColor + specular;
    }
    else {
        totalColor = ambient + diffuse + specular;
    }
    if(flagCubmap == true) {
        cubeMapTex = textureCube(cubeMap, R1);
        totalColor = ambient + cubeMapTex.rgb + specular;
    }
    if (flagShadowMap == true) {
        float depth = (shadow2DProj(shadowMap, shadowCoord).x == 1.0 ? 1.0 : 0.3);
        totalColor = totalColor * depth;
    }
    
    if (flagNormalMap == true) {
        vec3 r = reflect( -lightDir, N1 );
        float sDotN = abs( dot(lightDir, N1));
        diffuse = vec3(pcolor) * sDotN;
        if( sDotN > 0.0 ) {
            specular = vec3(1.0) * pow(max( dot(r,viewDir), 0.0 ), 20.0);
        }
        totalColor = ambient + diffuse + specular;
    }
    if (flagEdge == true) {
        totalColor = pass2().rgb;
    }

    
    gl_FragColor = vec4(totalColor, 1.0);
	
    if (flagAnimation1 == 1.0) {
        gl_FragColor = pcolor;
        gl_FragColor.a *= Transp;
    }
}
