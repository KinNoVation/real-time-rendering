//
//  main.h
//  lab4
//
//  Created by Yuxiang Kou on 10/29/12.
//  Copyright (c) 2012 Yuxiang Kou. All rights reserved.
//

#ifndef lab4_main_h
#define lab4_main_h
#include <GL/glew.h>
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lodepng.h"
#include <iostream>
#include "glm.hpp"
#include "matrix_transform.hpp"
#include "func_geometric.hpp"
#include "random.hpp"
#include "ply.h"
#define VBONUMBER 5     //here we have cube, cone, cylinder, sphere, bunny
#define WALLSCALAR 50.0
#define WALLTHICK 0.01

using namespace std;

unsigned char imageData[1048576];    //512x512x4, for the desktop surface
int imageWidth;
int imageHeight;
unsigned char imageData1[65536];    //128x128x4, px
int imageWidth1;
int imageHeight1;
unsigned char imageData2[65536];    //128x128x4, py
int imageWidth2;
int imageHeight2;
unsigned char imageData3[65536];    //128x128x4, pz
int imageWidth3;
int imageHeight3;
unsigned char imageData4[65536];    //128x128x4, nx
int imageWidth4;
int imageHeight4;
unsigned char imageData5[65536];    //128x128x4, ny
int imageWidth5;
int imageHeight5;
unsigned char imageData6[65536];    //128x128x4, nz
int imageWidth6;
int imageHeight6;
unsigned char imageDataNormalMap[4194304];  //256x256x4
int imageWidthNormalMap;
int imageHeightNormalMap;
//unsigned char* cubeMap[6];
//int cubeMapWidth[6];
//int cubeMapHeight[6];
GLuint textureId[7];
GLuint cubemap_texture;
GLuint normalMapHandle;

float cubeSize = 50;
float cubeColor[4] = {0.8588, 0.33333, 0.2549, 1};
float coneBaseRadius = 50;
float coneHeight = 50;
int coneSlices = 20;
int coneStacks = 20;
float coneColor[4] = {0.8588, 0.2549, 0.67843, 1};
float cylinderBaseRadius = 50;
float cylinderTopRadius = 50;
float cylinderHeight = 50;
int cylinderSlices = 20;
int cylinderStacks = 20;
float cylinderColor[4] = {0.61176, 0.8588, 0.2549, 1};
float cylinderColor1[4] = {219.0/255.0, 191.0/255.0, 65.0/255.0, 1};
float sphereRadius = 50;
int sphereSlices = 50;
int sphereStacks = 50;
float sphereColor[4] = {0.2549, 0.50588, 0.8588, 1};
float sphereColor1[4] = {63.0/255.0, 177.0/255.0, 152.0/255.0, 1};
float sphereColor2[4] = {250.0/255.0, 1.0, 196.0/255.0, 1};
float bunnyColor[4] = {255.0/255.0, 76.0/255.0, 185.0/255.0, 1};

glm::mat4 transformMatrix = glm::mat4(1.0f);

typedef struct {
    float location[4];      // (x, y, z, w)
    float normal[4];        // (x, y, z, w)
    float color[4];         // (r, g, b, w)
    float color1[4];
    float color2[4];
    float textureCoord[2];
    float tangent[4];       // also store particle parameters
} Vertex;
//NOTE: here we store different colors in each vertex, and we can pick colors via the offset value in glVertexAttribPointer() function.
//16: normal vector
//32: color
//48: color1
//64: color2

typedef struct Face {
    unsigned int count;
    unsigned int *vertices;
    float nx, ny, nz;
} Face;

char *string_list[] = {
    (char*)"x",
    (char*)"y",
    (char*)"z",
    (char*)"nx",
    (char*)"ny",
    (char*)"nz",
    (char*)"vertex_indices"
};

float cx, cy, cz;
float x_min, x_max, y_min, y_max, z_min, z_max;

Vertex** vertices = 0;
Face** faces = 0;
unsigned int vertexcount;
unsigned int facecount;
int vertexnormals = 0;
int facenormals = 0;

Vertex cubeVertices[8];
GLuint cubeIndices[26];
int cubeIndexTotal = 20;
//ATTENTION: here we need to make sure the size of cubeIndexTotal is correct, if we make it too large, the unused elements' values will be recognized as indices, then extra (useless) triangles (lines) will appear.

int coneIndexTotal;
int coneVertexTotal;
Vertex *coneVertices;
GLuint *coneIndices;

int sphereVertexTotal;
int sphereIndexTotal;
Vertex *sphereVertices;
GLuint *sphereIndices;

int nindices;       // the total number of indices of cylinder
Vertex *cyverts;   // cylinder vertices
GLuint *cindices;   // array that store the indices

int bunnyVertexTotal;
int bunnyIndexTotal;
Vertex *bunnyVertices;
GLuint *bunnyIndices;

GLuint vboHandle[VBONUMBER];   // a VBO that contains interleaved positions and colors
GLuint indexVBO[VBONUMBER];

// glut mouse control
int xform_mode = 0;
#define XFORM_NONE    0
#define XFORM_ROTATE  1
#define XFORM_SCALE 2

int press_x, press_y;
int release_x, release_y;
float x_angle = 0.0;
float y_angle = 0.0;
float scale_size = 1;
bool flagScene = true;
bool flagWireMode = false;

float liftLeftArmAngle = 0;
float liftRightArmAngle = 0;
int robotStep = 0;
//0: stand
//1: right step
//2: left step

GLuint programObject;
GLuint SetupGLSL(char*);

glm::mat4 projectionMatrix = glm::perspective(60.0f, 1.0f, 10.0f, 100.0f);
//glm::mat4 viewMatrix = glm::lookAt(glm::vec3(20.0, 10.0, 10.0),
//                                   glm::vec3(0.0, 5.0, 0.0),
//                                   glm::vec3(0.0, 1.0, 0.0));
glm::mat4 viewMatrixCamera = glm::lookAt(glm::vec3(20.0, 10.0, 10.0),
                                         glm::vec3(0.0, 5.0, 0.0),
                                         glm::vec3(0.0, 1.0, 0.0));
glm::mat4 viewMatrixLight = glm::lookAt(glm::vec3(15.0, 15.0, 15.0),
                                        glm::vec3(0.0, 0.0, 0.0),
                                        glm::vec3(0.0, 1.0, 0.0));
glm::mat4 viewMatrix = viewMatrixCamera;
glm::mat4 biasMatrix = glm::mat4(0.5, 0.0, 0.0, 0.0,
                                 0.0, 0.5, 0.0, 0.0,
                                 0.0, 0.0, 0.5, 0.0,
                                 0.5, 0.5, 0.5, 1.0);
//ATTENTION: the biasMatrix is mat4(0.5, 0.0, 0.0, 0.5,
//                                  0.0, 0.5, 0.0, 0.5,
//                                  0.0, 0.0, 0.5, 0.5,
//                                  0.0, 0.0, 0.0, 1.0);
//However, We must TRANSPOSE it!
glm::mat4 shadowMatrix = biasMatrix * projectionMatrix * viewMatrixLight;
//vec4 lightPosition = vec4(30.0, 80.0, 50.0, 1.0);
GLuint positionShader;
GLuint colorShader;
GLuint normalShader;
GLuint textureCoordShader;
GLuint tangentShader;

GLuint local2clipShader;
GLuint modelviewMatrixShader;
GLuint normalMatrixShader;
GLint textureShader;
GLint flagTextureShader;
GLint flagCubemapShader;
GLint cubeMapTextureShader;
GLuint inverseModelviewShader;

GLuint shadowMatrixShader;
GLuint depthTex;
GLuint shadowFBO;
GLint flagShadowMapShader;
GLint flagShadowMatrixShader;
GLuint shadowMapShader;

GLuint flagNormalMapShader;
GLuint normalMapShader;
GLuint flagCartoonShader;

bool flagEdgeDetection = false;
bool cartoon = false;
float ntime = 0.0;
GLuint timeShader;
GLuint flagEdgeShader;
glm::vec3 deskMove = glm::vec3(0.0, 0.0, -7.0);

/**** trackball ****/
int winWidth = 600;
int winHeight = 600;

float angle = 0.0;
//axis[3], trans[3];

bool recordView = false;

glm::vec3 lastPos = glm::vec3(0.0);
glm::vec3 axis = glm::vec3(0.0);
//float lastPos[3] = {0.0, 0.0, 0.0};
int curx, cury;
int startX, startY;

glm::mat4 trackball = glm::mat4(1.0f);

float animationTime = 0.0;
float timeInterval = 150;
GLuint initVelocityShader;
GLuint startTimeShader;
GLuint flagAnimationShader;
bool flagAnimation = false;
float particleStartTime = 0.0;
float particleRate = 0.00075f;

#endif
