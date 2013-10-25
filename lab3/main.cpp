
///////////////////////////////////////////////////////////////////////////////////////////////
//
//  CSE 5542 lab 4 - Yuxiang Kou
//  kou.29@osu.edu
//  614-906-6096
//
///////////////////////////////////////////////////////////////////////////////////////////////
#include "main.h"

void trackball_ptov(int x, int y, int width, int height, float v[3]) {
	float d, a;
	/* project x,y onto a hemisphere centered within width, height , note z is up here*/
    v[0] = x / width * 2.0 - 1.0;
    v[1] = y / height * 2.0 - 1.0;
//    if (d > 1.0) {
//        d = 1.0;
//    }
    v[2] = sqrt(1.0 - v[0] * v[0] - v[1] * v[1]);
	
    
//    v[0] = (2.0*x - width) / width;
//	v[1] = (height - 2.0F*y) / height;
//    d = sqrt(v[0]*v[0] + v[1]*v[1]);
//	v[2] = cos((M_PI/2.0) * ((d < 1.0) ? d : 1.0));
	a = 1.0 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	v[0] *= a;
	v[1] *= a;
	v[2] *= a;
}

void initCube() {
    //cube bottom center is local coordinates origin
    for (int i = 0; i < 8; i ++) {
        //first bottom 4 points count clockwize, then top 4 points count clockwize
        cubeVertices[i].location[0] = round(cosf(M_PI_2 * i + M_PI_4));
        cubeVertices[i].location[1] = round(sinf(M_PI_2 * i + M_PI_4));
        if (i < 4) {
            cubeVertices[i].location[2] = -1;
        }
        else {
            cubeVertices[i].location[2] = 1;
        }
        cubeVertices[i].location[3] = 1;
        
        cubeVertices[i].normal[0] = cubeVertices[i].location[0];
        cubeVertices[i].normal[1] = cubeVertices[i].location[1];
        cubeVertices[i].normal[2] = cubeVertices[i].location[2];
        cubeVertices[i].normal[3] = 0;
        
        cubeVertices[i].color[0] = cubeColor[0];
        cubeVertices[i].color[1] = cubeColor[1];
        cubeVertices[i].color[2] = cubeColor[2];
        cubeVertices[i].color[3] = cubeColor[3];
        
        cubeVertices[i].color1[0] = cylinderColor[0];
        cubeVertices[i].color1[1] = cylinderColor[1];
        cubeVertices[i].color1[2] = cylinderColor[2];
        cubeVertices[i].color1[3] = cylinderColor[3];
        
        cout << "cube " << i << " xyzw: (" << cubeVertices[i].location[0] << ", " << cubeVertices[i].location[1] << ", " << cubeVertices[i].location[2] << ", " << cubeVertices[i].location[3] << ") " << endl;
    }
    
    cubeVertices[0].tangent[0] = 0.5;
    cubeVertices[0].tangent[1] = 0.5;
    cubeVertices[0].tangent[2] = 1.0;
    cubeVertices[0].tangent[3] = 1.0;
    
    cubeVertices[1].tangent[0] = -0.5;
    cubeVertices[1].tangent[1] = 0.5;
    cubeVertices[1].tangent[2] = 1.0;
    cubeVertices[1].tangent[3] = 1.0;
    
    cubeVertices[2].tangent[0] = -0.5;
    cubeVertices[2].tangent[1] = -0.5;
    cubeVertices[2].tangent[2] = 1.0;
    cubeVertices[2].tangent[3] = 1.0;
    
    cubeVertices[3].tangent[0] = 0.5;
    cubeVertices[3].tangent[1] = -0.5;
    cubeVertices[3].tangent[2] = 1.0;
    cubeVertices[3].tangent[3] = 1.0;

    cubeVertices[4].tangent[0] = -0.5;
    cubeVertices[4].tangent[1] = -0.5;
    cubeVertices[4].tangent[2] = 1.0;
    cubeVertices[4].tangent[3] = 1.0;
    
    cubeVertices[5].tangent[0] = 0.5;
    cubeVertices[5].tangent[1] = -0.5;
    cubeVertices[5].tangent[2] = 1.0;
    cubeVertices[5].tangent[3] = 1.0;

    cubeVertices[6].tangent[0] = 0.5;
    cubeVertices[6].tangent[1] = 0.5;
    cubeVertices[6].tangent[2] = 1.0;
    cubeVertices[6].tangent[3] = 1.0;

    cubeVertices[7].tangent[0] = -0.5;
    cubeVertices[7].tangent[1] = 0.5;
    cubeVertices[7].tangent[2] = 1.0;
    cubeVertices[7].tangent[3] = 1.0;
    
    cubeVertices[0].textureCoord[0] = 0.0;
    cubeVertices[0].textureCoord[1] = 1.0;
    cubeVertices[1].textureCoord[0] = 0.0;
    cubeVertices[1].textureCoord[1] = 0.0;
    cubeVertices[2].textureCoord[0] = -10000.0;
    cubeVertices[2].textureCoord[1] = -10000.0;
    cubeVertices[3].textureCoord[0] = -10000.0;
    cubeVertices[3].textureCoord[1] = -10000.0;
    
    cubeVertices[4].textureCoord[0] = 1.0;
    cubeVertices[4].textureCoord[1] = 1.0;
    cubeVertices[5].textureCoord[0] = 1.0;
    cubeVertices[5].textureCoord[1] = 0.0;
    cubeVertices[6].textureCoord[0] = -10000.0;
    cubeVertices[6].textureCoord[1] = -10000.0;
    cubeVertices[7].textureCoord[0] = -10000.0;
    cubeVertices[7].textureCoord[1] = -10000.0;
    
    int n = 0;
    
    //these are the bottom
    cubeIndices[n++] = 1;
    cubeIndices[n++] = 0;
    cubeIndices[n++] = 2;
    cubeIndices[n++] = 3;
    cubeIndices[n++] = 3;
    cubeIndices[n++] = 0;
    
    //these are the side
    for (int j = 0; j < 1; j++) {
        for (int i = 0; i <= 4; i++) {
            int mi = i % 4;
            int idx = j * 4 + mi; // mesh[j][mi]
            int idx2 = (j + 1) * 4 + mi;
            cubeIndices[n++] = idx;
            cubeIndices[n++] = idx2;
        }
    }
    
    //these are the top
    cubeIndices[n++] = 4;
    cubeIndices[n++] = 5;
    cubeIndices[n++] = 7;
    cubeIndices[n++] = 6;
    
    //    for (int i = 0; i < cubeIndexTotal ; i ++) {
    //        cout << "cubeIndices[" << i << "]: " << cubeIndices[i] << endl;
    //    }
    
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * 8, cubeVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * cubeIndexTotal, cubeIndices, GL_STATIC_DRAW);  // load the index data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up
}
void initCone() {
    //cone bottom center is local coordinates origin
    coneVertexTotal = coneSlices * coneStacks;
    coneVertices = new Vertex[coneVertexTotal];
    float Dangle = 2 * M_PI / (float)(coneSlices);
    for (int j = 0; j < coneStacks; j ++) {
        for (int i = 0; i < coneSlices; i ++) {
            int idx = j * coneSlices + i; // mesh[j][i]
            //            cout << "idx: " << idx << endl;
            float angle = Dangle * i;
            coneVertices[idx].location[0] = cos(angle) * (coneStacks - 1 - j) / (coneStacks - 1);
            coneVertices[idx].location[1] = sin(angle) * (coneStacks - 1 - j) / (coneStacks - 1);
            coneVertices[idx].location[2] = ((float)j) / (float)(coneStacks - 1);
            coneVertices[idx].location[3] = 1.0;
            
            coneVertices[idx].normal[0] = coneVertices[idx].location[0];
            coneVertices[idx].normal[1] = coneVertices[idx].location[1];
            coneVertices[idx].normal[2] = 0.5;  //I calculated this value on papaer, the cone radious is 1, height is 1.
            //HOWEVER, we need to make sure the xyz components are compatible.
            coneVertices[idx].normal[3] = 0;
            
            coneVertices[idx].color[0] = coneColor[0];
            coneVertices[idx].color[1] = coneColor[1];
            coneVertices[idx].color[2] = coneColor[2];
            coneVertices[idx].color[3] = 1.0;
            //            cout << "cone " << idx << " xyzw: (" << coneVertices[idx].location[0] << ", " << coneVertices[idx].location[1] << ", " << coneVertices[idx].location[2] << ", " << coneVertices[idx].location[3] << ") " << endl;
        }
    }
    
    coneIndexTotal = coneStacks * 2 * coneSlices;   // index for the side
    coneIndexTotal += (coneSlices - 1) * 3 ;   // index for the bottom
    coneIndexTotal += 2;    // for the transition between bottom and side
    coneIndices = new GLuint[coneIndexTotal];
    int n = 0;
    
    //this part is the cone bottom indices
    for (int i = 2; i < coneSlices; i++) {
        coneIndices[n++] = 0;
        coneIndices[n++] = i - 1;
        coneIndices[n++] = i;
    }
    
    //this part is for safe transition between bottom and side
    coneIndices[n++] = 0;
    coneIndices[n++] = 0;
    
    //this part is for side
    for (int j = 0; j < coneStacks - 1; j++) {
        for (int i = 0; i <= coneSlices; i++) {
            int mi = i % coneSlices;
            int idx = j * coneSlices + mi; // mesh[j][mi]
            int idx2 = (j + 1) * coneSlices + mi;
            coneIndices[n++] = idx;
            coneIndices[n++] = idx2;
        }
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[2]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * coneVertexTotal, coneVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[2]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * coneIndexTotal, coneIndices, GL_STATIC_DRAW);  // load the index data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up
}
void initSphere() {
    //sphere center is local coordinates origin
    sphereVertexTotal = sphereSlices * sphereStacks;
    sphereVertices = new Vertex[sphereVertexTotal];
    //    float Dangle = 2 * M_PI / (float)(sphereSlices);
    for (int j = 0; j < sphereStacks ; j ++) {
        for (int i = 0; i < sphereSlices; i ++) {
            int idx = j * sphereSlices + i; // mesh[j][i]
            float angleTheta = M_PI / (float)(sphereStacks - 1) * j;    //here -1 is necessary, or a hole will appear at the bottom
            float anglePhi =  2 * M_PI / (float)(sphereSlices) * i;
            
            sphereVertices[idx].location[0] = sin(angleTheta) * cos(anglePhi);
            sphereVertices[idx].location[1] = sin(angleTheta) * sin(anglePhi);
            sphereVertices[idx].location[2] = cos(angleTheta);
            sphereVertices[idx].location[3] = 1.0;
            
            sphereVertices[idx].normal[0] = sphereVertices[idx].location[0];
            sphereVertices[idx].normal[1] = sphereVertices[idx].location[1];
            sphereVertices[idx].normal[2] = sphereVertices[idx].location[2];
            sphereVertices[idx].normal[3] = 0.0;
            
            sphereVertices[idx].color[0] = sphereColor[0];
            sphereVertices[idx].color[1] = sphereColor[1];
            sphereVertices[idx].color[2] = sphereColor[2];
            sphereVertices[idx].color[3] = 1.0;
            
            sphereVertices[idx].color1[0] = sphereColor1[0];
            sphereVertices[idx].color1[1] = sphereColor1[1];
            sphereVertices[idx].color1[2] = sphereColor1[2];
            sphereVertices[idx].color1[3] = sphereColor1[3];
            
            sphereVertices[idx].color2[0] = sphereColor2[0];
            sphereVertices[idx].color2[1] = sphereColor2[1];
            sphereVertices[idx].color2[2] = sphereColor2[2];
            sphereVertices[idx].color2[3] = sphereColor2[3];
        }
    }
    
    sphereIndexTotal = sphereStacks * 2 * sphereSlices - 2;   // index for the side
    sphereIndices = new GLuint[sphereIndexTotal];
    int n = 0;
    //this part is for side
    for (int j = 0; j < sphereStacks - 1; j++) {
        for (int i = 0; i <= sphereSlices; i++) {
            int mi = i % sphereSlices;
            int idx = j * sphereSlices + mi; // mesh[j][mi]
            int idx2 = (j + 1) * sphereSlices + mi;
            sphereIndices[n++] = idx;
            sphereIndices[n++] = idx2;
        }
    }
    
    //    for (int i = 0; i < sphereIndexTotal ; i ++) {
    //        cout << "sphereIndices[" << i << "]: " << sphereIndices[i] << endl;
    //    }
    
    glm::vec3 v(0.0f);
    float velocity, theta, phi;
    
    //    GLfloat *data = new GLfloat[nParticles * 3];
    for( GLuint i = 0; i<sphereVertexTotal; i++ ) {
        // Pick the direction of the velocity
        theta = glm::mix(0.0f, (float)M_PI / 6.0f, (float)glm::linearRand(0.0, 1.0));
        phi = glm::mix(0.0f, (float)((float)M_PI * 2.0), (float)glm::linearRand(0.0, 1.0));
        
        v.x = sinf(theta) * cosf(phi);
        v.y = cosf(theta);
        v.z = sinf(theta) * sinf(phi);
        
        // Scale to set the magnitude of the velocity (speed)
        velocity = glm::mix(1.25f,1.5f, glm::linearRand(0.0, 1.0));
        v = v * velocity;
        
        sphereVertices[i].tangent[0] = v.x;
        sphereVertices[i].tangent[1] = v.y;
        sphereVertices[i].tangent[2] = v.z;
        sphereVertices[i].tangent[3] = particleStartTime;
        particleStartTime += particleRate;
        
        //        cout << "sphereVertices[" << i << "]: (" << sphereVertices[i].tangent[0] << ", " << sphereVertices[i].tangent[1] << ", " << sphereVertices[i].tangent[2] << ", " << sphereVertices[i].tangent[3] << ")" << endl;
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[3]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * sphereVertexTotal, sphereVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[3]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * sphereIndexTotal, sphereIndices, GL_STATIC_DRAW);  // load the index data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up
}
void initCylinder() {
    //cylinder bottom center is local coordinates origin
    int nvertices = cylinderSlices * cylinderStacks;
    cyverts = new Vertex[nvertices];
    
    printf(" M PI = %f\n", M_PI);
    float Dangle = 2*M_PI/(float)(cylinderSlices);
    
    for (int j = 0; j < cylinderStacks; j++)
        for (int i= 0; i < cylinderSlices; i++) {
            int idx = j * cylinderSlices + i; // mesh[j][i]
            float angle = Dangle * i;
            cyverts[idx].location[0] = cos(angle);
            cyverts[idx].location[1] = sin(angle);
            cyverts[idx].location[2] = j*1.0/(float)(cylinderStacks-1);
            cyverts[idx].location[3] = 1.0;
            
            cyverts[idx].normal[0] = cyverts[idx].location[0];
            cyverts[idx].normal[1] = cyverts[idx].location[1];
            cyverts[idx].normal[2] = 0.0;
            cyverts[idx].normal[3] = 0.0;
            
            cyverts[idx].color[0] = cylinderColor[0];
            cyverts[idx].color[1] = cylinderColor[1];
            cyverts[idx].color[2] = cylinderColor[2];
            cyverts[idx].color[3] = 1.0;
            
            cyverts[idx].color1[0] = cylinderColor1[0];
            cyverts[idx].color1[1] = cylinderColor1[1];
            cyverts[idx].color1[2] = cylinderColor1[2];
            cyverts[idx].color1[3] = 1.0;
        }
    // now create the index array
    
    nindices = cylinderStacks * 2 * cylinderSlices;
    nindices += cylinderSlices * 6 - 8;
    cindices = new GLuint[nindices];
    int n = 0;
    
    //for the bottom
    for (int i = 2; i <= cylinderSlices; i ++) {
        cindices[n++] = 0;
        cindices[n++] = i - 1;
        cindices[n++] = i;
    }
    cindices[n++] = 0;
    
    //for the side
    for (int j = 0; j < cylinderStacks - 1; j++) {
        for (int i = 0; i <= cylinderSlices; i++) {
            int mi = i % cylinderSlices;
            int idx = j * cylinderSlices + mi; // mesh[j][mi]
            int idx2 = (j + 1) * cylinderSlices + mi;
            cindices[n++] = idx;
            cindices[n++] = idx2;
        }
    }
    
    //for the top
    GLuint temp = cindices[n - 1];
    cindices[n++] = temp;
    cindices[n++] = temp;
    for (int i = 1; i < cylinderSlices - 1; i ++) {
        cindices[n++] = temp;
        cindices[n++] = temp + i;
        cindices[n++] = temp + i + 1;
    }
    
    //    for (int i = 0; i < nindices ; i ++) {
    //        cout << "cylinderIndices[" << i << "]: " << cindices[i] << endl;
    //    }
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);   // bind the first handle
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * nvertices, cyverts, GL_STATIC_DRAW); // allocate space and copy the position data over
    glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[0]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*nindices, cindices, GL_STATIC_DRAW);  // load the index data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up
}
void initBunny() {
    bunnyVertexTotal = vertexcount;
    bunnyIndexTotal = facecount * 3;
    printf("facecount : %d\n", facecount);
    bunnyVertices = new Vertex[bunnyVertexTotal];
    bunnyIndices = new GLuint[bunnyIndexTotal];
    
    for (int i = 0; i < bunnyVertexTotal; i ++) {
        bunnyVertices[i].location[0] = vertices[i]->location[0];
        bunnyVertices[i].location[1] = vertices[i]->location[1];
        bunnyVertices[i].location[2] = vertices[i]->location[2];
        bunnyVertices[i].location[3] = 1.0;
        
        bunnyVertices[i].normal[0] = vertices[i]->normal[0];
        bunnyVertices[i].normal[1] = vertices[i]->normal[1];
        bunnyVertices[i].normal[2] = vertices[i]->normal[2];
        bunnyVertices[i].normal[3] = 0.0;
        //        printf("bunny faces %d: (%f, %f, %f)\n", i, bunnyVertices[i].normal[0], bunnyVertices[i].normal[1], bunnyVertices[i].normal[2]);
        
        
        bunnyVertices[i].color[0] = bunnyColor[0];
        bunnyVertices[i].color[1] = bunnyColor[1];
        bunnyVertices[i].color[2] = bunnyColor[2];
        bunnyVertices[i].color[3] = bunnyColor[3];
        
    }
    
    for (int i = 0; i < bunnyIndexTotal; i += 3) {
        bunnyIndices[i] = faces[i/3]->vertices[0];
        bunnyIndices[i+1] = faces[i/3]->vertices[1];
        bunnyIndices[i+2] = faces[i/3]->vertices[2];
        //        printf("bunny faces %d: (%d, %d, %d)\n", i/3, bunnyIndices[i], bunnyIndices[i+1], bunnyIndices[i+2]);
    }
    
    float scaleBunny = 0.0;
    
    for (int i = 0; i < bunnyVertexTotal; i ++) {
        scaleBunny = glm::max(bunnyVertices[i].location[0], scaleBunny);
        scaleBunny = glm::max(bunnyVertices[i].location[1], scaleBunny);
        scaleBunny = glm::max(bunnyVertices[i].location[2], scaleBunny);
    }
    
    for (int i = 0; i < bunnyVertexTotal; i ++) {
        bunnyVertices[i].location[0] = bunnyVertices[i].location[0] / scaleBunny;
        bunnyVertices[i].location[1] = bunnyVertices[i].location[1] / scaleBunny;
        bunnyVertices[i].location[2] = bunnyVertices[i].location[2] / scaleBunny;
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[4]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * bunnyVertexTotal, bunnyVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[4]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * bunnyIndexTotal, bunnyIndices, GL_STATIC_DRAW);  // load the index data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up
}

void initTexture() {
    //set up texture mapping:
    //textureId[0] is for desktop image
    //textureId[1 - 7] is for walls
    glGenTextures(7, textureId); // generate one texture handle
//    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureId[0]); // bind this handle to a 2D texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidth, imageHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData); // usage same as described before
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    
//    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, textureId[1]); // bind this handle to a 2D texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidth1, imageHeight1, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData1); // usage same as described before
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    
    glBindTexture(GL_TEXTURE_2D, textureId[2]); // bind this handle to a 2D texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidth2, imageHeight2, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData2); // usage same as described before
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    
    glBindTexture(GL_TEXTURE_2D, textureId[3]); // bind this handle to a 2D texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidth3, imageHeight3, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData3); // usage same as described before

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    
    glBindTexture(GL_TEXTURE_2D, textureId[4]); // bind this handle to a 2D texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidth4, imageHeight4, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData4); // usage same as described before
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    
    glBindTexture(GL_TEXTURE_2D, textureId[5]); // bind this handle to a 2D texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidth5, imageHeight5, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData5); // usage same as described before

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    
    glBindTexture(GL_TEXTURE_2D, textureId[6]); // bind this handle to a 2D texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidth6, imageHeight6, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData6); // usage same as described before

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    
    glActiveTexture(GL_TEXTURE1);
    glGenTextures(1, &normalMapHandle);
    glBindTexture(GL_TEXTURE_2D, normalMapHandle); // bind this handle to a 2D texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidthNormalMap, imageHeightNormalMap, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageDataNormalMap); // usage same as described before
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glActiveTexture(GL_TEXTURE0);

    glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X, 0, GL_RGBA, imageWidth1, imageHeight1, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData1);
    glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Y, 0, GL_RGBA, imageWidth2, imageHeight2, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData2);
    glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Z, 0, GL_RGBA, imageWidth3, imageHeight3, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData3);
    glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 0, GL_RGBA, imageWidth4, imageHeight4, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData4);
    glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 0, GL_RGBA, imageWidth5, imageHeight5, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData5);
    glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, 0, GL_RGBA, imageWidth6, imageHeight6, 0, GL_RGBA, GL_UNSIGNED_BYTE, &imageData6);
    
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER,
                    GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);
}
void initVBO() {
    glGenBuffers(VBONUMBER, vboHandle);   // create an interleaved VBO object
    glGenBuffers(VBONUMBER, indexVBO);
    
    initCylinder();
    initCube();
    initCone();
    initSphere();
    initBunny();
}
void initShader() {
    positionShader = glGetAttribLocation(programObject, "position");
    colorShader = glGetAttribLocation(programObject, "color1");
    normalShader = glGetAttribLocation(programObject, "normal1");
    textureCoordShader = glGetAttribLocation(programObject, "textureCoord");
    
    local2clipShader = glGetUniformLocation(programObject, "local2clip"); // get the location of the uniform variable in the shader
    normalMatrixShader = glGetUniformLocation(programObject, "normalMatrix");
    modelviewMatrixShader = glGetUniformLocation(programObject, "modelviewMatrix");
    textureShader = glGetUniformLocation(programObject, "tex1");
    flagTextureShader = glGetUniformLocation(programObject, "flagTexture");
    cubeMapTextureShader = glGetUniformLocation(programObject, "cubeMap");
    inverseModelviewShader = glGetUniformLocation(programObject, "inverseModelview");
    flagCubemapShader = glGetUniformLocation(programObject, "flagCubmap");
    flagShadowMapShader = glGetUniformLocation(programObject, "flagShadowMap");
    flagShadowMatrixShader = glGetUniformLocation(programObject, "flagShadowMatrix");
    shadowMatrixShader = glGetUniformLocation(programObject, "shadowMatrix");
    shadowMapShader = glGetUniformLocation(programObject, "shadowMap");
    
    flagNormalMapShader = glGetUniformLocation(programObject, "flagNormalMap");
    normalMapShader = glGetUniformLocation(programObject, "normalMap");
    tangentShader = glGetAttribLocation(programObject, "tangent");
    
    flagCartoonShader = glGetUniformLocation(programObject, "flagCartoon");
    
    timeShader = glGetUniformLocation(programObject, "Time");
    flagEdgeShader = glGetUniformLocation(programObject, "flagEdge");
    flagAnimationShader = glGetUniformLocation(programObject, "flagAnimation");
    
    //    cout << "flagTextureShader: " << flagTextureShader << endl;
    
    glEnableVertexAttribArray(positionShader);
    glEnableVertexAttribArray(colorShader);
    glEnableVertexAttribArray(normalShader);
    glEnableVertexAttribArray(textureCoordShader);
    glEnableVertexAttribArray(tangentShader);
}
void initFBO() {
    glActiveTexture(GL_TEXTURE5);
    glGenFramebuffers(1, &shadowFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, shadowFBO);
    
    GLfloat border[] = {1.0f, 0.0f, 0.0f, 0.0f};
    
    glGenTextures(1, &depthTex);
    glActiveTexture(GL_TEXTURE5);
    glBindTexture(GL_TEXTURE_2D, depthTex);
    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, 1024, 1024, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, border);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LESS);
    
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);
    
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTex, 0);
    // check FBO status
    GLenum FBOstatus = glCheckFramebufferStatusEXT (GL_FRAMEBUFFER);
    if(FBOstatus != GL_FRAMEBUFFER_COMPLETE)
        printf("GL_FRAMEBUFFER_COMPLETE failed, CANNOT use FBO\n");
    glBindFramebuffer(GL_FRAMEBUFFER, 0); // go back to the default framebuffer
    glActiveTexture(GL_TEXTURE0);
}
void initParticle() {
    glm::vec3 v(0.0f);
    float velocity, theta, phi;

//    GLfloat *data = new GLfloat[nParticles * 3];
    for( GLuint i = 0; i<sphereVertexTotal; i++ ) {
        // Pick the direction of the velocity
        theta = glm::mix(0.0f, (float)M_PI / 6.0f, (float)glm::linearRand(0.0, 1.0));
        phi = glm::mix(0.0f, (float)((float)M_PI * 2.0), (float)glm::linearRand(0.0, 1.0));
        
        v.x = sinf(theta) * cosf(phi);
        v.y = cosf(theta);
        v.z = sinf(theta) * sinf(phi);
        
        // Scale to set the magnitude of the velocity (speed)
        velocity = glm::mix(1.25f,1.5f, glm::linearRand(0.0, 1.0));
        v = v * velocity;
        
        sphereVertices[i].tangent[0] = v.x;
        sphereVertices[i].tangent[1] = v.y;
        sphereVertices[i].tangent[2] = v.z;
        sphereVertices[i].tangent[3] = particleStartTime;
        particleStartTime += particleRate;
        
//        cout << "sphereVertices[" << i << "]: (" << sphereVertices[i].tangent[0] << ", " << sphereVertices[i].tangent[1] << ", " << sphereVertices[i].tangent[2] << ", " << sphereVertices[i].tangent[3] << ")" << endl;
    }
//    glBindBuffer(GL_ARRAY_BUFFER,initVel);
//    glBufferSubData(GL_ARRAY_BUFFER, 0,
//                    nParticles * 3 * sizeof(float), data);
    
//    float * data = new GLfloat[nParticles];
    
//    for( unsigned int i = 0; i<sphereVertexTotal; i++ ) {
//
//    }
//    glBindBuffer(GL_ARRAY_BUFFER,startTime);
//    glBufferSubData(GL_ARRAY_BUFFER, 0,
//                    nParticles * sizeof(float), data);
}

void displayCube() {
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    glVertexAttribPointer(tangentShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+88);
    
    //desk leg 1
    glm::mat4 leg1Matrix = glm::mat4(1.0f);
    
    leg1Matrix = glm::rotate(leg1Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg1Matrix = glm::rotate(leg1Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(scale_size, scale_size, scale_size));
    
    //    leg1Matrix = glm::translate(leg1Matrix, glm::vec3(-2.0, 2.0, 2.0));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    leg1Matrix = viewMatrix * leg1Matrix;
    glm::mat4 leg1ModelviewMatrix = leg1Matrix;
    glm::mat3 leg1NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg1Matrix)));
    leg1Matrix = projectionMatrix * leg1Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg1ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg1NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg1Matrix[0][0]);
    
    //    glm::mat3 inv = glm::inverse(glm::mat3(leg1ModelviewMatrix));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[1]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
}
void displayCylinder() {
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[0]);
    
    glEnableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
    //    glEnableClientState(GL_COLOR_ARRAY); // enable the color array on the client side
    
    // Tells OpenGL how to walk through the VBOs, i.e., how the data are packed
    // number of coordinates per vertex (4 here), type of the coordinates,
    // stride between consecutive vertices, and pointers to the first coordinate
    
    glVertexPointer(4,GL_FLOAT, sizeof(Vertex),(char*) NULL+0);
    glColorPointer(4,GL_FLOAT,  sizeof(Vertex),(char*) NULL+32);
    
    // Now we are ready to draw, using the triangle indices in the buffer
    // The following code uses GLM, the OpenGL mathematics library. This is to replace
    // the OpenGL2.x matrix functions which have been deprecated.
    // use the GLM helper functions to construct projection and modelview matrices
    
    //the base of lamp
    glPushMatrix();
    glRotatef(x_angle, 0, 1,0);
    glRotatef(y_angle, 1,0,0);
    glScalef(scale_size, scale_size, scale_size);
    glTranslatef(-1, 5.5, 1);
    glRotatef(90, 1, 0, 0);
    glScalef(1, 1, 0.5);
    glColor3f(1, 1, 0);
    glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0);
    glPopMatrix();
    
    //the middle of lamp
    glPushMatrix();
    glRotatef(x_angle, 0, 1,0);
    glRotatef(y_angle, 1,0,0);
    glScalef(scale_size, scale_size, scale_size);
    glTranslatef(-1, 7, 1);
    glRotatef(90, 1, 0, 0);
    glScalef(0.1, 0.1, 1.5);
    glColor3f(1, 0, 1);
    glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0);
    glPopMatrix();
    
    glDisableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
    glDisableClientState(GL_COLOR_ARRAY); // enable the color array on the client side
}
void displayCone() {
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[2]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[2]);
    
    glEnableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
    //    glEnableClientState(GL_COLOR_ARRAY); // enable the color array on the client side
    
    // Tells OpenGL how to walk through the VBOs, i.e., how the data are packed
    // number of coordinates per vertex (4 here), type of the coordinates,
    // stride between consecutive vertices, and pointers to the first coordinate
    
    glVertexPointer(4, GL_FLOAT, sizeof(Vertex), (char*) NULL + 0);
    glColorPointer(4, GL_FLOAT,  sizeof(Vertex), (char*) NULL + 32);
    
    // Now we are ready to draw, using the triangle indices in the buffer
    //
    // The following code uses GLM, the OpenGL mathematics library. This is to replace
    // the OpenGL2.x matrix functions which have been deprecated.
    //
    //
    // use the GLM helper functions to construct projection and modelview matrices
    
    //top of lamp
    glPushMatrix();
    glRotatef(x_angle, 0, 1,0);
    glRotatef(y_angle, 1,0,0);
    glScalef(scale_size, scale_size, scale_size);
    glTranslatef(0, 6.5, 0);
    glRotatef(-45, 0, 1, 0);
    glRotatef(-20, 1, 0, 0);
    glScalef(0.5, 0.5, 1.5);
    glColor3f(1, 1.5, 0.5);
    glDrawElements(GL_TRIANGLE_STRIP, coneIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    glPopMatrix();
    
    glDisableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
    glDisableClientState(GL_COLOR_ARRAY); // enable the color array on the client side
}
void displaySphere() {
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[3]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[3]);
    
    glEnableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
    //    glEnableClientState(GL_COLOR_ARRAY); // enable the color array on the client side
    
    // Tells OpenGL how to walk through the VBOs, i.e., how the data are packed
    // number of coordinates per vertex (4 here), type of the coordinates,
    // stride between consecutive vertices, and pointers to the first coordinate
    
    glVertexPointer(4, GL_FLOAT, sizeof(Vertex), (char*) NULL + 0);
    glColorPointer(4, GL_FLOAT, sizeof(Vertex), (char*) NULL + 32);
    
    // Now we are ready to draw, using the triangle indices in the buffer
    //
    // The following code uses GLM, the OpenGL mathematics library. This is to replace
    // the OpenGL2.x matrix functions which have been deprecated.
    //
    //
    // use the GLM helper functions to construct projection and modelview matrices
    
    //the sphere on the desk
    glPushMatrix();
    glRotatef(x_angle, 0, 1, 0);
    glRotatef(y_angle, 1, 0, 0);
    glScalef(scale_size, scale_size, scale_size);
    glTranslatef(1.2, 5.5, -1.2);
    glScalef(1, 0.7, 1);
    glColor3f(0, 1, 0);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    glPopMatrix();
    
    glDisableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
    glDisableClientState(GL_COLOR_ARRAY); // enable the color array on the client side
}

void displayWalls() {
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //wall px
    glm::mat4 pxMatrix = glm::mat4(1.0f);
    
    pxMatrix = glm::rotate(pxMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pxMatrix = glm::rotate(pxMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pxMatrix = glm::scale(pxMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pxMatrix = glm::translate(pxMatrix, glm::vec3(-(WALLSCALAR + WALLTHICK), WALLSCALAR, 0.0));
    pxMatrix = glm::rotate(pxMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    pxMatrix = glm::rotate(pxMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    pxMatrix = glm::scale(pxMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    pxMatrix = glm::rotate(pxMatrix, 180.0f, glm::vec3(1.0, 0.0, 0.0));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    pxMatrix = viewMatrix * pxMatrix;
    glm::mat4 pxModelviewMatrix = pxMatrix;
    glm::mat3 pxNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pxMatrix)));
    pxMatrix = projectionMatrix * pxMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pxModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pxNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pxMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[1]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall py
    glm::mat4 pyMatrix = glm::mat4(1.0f);
    
    pyMatrix = glm::rotate(pyMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pyMatrix = glm::rotate(pyMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pyMatrix = glm::scale(pyMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pyMatrix = glm::translate(pyMatrix, glm::vec3(0.0, WALLSCALAR * 2 + WALLTHICK, 0.0));
    pyMatrix = glm::rotate(pyMatrix, 180.0f, glm::vec3(0.0, 0.0, 1.0));
    pyMatrix = glm::rotate(pyMatrix, -90.0f, glm::vec3(0.0, 1.0, 0.0));
    pyMatrix = glm::scale(pyMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    pyMatrix = viewMatrix * pyMatrix;
    glm::mat4 pyModelviewMatrix = pyMatrix;
    glm::mat3 pyNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pyMatrix)));
    pyMatrix = projectionMatrix * pyMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pyModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pyNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pyMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[2]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall pz
    glm::mat4 pzMatrix = glm::mat4(1.0f);
    
    pzMatrix = glm::rotate(pzMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pzMatrix = glm::rotate(pzMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pzMatrix = glm::scale(pzMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pzMatrix = glm::translate(pzMatrix, glm::vec3(0.0, WALLSCALAR, -WALLSCALAR - WALLTHICK));
    pzMatrix = glm::rotate(pzMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    pzMatrix = glm::rotate(pzMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    pzMatrix = glm::rotate(pzMatrix, -90.0f, glm::vec3(1.0, 0.0, 0.0));
    pzMatrix = glm::scale(pzMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    pzMatrix = viewMatrix * pzMatrix;
    glm::mat4 pzModelviewMatrix = pzMatrix;
    glm::mat3 pzNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pzMatrix)));
    pzMatrix = projectionMatrix * pzMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pzModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pzNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pzMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[3]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall nx
    glm::mat4 nxMatrix = glm::mat4(1.0f);
    
    nxMatrix = glm::rotate(nxMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nxMatrix = glm::rotate(nxMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nxMatrix = glm::scale(nxMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nxMatrix = glm::translate(nxMatrix, glm::vec3(WALLSCALAR + WALLTHICK, WALLSCALAR, 0.0));
    nxMatrix = glm::rotate(nxMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    nxMatrix = glm::rotate(nxMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    //    nxMatrix = glm::rotate(nxMatrix, 180.0f, glm::vec3(1.0, 0.0, 0.0));
    nxMatrix = glm::scale(nxMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    nxMatrix = viewMatrix * nxMatrix;
    glm::mat4 nxModelviewMatrix = nxMatrix;
    glm::mat3 nxNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nxMatrix)));
    nxMatrix = projectionMatrix * nxMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nxModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nxNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nxMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[4]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall ny
    glm::mat4 nyMatrix = glm::mat4(1.0f);
    
    nyMatrix = glm::rotate(nyMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nyMatrix = glm::rotate(nyMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nyMatrix = glm::scale(nyMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nyMatrix = glm::translate(nyMatrix, glm::vec3(0.0, -WALLTHICK, 0.0));
    //    nyMatrix = glm::rotate(nyMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    nyMatrix = glm::rotate(nyMatrix, 90.0f, glm::vec3(0.0, 1.0, 0.0));
    nyMatrix = glm::scale(nyMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    nyMatrix = viewMatrix * nyMatrix;
    glm::mat4 nyModelviewMatrix = nyMatrix;
    glm::mat3 nyNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nyMatrix)));
    nyMatrix = projectionMatrix * nyMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nyModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nyNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nyMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[5]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall nz
    glm::mat4 nzMatrix = glm::mat4(1.0f);
    
    nzMatrix = glm::rotate(nzMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nzMatrix = glm::rotate(nzMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nzMatrix = glm::scale(nzMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nzMatrix = glm::translate(nzMatrix, glm::vec3(0.0, WALLSCALAR, WALLSCALAR + WALLTHICK));
    nzMatrix = glm::rotate(nzMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    nzMatrix = glm::rotate(nzMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    nzMatrix = glm::rotate(nzMatrix, 90.0f, glm::vec3(1.0, 0.0, 0.0));
    nzMatrix = glm::scale(nzMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    nzMatrix = viewMatrix * nzMatrix;
    glm::mat4 nzModelviewMatrix = nzMatrix;
    glm::mat3 nzNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nzMatrix)));
    nzMatrix = projectionMatrix * nzMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nzModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nzNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nzMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[6]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
}
void displayWalls1() {
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //wall px
    glm::mat4 pxMatrix = glm::mat4(1.0f);
    
    pxMatrix = glm::rotate(pxMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pxMatrix = glm::rotate(pxMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pxMatrix = glm::scale(pxMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pxMatrix = glm::translate(pxMatrix, glm::vec3(-(WALLSCALAR + WALLTHICK), WALLSCALAR, 0.0));
    pxMatrix = glm::rotate(pxMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    pxMatrix = glm::rotate(pxMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    pxMatrix = glm::scale(pxMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    pxMatrix = glm::rotate(pxMatrix, 180.0f, glm::vec3(1.0, 0.0, 0.0));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    pxMatrix = viewMatrix * pxMatrix;
    glm::mat4 pxModelviewMatrix = pxMatrix;
    glm::mat3 pxNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pxMatrix)));
    pxMatrix = projectionMatrix * pxMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pxModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pxNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pxMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[1]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall py
    glm::mat4 pyMatrix = glm::mat4(1.0f);
    
    pyMatrix = glm::rotate(pyMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pyMatrix = glm::rotate(pyMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pyMatrix = glm::scale(pyMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pyMatrix = glm::translate(pyMatrix, glm::vec3(0.0, WALLSCALAR * 2 + WALLTHICK, 0.0));
    pyMatrix = glm::rotate(pyMatrix, 180.0f, glm::vec3(0.0, 0.0, 1.0));
    pyMatrix = glm::rotate(pyMatrix, -90.0f, glm::vec3(0.0, 1.0, 0.0));
    pyMatrix = glm::scale(pyMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    pyMatrix = viewMatrix * pyMatrix;
    glm::mat4 pyModelviewMatrix = pyMatrix;
    glm::mat3 pyNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pyMatrix)));
    pyMatrix = projectionMatrix * pyMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pyModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pyNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pyMatrix[0][0]);
    
    glBindTexture(GL_TEXTURE_2D, textureId[2]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
//    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall pz
    glm::mat4 pzMatrix = glm::mat4(1.0f);
    
    pzMatrix = glm::rotate(pzMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pzMatrix = glm::rotate(pzMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pzMatrix = glm::scale(pzMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pzMatrix = glm::translate(pzMatrix, glm::vec3(0.0, WALLSCALAR, -WALLSCALAR - WALLTHICK));
    pzMatrix = glm::rotate(pzMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    pzMatrix = glm::rotate(pzMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    pzMatrix = glm::rotate(pzMatrix, -90.0f, glm::vec3(1.0, 0.0, 0.0));
    pzMatrix = glm::scale(pzMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    pzMatrix = viewMatrix * pzMatrix;
    glm::mat4 pzModelviewMatrix = pzMatrix;
    glm::mat3 pzNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pzMatrix)));
    pzMatrix = projectionMatrix * pzMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pzModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pzNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pzMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[3]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall nx
    glm::mat4 nxMatrix = glm::mat4(1.0f);
    
    nxMatrix = glm::rotate(nxMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nxMatrix = glm::rotate(nxMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nxMatrix = glm::scale(nxMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nxMatrix = glm::translate(nxMatrix, glm::vec3(WALLSCALAR + WALLTHICK, WALLSCALAR, 0.0));
    nxMatrix = glm::rotate(nxMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    nxMatrix = glm::rotate(nxMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    //    nxMatrix = glm::rotate(nxMatrix, 180.0f, glm::vec3(1.0, 0.0, 0.0));
    nxMatrix = glm::scale(nxMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    nxMatrix = viewMatrix * nxMatrix;
    glm::mat4 nxModelviewMatrix = nxMatrix;
    glm::mat3 nxNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nxMatrix)));
    nxMatrix = projectionMatrix * nxMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nxModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nxNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nxMatrix[0][0]);
    glUniform1i(flagNormalMapShader, false);

    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[4]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    glUniform1i(flagNormalMapShader, false);

    //wall ny: floor
    glm::mat4 nyMatrix = glm::mat4(1.0f);
    
    nyMatrix = glm::rotate(nyMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nyMatrix = glm::rotate(nyMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nyMatrix = glm::scale(nyMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nyMatrix = glm::translate(nyMatrix, glm::vec3(0.0, -WALLTHICK, 0.0));
    nyMatrix = glm::rotate(nyMatrix, 90.0f, glm::vec3(0.0, 1.0, 0.0));
    nyMatrix = glm::scale(nyMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    glm::mat4 nyShadowMatrix = shadowMatrix * nyMatrix;
    nyMatrix = viewMatrix * nyMatrix;
    glm::mat4 nyModelviewMatrix = nyMatrix;
    glm::mat3 nyNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nyMatrix)));
    nyMatrix = projectionMatrix * nyMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nyModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nyNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nyMatrix[0][0]);
    glUniformMatrix4fv(shadowMatrixShader, 1, GL_FALSE, &nyShadowMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glActiveTexture(GL_TEXTURE0);
    
    if (flagEdgeDetection == true) {
        glBindTexture(GL_TEXTURE_2D, depthTex); // bind this handle to a 2D texture
    }
    else {
        glBindTexture(GL_TEXTURE_2D, textureId[5]); // bind this handle to a 2D texture
    }
    glActiveTexture(GL_TEXTURE5);
    glBindTexture(GL_TEXTURE_2D, depthTex);
    
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(flagTextureShader, true);
    glUniform1i(flagShadowMapShader, true);
    glUniform1i(flagNormalMapShader, false);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagNormalMapShader, false);
    glUniform1i(flagShadowMapShader, false);
    glUniform1i(flagTextureShader, false);
    //TODO: shadow not only on ground, but also on desktop.
    //wall nz
    glm::mat4 nzMatrix = glm::mat4(1.0f);
    
    nzMatrix = glm::rotate(nzMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nzMatrix = glm::rotate(nzMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nzMatrix = glm::scale(nzMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nzMatrix = glm::translate(nzMatrix, glm::vec3(0.0, WALLSCALAR, WALLSCALAR + WALLTHICK));
    nzMatrix = glm::rotate(nzMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    nzMatrix = glm::rotate(nzMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    nzMatrix = glm::rotate(nzMatrix, 90.0f, glm::vec3(1.0, 0.0, 0.0));
    nzMatrix = glm::scale(nzMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    nzMatrix = viewMatrix * nzMatrix;
    glm::mat4 nzModelviewMatrix = nzMatrix;
    glm::mat3 nzNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nzMatrix)));
    nzMatrix = projectionMatrix * nzMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nzModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nzNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nzMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[6]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
}
void displayEdge() {
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
/* other walls
    //wall px
    glm::mat4 pxMatrix = glm::mat4(1.0f);
    
    pxMatrix = glm::rotate(pxMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pxMatrix = glm::rotate(pxMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pxMatrix = glm::scale(pxMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pxMatrix = glm::translate(pxMatrix, glm::vec3(-(WALLSCALAR + WALLTHICK), WALLSCALAR, 0.0));
    pxMatrix = glm::rotate(pxMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    pxMatrix = glm::rotate(pxMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    pxMatrix = glm::scale(pxMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    pxMatrix = glm::rotate(pxMatrix, 180.0f, glm::vec3(1.0, 0.0, 0.0));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    pxMatrix = viewMatrix * pxMatrix;
    glm::mat4 pxModelviewMatrix = pxMatrix;
    glm::mat3 pxNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pxMatrix)));
    pxMatrix = projectionMatrix * pxMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pxModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pxNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pxMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[1]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall py
    glm::mat4 pyMatrix = glm::mat4(1.0f);
    
    pyMatrix = glm::rotate(pyMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pyMatrix = glm::rotate(pyMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pyMatrix = glm::scale(pyMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pyMatrix = glm::translate(pyMatrix, glm::vec3(0.0, WALLSCALAR * 2 + WALLTHICK, 0.0));
    pyMatrix = glm::rotate(pyMatrix, 180.0f, glm::vec3(0.0, 0.0, 1.0));
    pyMatrix = glm::rotate(pyMatrix, -90.0f, glm::vec3(0.0, 1.0, 0.0));
    pyMatrix = glm::scale(pyMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    pyMatrix = viewMatrix * pyMatrix;
    glm::mat4 pyModelviewMatrix = pyMatrix;
    glm::mat3 pyNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pyMatrix)));
    pyMatrix = projectionMatrix * pyMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pyModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pyNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pyMatrix[0][0]);
    
    glBindTexture(GL_TEXTURE_2D, textureId[2]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    //    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall pz
    glm::mat4 pzMatrix = glm::mat4(1.0f);
    
    pzMatrix = glm::rotate(pzMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    pzMatrix = glm::rotate(pzMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    pzMatrix = glm::scale(pzMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    pzMatrix = glm::translate(pzMatrix, glm::vec3(0.0, WALLSCALAR, -WALLSCALAR - WALLTHICK));
    pzMatrix = glm::rotate(pzMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    pzMatrix = glm::rotate(pzMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    pzMatrix = glm::rotate(pzMatrix, -90.0f, glm::vec3(1.0, 0.0, 0.0));
    pzMatrix = glm::scale(pzMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    pzMatrix = viewMatrix * pzMatrix;
    glm::mat4 pzModelviewMatrix = pzMatrix;
    glm::mat3 pzNormalMatrix = glm::transpose(glm::inverse(glm::mat3(pzMatrix)));
    pzMatrix = projectionMatrix * pzMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &pzModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &pzNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &pzMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[3]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    //wall nx
    glm::mat4 nxMatrix = glm::mat4(1.0f);
    
    nxMatrix = glm::rotate(nxMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nxMatrix = glm::rotate(nxMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nxMatrix = glm::scale(nxMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nxMatrix = glm::translate(nxMatrix, glm::vec3(WALLSCALAR + WALLTHICK, WALLSCALAR, 0.0));
    nxMatrix = glm::rotate(nxMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    nxMatrix = glm::rotate(nxMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    //    nxMatrix = glm::rotate(nxMatrix, 180.0f, glm::vec3(1.0, 0.0, 0.0));
    nxMatrix = glm::scale(nxMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    //    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(1.0, 0.01, 0.3));
    nxMatrix = viewMatrix * nxMatrix;
    glm::mat4 nxModelviewMatrix = nxMatrix;
    glm::mat3 nxNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nxMatrix)));
    nxMatrix = projectionMatrix * nxMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nxModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nxNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nxMatrix[0][0]);
    glUniform1i(flagNormalMapShader, false);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[4]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    glUniform1i(flagNormalMapShader, false);
*/
    //wall ny: floor
    glm::mat4 nyMatrix = glm::mat4(1.0f);
    
    nyMatrix = glm::rotate(nyMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nyMatrix = glm::rotate(nyMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nyMatrix = glm::scale(nyMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nyMatrix = glm::translate(nyMatrix, glm::vec3(0.0, 5.0, 0.0));
    nyMatrix = glm::rotate(nyMatrix, 90.0f, glm::vec3(0.0, 0.0, -1.0));
    nyMatrix = glm::rotate(nyMatrix, 26.56f, glm::vec3(1.0, 0.0, 0.0));
    nyMatrix = glm::rotate(nyMatrix, 11.537f, glm::vec3(0.0, 0.0, 1.0));
    nyMatrix = glm::rotate(nyMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    nyMatrix = glm::translate(nyMatrix, glm::vec3(0.0, -2.5, 0.0));
    nyMatrix = glm::scale(nyMatrix, glm::vec3(14.1421, 1.0, 14.1421));
    glm::mat4 nyShadowMatrix = shadowMatrix * nyMatrix;
    nyMatrix = viewMatrix * nyMatrix;
    glm::mat4 nyModelviewMatrix = nyMatrix;
    glm::mat3 nyNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nyMatrix)));
    nyMatrix = projectionMatrix * nyMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nyModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nyNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nyMatrix[0][0]);
    glUniformMatrix4fv(shadowMatrixShader, 1, GL_FALSE, &nyShadowMatrix[0][0]);
    glActiveTexture(GL_TEXTURE0);
    if (flagEdgeDetection == true) {
        glBindTexture(GL_TEXTURE_2D, depthTex); // bind this handle to a 2D texture
    }
    else {
        glBindTexture(GL_TEXTURE_2D, textureId[5]); // bind this handle to a 2D texture
    }
    glActiveTexture(GL_TEXTURE5);
    glBindTexture(GL_TEXTURE_2D, depthTex);
    
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(flagTextureShader, true);
    glUniform1i(flagShadowMapShader, true);
    glUniform1i(flagNormalMapShader, false);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagNormalMapShader, false);
    glUniform1i(flagShadowMapShader, false);
    glUniform1i(flagTextureShader, false);
    
/* other walls
    //TODO: shadow not only on ground, but also on desktop.
    //wall nz
    glm::mat4 nzMatrix = glm::mat4(1.0f);
    
    nzMatrix = glm::rotate(nzMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    nzMatrix = glm::rotate(nzMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    nzMatrix = glm::scale(nzMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    nzMatrix = glm::translate(nzMatrix, glm::vec3(0.0, WALLSCALAR, WALLSCALAR + WALLTHICK));
    nzMatrix = glm::rotate(nzMatrix, 90.0f, glm::vec3(0.0, 0.0, 1.0));
    nzMatrix = glm::rotate(nzMatrix, 180.0f, glm::vec3(0.0, 1.0, 0.0));
    nzMatrix = glm::rotate(nzMatrix, 90.0f, glm::vec3(1.0, 0.0, 0.0));
    nzMatrix = glm::scale(nzMatrix, glm::vec3(WALLSCALAR, WALLTHICK, WALLSCALAR));
    nzMatrix = viewMatrix * nzMatrix;
    glm::mat4 nzModelviewMatrix = nzMatrix;
    glm::mat3 nzNormalMatrix = glm::transpose(glm::inverse(glm::mat3(nzMatrix)));
    nzMatrix = projectionMatrix * nzMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &nzModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &nzNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &nzMatrix[0][0]);
    
    //    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(viewMatrix)));
    //    glUniformMatrix4fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    glBindTexture(GL_TEXTURE_2D, textureId[6]); // bind this handle to a 2D texture
    glUniform1i(flagTextureShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
 */
}

void displayScene() {
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+48);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //floor
    glm::mat4 floorMatrix = glm::mat4(1.0f);
    
    floorMatrix = glm::rotate(floorMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    floorMatrix = glm::rotate(floorMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    floorMatrix = glm::scale(floorMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    //    floorMatrix = glm::translate(floorMatrix, glm::vec3(0.0, -1.0, -150.0));
    //    floorMatrix = glm::scale(floorMatrix, glm::vec3(260.0, 1.0, 260.0));
    floorMatrix = glm::scale(floorMatrix, glm::vec3(260.0, 0.01, 260.0));
    floorMatrix = viewMatrix * floorMatrix;
    glm::mat4 floorModelviewMatrix = floorMatrix;
    glm::mat3 floorNormalMatrix = glm::transpose(glm::inverse(glm::mat3(floorMatrix)));
    floorMatrix = projectionMatrix * floorMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &floorModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &floorNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &floorMatrix[0][0]);
    //    glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0);
    
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+48);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    glVertexAttribPointer(tangentShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 88);

    //desk leg 1
    glm::mat4 leg1Matrix = glm::mat4(1.0f);
    
    leg1Matrix = glm::rotate(leg1Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg1Matrix = glm::rotate(leg1Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(scale_size, scale_size, scale_size));
    
    leg1Matrix = glm::translate(leg1Matrix, glm::vec3(-2.0, 2.0, 2.0) + deskMove);
    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(0.3, 2.0, 0.3));
    leg1Matrix = viewMatrix * leg1Matrix;
    glm::mat4 leg1ModelviewMatrix = leg1Matrix;
    glm::mat3 leg1NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg1Matrix)));
    leg1Matrix = projectionMatrix * leg1Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg1ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg1NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg1Matrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //desk leg 2
    glm::mat4 leg2Matrix = glm::mat4(1.0f);
    
    leg2Matrix = glm::rotate(leg2Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg2Matrix = glm::rotate(leg2Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg2Matrix = glm::scale(leg2Matrix, glm::vec3(scale_size, scale_size, scale_size));
    
    leg2Matrix = glm::translate(leg2Matrix, glm::vec3(-2.0, 2.0, -2.0) + deskMove);
    leg2Matrix = glm::scale(leg2Matrix, glm::vec3(0.3, 2.0, 0.3));
    leg2Matrix = viewMatrix * leg2Matrix;
    glm::mat4 leg2ModelviewMatrix = leg2Matrix;
    glm::mat3 leg2NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg2Matrix)));
    leg2Matrix = projectionMatrix * leg2Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg2ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg2NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg2Matrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //desk leg 3
    glm::mat4 leg3Matrix = glm::mat4(1.0f);
    
    leg3Matrix = glm::rotate(leg3Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg3Matrix = glm::rotate(leg3Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg3Matrix = glm::scale(leg3Matrix, glm::vec3(scale_size, scale_size, scale_size));
    
    leg3Matrix = glm::translate(leg3Matrix, glm::vec3(2.0, 2.0, -2.0) + deskMove);
    leg3Matrix = glm::scale(leg3Matrix, glm::vec3(0.3, 2.0, 0.3));
    leg3Matrix = viewMatrix * leg3Matrix;
    glm::mat4 leg3ModelviewMatrix = leg3Matrix;
    glm::mat3 leg3NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg3Matrix)));
    leg3Matrix = projectionMatrix * leg3Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg3ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg3NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg3Matrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //desk leg 4
    glm::mat4 leg4Matrix = glm::mat4(1.0f);
    
    leg4Matrix = glm::rotate(leg4Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg4Matrix = glm::rotate(leg4Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg4Matrix = glm::scale(leg4Matrix, glm::vec3(scale_size, scale_size, scale_size));
    
    leg4Matrix = glm::translate(leg4Matrix, glm::vec3(2.0, 2.0, 2.0) + deskMove);
    leg4Matrix = glm::scale(leg4Matrix, glm::vec3(0.3, 2.0, 0.3));
    leg4Matrix = viewMatrix * leg4Matrix;
    glm::mat4 leg4ModelviewMatrix = leg4Matrix;
    glm::mat3 leg4NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg4Matrix)));
    leg4Matrix = projectionMatrix * leg4Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg4ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg4NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg4Matrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //desk top
    glm::mat4 deskTopMatrix = glm::mat4(1.0f);
    
    deskTopMatrix = glm::rotate(deskTopMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    deskTopMatrix = glm::rotate(deskTopMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    deskTopMatrix = glm::scale(deskTopMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    deskTopMatrix = glm::translate(deskTopMatrix, glm::vec3(0.0, 4.5, 0.0) + deskMove);
    deskTopMatrix = glm::scale(deskTopMatrix, glm::vec3(3.0, 0.5, 3.0));
    deskTopMatrix = viewMatrix * deskTopMatrix;
    glm::mat4 deskTopModelviewMatrix = deskTopMatrix;
    glm::mat3 deskTopNormalMatrix = glm::transpose(glm::inverse(glm::mat3(deskTopMatrix)));
    deskTopMatrix = projectionMatrix * deskTopMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &deskTopModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &deskTopNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &deskTopMatrix[0][0]);
    glUniform1i(flagTextureShader, true);
    glBindTexture(GL_TEXTURE_2D, textureId[0]); // bind this handle to a 2D texture
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    glUniform1i(flagTextureShader, false);
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[0]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+48);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    
    //the base of lamp
    glm::mat4 lampBaseMatrix = glm::mat4(1.0f);
    
    lampBaseMatrix = glm::rotate(lampBaseMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    lampBaseMatrix = glm::rotate(lampBaseMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    lampBaseMatrix = glm::scale(lampBaseMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    lampBaseMatrix = glm::translate(lampBaseMatrix, glm::vec3(-1.0, 5.5, 1.0) + deskMove);
    lampBaseMatrix = glm::rotate(lampBaseMatrix, 90.0f, glm::vec3(1.0, 0.0, 0.0));
    lampBaseMatrix = glm::scale(lampBaseMatrix, glm::vec3(1.0, 1.0, 0.5));
    lampBaseMatrix = viewMatrix * lampBaseMatrix;
    glm::mat4 lampBaseModelviewMatrix = lampBaseMatrix;
    glm::mat3 lampBaseNormalMatrix = glm::transpose(glm::inverse(glm::mat3(lampBaseMatrix)));
    lampBaseMatrix = projectionMatrix * lampBaseMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &lampBaseModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &lampBaseNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &lampBaseMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //the middle of lamp
    glm::mat4 lampMiddleMatrix = glm::mat4(1.0f);
    
    lampMiddleMatrix = glm::rotate(lampMiddleMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    lampMiddleMatrix = glm::rotate(lampMiddleMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    lampMiddleMatrix = glm::scale(lampMiddleMatrix,glm::vec3(scale_size, scale_size, scale_size));
    
    lampMiddleMatrix = glm::translate(lampMiddleMatrix, glm::vec3(-1.0, 7.0, 1.0) + deskMove);
    lampMiddleMatrix = glm::rotate(lampMiddleMatrix, 90.0f, glm::vec3(1.0, 0.0, 0.0));
    lampMiddleMatrix = glm::scale(lampMiddleMatrix, glm::vec3(0.1, 0.1, 1.5));
    lampMiddleMatrix = viewMatrix * lampMiddleMatrix;
    glm::mat4 lampMiddleModelviewMatrix = lampMiddleMatrix;
    glm::mat3 lampMiddleNormalMatrix = glm::transpose(glm::inverse(glm::mat3(lampMiddleMatrix)));
    lampMiddleMatrix = projectionMatrix * lampMiddleMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &lampMiddleModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &lampMiddleNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &lampMiddleMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0);
    
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[2]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[2]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //top of lamp
    glm::mat4 lampTopMatrix = glm::mat4(1.0f);
    
    lampTopMatrix = glm::rotate(lampTopMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    lampTopMatrix = glm::rotate(lampTopMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    lampTopMatrix = glm::scale(lampTopMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    lampTopMatrix = glm::translate(lampTopMatrix, glm::vec3(0.0, 6.5, 0.0) + deskMove);
    lampTopMatrix = glm::rotate(lampTopMatrix, -45.0f, glm::vec3(0.0, 1.0, 0.0));
    lampTopMatrix = glm::rotate(lampTopMatrix, -20.0f, glm::vec3(1.0, 0.0, 0.0));
    lampTopMatrix = glm::scale(lampTopMatrix, glm::vec3(0.5, 0.5, 1.5));
    lampTopMatrix = viewMatrix * lampTopMatrix;
    glm::mat4 lampTopModelviewMatrix = lampTopMatrix;
    glm::mat3 lampTopNormalMatrix = glm::transpose(glm::inverse(glm::mat3(lampTopMatrix)));
    lampTopMatrix = projectionMatrix * lampTopMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &lampTopModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &lampTopNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &lampTopMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, coneIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[3]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[3]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+48);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //the pot on the desk
    glm::mat4 potMatrix = glm::mat4(1.0f);
    
    potMatrix = glm::rotate(potMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    potMatrix = glm::rotate(potMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    potMatrix = glm::scale(potMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    potMatrix = glm::translate(potMatrix, glm::vec3(0.0, 5.0, 0.0));
    potMatrix = glm::scale(potMatrix, glm::vec3(2.0, 2.0, 2.0));
    potMatrix = viewMatrix * potMatrix;
    glm::mat4 potModelviewMatrix = potMatrix;
    glm::mat3 potNormalMatrix = glm::transpose(glm::inverse(glm::mat3(potMatrix)));
    potMatrix = projectionMatrix * potMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &potModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &potNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &potMatrix[0][0]);
    
    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(potModelviewMatrix)));
    glUniformMatrix3fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    
    glUniform1i(flagCubemapShader, true);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    glUniform1i(flagCubemapShader, false);
}
void displayScene1() {
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+48);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //floor
    glm::mat4 floorMatrix = glm::mat4(1.0f);
    
    floorMatrix = glm::rotate(floorMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    floorMatrix = glm::rotate(floorMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    floorMatrix = glm::scale(floorMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    //    floorMatrix = glm::translate(floorMatrix, glm::vec3(0.0, -1.0, -150.0));
    //    floorMatrix = glm::scale(floorMatrix, glm::vec3(260.0, 1.0, 260.0));
    floorMatrix = glm::scale(floorMatrix, glm::vec3(260.0, 0.01, 260.0));
    floorMatrix = viewMatrix * floorMatrix;
    glm::mat4 floorModelviewMatrix = floorMatrix;
    glm::mat3 floorNormalMatrix = glm::transpose(glm::inverse(glm::mat3(floorMatrix)));
    floorMatrix = projectionMatrix * floorMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &floorModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &floorNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &floorMatrix[0][0]);
    //    glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0);
    
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[1]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+48);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //desk leg 1
    glm::mat4 leg1Matrix = glm::mat4(1.0f);
    
    leg1Matrix = glm::rotate(leg1Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg1Matrix = glm::rotate(leg1Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(scale_size, scale_size, scale_size));
    

    
    leg1Matrix = glm::translate(leg1Matrix, glm::vec3(-2.0, 2.0, 2.0) + deskMove);
    leg1Matrix = glm::scale(leg1Matrix, glm::vec3(0.3, 2.0, 0.3));
    leg1Matrix = viewMatrix * leg1Matrix;
//    leg1Matrix = trackball * leg1Matrix;
    glm::mat4 leg1ModelviewMatrix = leg1Matrix;
    glm::mat3 leg1NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg1Matrix)));
    leg1Matrix = projectionMatrix * leg1Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg1ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg1NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg1Matrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //desk leg 2
    glm::mat4 leg2Matrix = glm::mat4(1.0f);
    
    leg2Matrix = glm::rotate(leg2Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg2Matrix = glm::rotate(leg2Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg2Matrix = glm::scale(leg2Matrix, glm::vec3(scale_size, scale_size, scale_size));
    
    leg2Matrix = glm::translate(leg2Matrix, glm::vec3(-2.0, 2.0, -2.0) + deskMove);
    leg2Matrix = glm::scale(leg2Matrix, glm::vec3(0.3, 2.0, 0.3));
    leg2Matrix = viewMatrix * leg2Matrix;
    glm::mat4 leg2ModelviewMatrix = leg2Matrix;
    glm::mat3 leg2NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg2Matrix)));
    leg2Matrix = projectionMatrix * leg2Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg2ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg2NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg2Matrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //desk leg 3
    glm::mat4 leg3Matrix = glm::mat4(1.0f);
    
    leg3Matrix = glm::rotate(leg3Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg3Matrix = glm::rotate(leg3Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg3Matrix = glm::scale(leg3Matrix, glm::vec3(scale_size, scale_size, scale_size));
    
    leg3Matrix = glm::translate(leg3Matrix, glm::vec3(2.0, 2.0, -2.0) + deskMove);
    leg3Matrix = glm::scale(leg3Matrix, glm::vec3(0.3, 2.0, 0.3));
    leg3Matrix = viewMatrix * leg3Matrix;
    glm::mat4 leg3ModelviewMatrix = leg3Matrix;
    glm::mat3 leg3NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg3Matrix)));
    leg3Matrix = projectionMatrix * leg3Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg3ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg3NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg3Matrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //desk leg 4
    glm::mat4 leg4Matrix = glm::mat4(1.0f);
    
    leg4Matrix = glm::rotate(leg4Matrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    leg4Matrix = glm::rotate(leg4Matrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    leg4Matrix = glm::scale(leg4Matrix, glm::vec3(scale_size, scale_size, scale_size));
    
    leg4Matrix = glm::translate(leg4Matrix, glm::vec3(2.0, 2.0, 2.0) + deskMove);
    leg4Matrix = glm::scale(leg4Matrix, glm::vec3(0.3, 2.0, 0.3));
    leg4Matrix = viewMatrix * leg4Matrix;
    glm::mat4 leg4ModelviewMatrix = leg4Matrix;
    glm::mat3 leg4NormalMatrix = glm::transpose(glm::inverse(glm::mat3(leg4Matrix)));
    leg4Matrix = projectionMatrix * leg4Matrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leg4ModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leg4NormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leg4Matrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //desk top
    glm::mat4 deskTopMatrix = glm::mat4(1.0f);
    
    deskTopMatrix = glm::rotate(deskTopMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    deskTopMatrix = glm::rotate(deskTopMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    deskTopMatrix = glm::scale(deskTopMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    deskTopMatrix = glm::translate(deskTopMatrix, glm::vec3(0.0, 4.5, 0.0) + deskMove);
    deskTopMatrix = glm::scale(deskTopMatrix, glm::vec3(3.0, 0.5, 3.0));
    glm::mat4 deskTopShadowMatrix = shadowMatrix * deskTopMatrix;
    deskTopMatrix = viewMatrix * deskTopMatrix;
    glm::mat4 deskTopModelviewMatrix = deskTopMatrix;
    glm::mat3 deskTopNormalMatrix = glm::transpose(glm::inverse(glm::mat3(deskTopMatrix)));
    deskTopMatrix = projectionMatrix * deskTopMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &deskTopModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &deskTopNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &deskTopMatrix[0][0]);
    glUniformMatrix4fv(shadowMatrixShader, 1, GL_FALSE, &deskTopShadowMatrix[0][0]);
    glUniform1i(flagTextureShader, true);
    glBindTexture(GL_TEXTURE_2D, textureId[0]); // bind this handle to a 2D texture
    
    glActiveTexture(GL_TEXTURE5);
    glUniform1i(flagShadowMapShader, false);
    glActiveTexture(GL_TEXTURE0);
    //TODO: not only display bump mapping on one surface of cube, but also on other primitives
    glUniform1i(flagNormalMapShader, false);

    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    glUniform1i(flagNormalMapShader, false);
    
    glUniform1i(flagShadowMapShader, false);
    
    glUniform1i(flagTextureShader, false);
    
    
    
    //desk top: bump mapping
    deskTopMatrix = glm::mat4(1.0f);
    
    deskTopMatrix = glm::rotate(deskTopMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    deskTopMatrix = glm::rotate(deskTopMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    deskTopMatrix = glm::scale(deskTopMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    deskTopMatrix = glm::translate(deskTopMatrix, glm::vec3(0.0, 5.01, 0.0) + deskMove);
    deskTopMatrix = glm::scale(deskTopMatrix, glm::vec3(3.0, 0.01, 3.0));
    deskTopShadowMatrix = shadowMatrix * deskTopMatrix;
    deskTopMatrix = viewMatrix * deskTopMatrix;
    deskTopModelviewMatrix = deskTopMatrix;
    deskTopNormalMatrix = glm::transpose(glm::inverse(glm::mat3(deskTopMatrix)));
    deskTopMatrix = projectionMatrix * deskTopMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &deskTopModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &deskTopNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &deskTopMatrix[0][0]);
    glUniformMatrix4fv(shadowMatrixShader, 1, GL_FALSE, &deskTopShadowMatrix[0][0]);
    glUniform1i(flagTextureShader, true);
    glBindTexture(GL_TEXTURE_2D, textureId[0]); // bind this handle to a 2D texture
    
    glActiveTexture(GL_TEXTURE5);
    glUniform1i(flagShadowMapShader, false);
    glActiveTexture(GL_TEXTURE0);
    //TODO: not only display bump mapping on one surface of cube, but also on other primitives
    glUniform1i(flagNormalMapShader, true);
    
    glDrawElements(GL_TRIANGLE_STRIP, cubeIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
    
    glUniform1i(flagNormalMapShader, false);
    
    glUniform1i(flagShadowMapShader, false);
    
    glUniform1i(flagTextureShader, false);
    
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[0]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+48);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    
    //the base of lamp
    glm::mat4 lampBaseMatrix = glm::mat4(1.0f);
    
    lampBaseMatrix = glm::rotate(lampBaseMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    lampBaseMatrix = glm::rotate(lampBaseMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    lampBaseMatrix = glm::scale(lampBaseMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    lampBaseMatrix = glm::translate(lampBaseMatrix, glm::vec3(-1.0, 5.5, 1.0) + deskMove);
    lampBaseMatrix = glm::rotate(lampBaseMatrix, 90.0f, glm::vec3(1.0, 0.0, 0.0));
    lampBaseMatrix = glm::scale(lampBaseMatrix, glm::vec3(1.0, 1.0, 0.5));
    lampBaseMatrix = viewMatrix * lampBaseMatrix;
    glm::mat4 lampBaseModelviewMatrix = lampBaseMatrix;
    glm::mat3 lampBaseNormalMatrix = glm::transpose(glm::inverse(glm::mat3(lampBaseMatrix)));
    lampBaseMatrix = projectionMatrix * lampBaseMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &lampBaseModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &lampBaseNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &lampBaseMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0);
    
    //the middle of lamp
    glm::mat4 lampMiddleMatrix = glm::mat4(1.0f);
    
    lampMiddleMatrix = glm::rotate(lampMiddleMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    lampMiddleMatrix = glm::rotate(lampMiddleMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    lampMiddleMatrix = glm::scale(lampMiddleMatrix,glm::vec3(scale_size, scale_size, scale_size));
    
    lampMiddleMatrix = glm::translate(lampMiddleMatrix, glm::vec3(-1.0, 7.0, 1.0) + deskMove);
    lampMiddleMatrix = glm::rotate(lampMiddleMatrix, 90.0f, glm::vec3(1.0, 0.0, 0.0));
    lampMiddleMatrix = glm::scale(lampMiddleMatrix, glm::vec3(0.1, 0.1, 1.5));
    lampMiddleMatrix = viewMatrix * lampMiddleMatrix;
    glm::mat4 lampMiddleModelviewMatrix = lampMiddleMatrix;
    glm::mat3 lampMiddleNormalMatrix = glm::transpose(glm::inverse(glm::mat3(lampMiddleMatrix)));
    lampMiddleMatrix = projectionMatrix * lampMiddleMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &lampMiddleModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &lampMiddleNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &lampMiddleMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0);
    
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[2]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[2]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //top of lamp
    glm::mat4 lampTopMatrix = glm::mat4(1.0f);
    
    lampTopMatrix = glm::rotate(lampTopMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    lampTopMatrix = glm::rotate(lampTopMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    lampTopMatrix = glm::scale(lampTopMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    lampTopMatrix = glm::translate(lampTopMatrix, glm::vec3(0.0, 6.5, 0.0) + deskMove);
    lampTopMatrix = glm::rotate(lampTopMatrix, -45.0f, glm::vec3(0.0, 1.0, 0.0));
    lampTopMatrix = glm::rotate(lampTopMatrix, -20.0f, glm::vec3(1.0, 0.0, 0.0));
    lampTopMatrix = glm::scale(lampTopMatrix, glm::vec3(0.5, 0.5, 1.5));
    lampTopMatrix = viewMatrix * lampTopMatrix;
    glm::mat4 lampTopModelviewMatrix = lampTopMatrix;
    glm::mat3 lampTopNormalMatrix = glm::transpose(glm::inverse(glm::mat3(lampTopMatrix)));
    lampTopMatrix = projectionMatrix * lampTopMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &lampTopModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &lampTopNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &lampTopMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, coneIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[3]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[3]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    //the pot on the desk
    glm::mat4 potMatrix = glm::mat4(1.0f);
    
    potMatrix = glm::rotate(potMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    potMatrix = glm::rotate(potMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    potMatrix = glm::scale(potMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    potMatrix = glm::translate(potMatrix, glm::vec3(0.0, 5.0, 0.0));
    potMatrix = glm::scale(potMatrix, glm::vec3(2.0, 2.0, 2.0));
    potMatrix = viewMatrix * potMatrix;
    glm::mat4 potModelviewMatrix = potMatrix;
    glm::mat3 potNormalMatrix = glm::transpose(glm::inverse(glm::mat3(potMatrix)));
    potMatrix = projectionMatrix * potMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &potModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &potNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &potMatrix[0][0]);
    
    glm::mat3 inv = glm::transpose(glm::inverse(glm::mat3(potModelviewMatrix)));
    glUniformMatrix3fv(inverseModelviewShader, 1, GL_FALSE, &inv[0][0]);
    
    if (flagAnimation == false) {
        glUniform1i(flagCubemapShader, true);
        glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
        glUniform1i(flagCubemapShader, false);
    }
    else {
        glUniform1i(flagAnimationShader, flagAnimation);
        glDrawElements(GL_POINTS, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
        glUniform1i(flagAnimationShader, false);
    }
}
void displayRobot() {
    
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[3]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[3]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    glVertexAttribPointer(tangentShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+88);
    
    //torso green center(0, 5, 6)
    glm::mat4 baseMatrix = glm::mat4(1.0f);
    
    //left and right mouse control
    baseMatrix = glm::rotate(baseMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    baseMatrix = glm::rotate(baseMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    baseMatrix = glm::scale(baseMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    baseMatrix = glm::translate(baseMatrix, glm::vec3(0.0, 5.0, 6.0));    //torso center
    baseMatrix = baseMatrix * transformMatrix;
    //use this transformMatrix to move and rotate robot
    
    glm::mat4 torsoMatrix = baseMatrix;
    torsoMatrix = glm::scale(torsoMatrix, glm::vec3(0.5, 1.0, 0.75));
    torsoMatrix = viewMatrix * torsoMatrix;
    glm::mat4 torsoModelviewMatrix = torsoMatrix;
    glm::mat3 torsoNormalMatrix = glm::transpose(glm::inverse(glm::mat3(torsoMatrix)));
    torsoMatrix = projectionMatrix * torsoMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &torsoModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &torsoNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &torsoMatrix[0][0]);
//TODO: change back
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    glUniform1i(flagAnimationShader, flagAnimation);
//    glDrawElements(GL_POINTS, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    glUniform1i(flagAnimationShader, false);

    //shoulder blue center(0, 6.25, 0)
    baseMatrix = glm::translate(baseMatrix, glm::vec3(0.0, 1.25, 0.0));
    glm::mat4 shoulderMatrix = baseMatrix;
    shoulderMatrix = glm::scale(shoulderMatrix, glm::vec3(0.5, 0.25, 1.0));
    shoulderMatrix = viewMatrix * shoulderMatrix;
    glm::mat4 shoulderModelviewMatrix = shoulderMatrix;
    glm::mat3 shoulderNormalMatrix = glm::transpose(glm::inverse(glm::mat3(shoulderMatrix)));
    shoulderMatrix = projectionMatrix * shoulderMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &shoulderModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &shoulderNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &shoulderMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //right upper arm stone-blue center(0, 5.5, 1)
    glm::mat4 rightArmMatrix = baseMatrix;
    if (robotStep == 1) {   //robot walk animation
        rightArmMatrix = glm::rotate(rightArmMatrix, -10.0f, glm::vec3(0.0, 0.0, 1.0));
    }
    if (robotStep == 2) {
        rightArmMatrix = glm::rotate(rightArmMatrix, 10.0f, glm::vec3(0.0, 0.0, 1.0));
    }
    rightArmMatrix = glm::rotate(rightArmMatrix, liftRightArmAngle, glm::vec3(0.0, 0.0, 1.0));
    rightArmMatrix = glm::translate(rightArmMatrix, glm::vec3(0.0, -0.75, 1.0));
    glm::mat4 rightUpperArmMatrix = glm::scale(rightArmMatrix, glm::vec3(0.25, 0.75, 0.25));
    rightUpperArmMatrix = viewMatrix * rightUpperArmMatrix;
    glm::mat4 rightUpperArmModelviewMatrix = rightUpperArmMatrix;
    glm::mat3 rightUpperArmNormalMatrix = glm::transpose(glm::inverse(glm::mat3(rightUpperArmMatrix)));
    rightUpperArmMatrix = projectionMatrix * rightUpperArmMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &rightUpperArmModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &rightUpperArmNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &rightUpperArmMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //right lower arm stone-blue center(0, 4, 1)
    rightArmMatrix = glm::translate(rightArmMatrix, glm::vec3(0.0, -1.5, 0.0));
    glm::mat4 rightLowerArmMatrix = rightArmMatrix;
    rightLowerArmMatrix = glm::scale(rightLowerArmMatrix, glm::vec3(0.25, 0.75, 0.25));
    rightLowerArmMatrix = viewMatrix * rightLowerArmMatrix;
    glm::mat4 rightLowerArmModelviewMatrix = rightLowerArmMatrix;
    glm::mat3 rightLowerArmNormalMatrix = glm::transpose(glm::inverse(glm::mat3(rightLowerArmMatrix)));
    rightLowerArmMatrix = projectionMatrix * rightLowerArmMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &rightLowerArmModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &rightLowerArmNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &rightLowerArmMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //right hand stone-blue center(0, 3.5, 1)
    rightArmMatrix = glm::translate(rightArmMatrix, glm::vec3(0.0, -1.0, 0.0));
    glm::mat4 rightHandMatrix = rightArmMatrix;
    rightHandMatrix = glm::scale(rightHandMatrix, glm::vec3(0.5, 0.25, 0.25));
    rightHandMatrix = viewMatrix * rightHandMatrix;
    glm::mat4 rightHandModelviewMatrix = rightHandMatrix;
    glm::mat3 rightHandNormalMatrix = glm::transpose(glm::inverse(glm::mat3(rightHandMatrix)));
    rightHandMatrix = projectionMatrix * rightHandMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &rightHandModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &rightHandNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &rightHandMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //left upper arm stone-blue center(0, 6, -1)
    glm::mat4 leftArmMatrix = baseMatrix;
    if (robotStep == 1) {   //robot walk animation
        leftArmMatrix = glm::rotate(leftArmMatrix, 10.0f, glm::vec3(0.0, 0.0, 1.0));
    }
    if (robotStep == 2) {
        leftArmMatrix = glm::rotate(leftArmMatrix, -10.0f, glm::vec3(0.0, 0.0, 1.0));
    }
    leftArmMatrix = glm::rotate(leftArmMatrix, liftLeftArmAngle, glm::vec3(0.0, 0.0, -1.0));
    leftArmMatrix = glm::translate(leftArmMatrix, glm::vec3(0.0, -0.75, -1.0));
    glm::mat4 leftUpperArmMatrix = glm::scale(leftArmMatrix, glm::vec3(0.25, 0.75, 0.25));
    leftUpperArmMatrix = viewMatrix * leftUpperArmMatrix;
    glm::mat4 leftUpperArmModelviewMatrix = leftUpperArmMatrix;
    glm::mat3 leftUpperArmNormalMatrix = glm::transpose(glm::inverse(glm::mat3(leftUpperArmMatrix)));
    leftUpperArmMatrix = projectionMatrix * leftUpperArmMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leftUpperArmModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leftUpperArmNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leftUpperArmMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //left lower arm stone-blue center(0, 4.5, -1)
    leftArmMatrix = glm::translate(leftArmMatrix, glm::vec3(0.0, -1.5, 0.0));
    glm::mat4 leftLowerArmMatrix = glm::scale(leftArmMatrix, glm::vec3(0.25, 0.75, 0.25));
    leftLowerArmMatrix = viewMatrix * leftLowerArmMatrix;
    glm::mat4 leftLowerArmModelviewMatrix = leftLowerArmMatrix;
    glm::mat3 leftLowerArmNormalMatrix = glm::transpose(glm::inverse(glm::mat3(leftLowerArmMatrix)));
    leftLowerArmMatrix = projectionMatrix * leftLowerArmMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leftLowerArmModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leftLowerArmNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leftLowerArmMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //left hand stone-blue center(0, 3.5, -1)
    leftArmMatrix = glm::translate(leftArmMatrix, glm::vec3(0.0, -1.0, 0.0));
    glm::mat4 leftHandMatrix = glm::scale(leftArmMatrix, glm::vec3(0.5, 0.25, 0.25));
    leftHandMatrix = viewMatrix * leftHandMatrix;
    glm::mat4 leftHandModelviewMatrix = leftHandMatrix;
    glm::mat3 leftHandNormalMatrix = glm::transpose(glm::inverse(glm::mat3(leftHandMatrix)));
    leftHandMatrix = projectionMatrix * leftHandMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leftHandModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leftHandNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leftHandMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    //finish left arm
    
    //neck red center(0, 6.75, 0)
    baseMatrix = glm::translate(baseMatrix, glm::vec3(0.0, 0.5, 0.0));
    glm::mat4 neckMatrix = glm::scale(baseMatrix, glm::vec3(0.25, 0.25, 0.25));
    neckMatrix = viewMatrix * neckMatrix;
    glm::mat4 neckModelviewMatrix = neckMatrix;
    glm::mat3 neckNormalMatrix = glm::transpose(glm::inverse(glm::mat3(neckMatrix)));
    neckMatrix = projectionMatrix * neckMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &neckModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &neckNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &neckMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //head pink center(0, 7.5, 0)
    baseMatrix = glm::translate(baseMatrix, glm::vec3(0.0, 0.75, 0.0));
    glm::mat4 headMatrix = glm::scale(baseMatrix, glm::vec3(0.5, 0.5, 0.5));
    headMatrix = viewMatrix * headMatrix;
    glm::mat4 headModelviewMatrix = headMatrix;
    glm::mat3 headNormalMatrix = glm::transpose(glm::inverse(glm::mat3(headMatrix)));
    headMatrix = projectionMatrix * headMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &headModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &headNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &headMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+64);
    
    //face yellow center(0.6, 7.5, 0)
    glm::mat4 faceMatrix = glm::translate(baseMatrix, glm::vec3(0.4, 0.1, 0.1));
    glm::mat4 leftEyeMatrix = glm::scale(faceMatrix, glm::vec3(0.1, 0.1, 0.1));
    leftEyeMatrix = viewMatrix * leftEyeMatrix;
    glm::mat4 leftEyeModelviewMatrix = leftEyeMatrix;
    glm::mat3 leftEyeNormalMatrix = glm::transpose(glm::inverse(glm::mat3(leftEyeMatrix)));
    leftEyeMatrix = projectionMatrix * leftEyeMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leftEyeModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leftEyeNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leftEyeMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    faceMatrix = glm::translate(faceMatrix, glm::vec3(0.0, 0.0, -0.2));
    glm::mat4 rightEyeMatrix = glm::scale(faceMatrix, glm::vec3(0.1, 0.1, 0.1));
    rightEyeMatrix = viewMatrix * rightEyeMatrix;
    glm::mat4 rightEyeModelviewMatrix = rightEyeMatrix;
    glm::mat3 rightEyeNormalMatrix = glm::transpose(glm::inverse(glm::mat3(rightEyeMatrix)));
    rightEyeMatrix = projectionMatrix * rightEyeMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &rightEyeModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &rightEyeNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &rightEyeMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    faceMatrix = glm::translate(faceMatrix, glm::vec3(0.0, -0.3, 0.1));
    glm::mat4 mouseMatrix = glm::scale(faceMatrix, glm::vec3(0.1, 0.1, 0.2));
    mouseMatrix = viewMatrix * mouseMatrix;
    glm::mat4 mouseModelviewMatrix = mouseMatrix;
    glm::mat3 mouseNormalMatrix = glm::transpose(glm::inverse(glm::mat3(mouseMatrix)));
    mouseMatrix = projectionMatrix * mouseMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &mouseModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &mouseNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &mouseMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    
    //waist yellow center(0, 3.75, 0)
    baseMatrix = glm::translate(baseMatrix, glm::vec3(0.0, -3.75, 0.0));
    glm::mat4 waistMatrix = glm::scale(baseMatrix, glm::vec3(0.25, 0.25, 0.75));
    waistMatrix = viewMatrix * waistMatrix;
    glm::mat4 waistModelviewMatrix = waistMatrix;
    glm::mat3 waistNormalMatrix = glm::transpose(glm::inverse(glm::mat3(waistMatrix)));
    waistMatrix = projectionMatrix * waistMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &waistModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &waistNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &waistMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //start right leg, base on waist
    //right upper leg blue center(0, 2.75, 0.5)
    glm::mat4 rightLegMatrix = baseMatrix;
    if (robotStep == 1) {   //robot walk animation
        rightLegMatrix = glm::rotate(rightLegMatrix, 20.0f, glm::vec3(0.0, 0.0, 1.0));
    }
    if (robotStep == 2) {
        rightLegMatrix = glm::rotate(rightLegMatrix, -10.0f, glm::vec3(0.0, 0.0, 1.0));
    }
    rightLegMatrix = glm::translate(rightLegMatrix, glm::vec3(0.0, -1.0, 0.5));
    glm::mat4 rightUpperLegMatrix = glm::scale(rightLegMatrix, glm::vec3(0.25, 0.75, 0.25));
    rightUpperLegMatrix = viewMatrix * rightUpperLegMatrix;
    glm::mat4 rightUpperLegModelviewMatrix = rightUpperLegMatrix;
    glm::mat3 rightUpperLegNormalMatrix = glm::transpose(glm::inverse(glm::mat3(rightUpperLegMatrix)));
    rightUpperLegMatrix = projectionMatrix * rightUpperLegMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &rightUpperLegModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &rightUpperLegNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &rightUpperLegMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //right lower leg blue center(0, 1.25, 0.5)
    if (robotStep == 1) {
        rightLegMatrix = glm::translate(rightLegMatrix, glm::vec3(0.0, -0.75, 0.0));
        rightLegMatrix = glm::rotate(rightLegMatrix, -10.0f, glm::vec3(0.0, 0.0, 1.0));
        rightLegMatrix = glm::translate(rightLegMatrix, glm::vec3(0.0, -0.75, 0.0));
    }
    if (robotStep == 2) {
        rightLegMatrix = glm::translate(rightLegMatrix, glm::vec3(0.0, -0.75, 0.0));
        rightLegMatrix = glm::rotate(rightLegMatrix, 5.0f, glm::vec3(0.0, 0.0, 1.0));
        rightLegMatrix = glm::translate(rightLegMatrix, glm::vec3(0.0, -0.75, 0.0));
    }
    if (robotStep == 0) {
        rightLegMatrix = glm::translate(rightLegMatrix, glm::vec3(0.0, -1.5, 0.0));
    }
    glm::mat4 rightLowerLegMatrix = glm::scale(rightLegMatrix, glm::vec3(0.25, 0.75, 0.25));
    rightLowerLegMatrix = viewMatrix * rightLowerLegMatrix;
    glm::mat4 rightLowerLegModelviewMatrix = rightLowerLegMatrix;
    glm::mat3 rightLowerLegNormalMatrix = glm::transpose(glm::inverse(glm::mat3(rightLowerLegMatrix)));
    rightLowerLegMatrix = projectionMatrix * rightLowerLegMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &rightLowerLegModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &rightLowerLegNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &rightLowerLegMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //right foot blue center(0, 0.25, 0.5)
    rightLegMatrix = glm::translate(rightLegMatrix, glm::vec3(0.0, -1.0, 0.0));
    glm::mat4 rightFootMatrix = glm::scale(rightLegMatrix, glm::vec3(0.5, 0.25, 0.25));
    rightFootMatrix = viewMatrix * rightFootMatrix;
    glm::mat4 rightFootModelviewMatrix = rightFootMatrix;
    glm::mat3 rightFootNormalMatrix = glm::transpose(glm::inverse(glm::mat3(rightFootMatrix)));
    rightFootMatrix = projectionMatrix * rightFootMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &rightFootModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &rightFootNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &rightFootMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //start left leg, base on waist
    //left upper leg blue center(0, 2.75, -0.5)
    glm::mat4 leftLegMatrix = baseMatrix;
    if (robotStep == 1) {   //robot walk animation
        leftLegMatrix = glm::rotate(leftLegMatrix, -10.0f, glm::vec3(0.0, 0.0, 1.0));
    }
    if (robotStep == 2) {
        leftLegMatrix = glm::rotate(leftLegMatrix, 20.0f, glm::vec3(0.0, 0.0, 1.0));
    }
    leftLegMatrix = glm::translate(leftLegMatrix, glm::vec3(0.0, -1.0, -0.5));
    glm::mat4 leftUpperLegMatrix = glm::scale(leftLegMatrix, glm::vec3(0.25, 0.75, 0.25));
    leftUpperLegMatrix = viewMatrix * leftUpperLegMatrix;
    glm::mat4 leftUpperLegModelviewMatrix = leftUpperLegMatrix;
    glm::mat3 leftUpperLegNormalMatrix = glm::transpose(glm::inverse(glm::mat3(leftUpperLegMatrix)));
    leftUpperLegMatrix = projectionMatrix * leftUpperLegMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leftUpperLegModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leftUpperLegNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leftUpperLegMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //left lower leg blue center(0, 1.25, -0.5)
    if (robotStep == 1) {
        leftLegMatrix = glm::translate(leftLegMatrix, glm::vec3(0.0, -0.75, 0.0));
        leftLegMatrix = glm::rotate(leftLegMatrix, 5.0f, glm::vec3(0.0, 0.0, 1.0));
        leftLegMatrix = glm::translate(leftLegMatrix, glm::vec3(0.0, -0.75, 0.0));
    }
    if (robotStep == 2) {
        leftLegMatrix = glm::translate(leftLegMatrix, glm::vec3(0.0, -0.75, 0.0));
        leftLegMatrix = glm::rotate(leftLegMatrix, -10.0f, glm::vec3(0.0, 0.0, 1.0));
        leftLegMatrix = glm::translate(leftLegMatrix, glm::vec3(0.0, -0.75, 0.0));
    }
    if (robotStep == 0) {
        leftLegMatrix = glm::translate(leftLegMatrix, glm::vec3(0.0, -1.5, 0.0));
    }
    glm::mat4 leftLowerLegMatrix = glm::scale(leftLegMatrix, glm::vec3(0.25, 0.75, 0.25));
    leftLowerLegMatrix = viewMatrix * leftLowerLegMatrix;
    glm::mat4 leftLowerLegModelviewMatrix = leftLowerLegMatrix;
    glm::mat3 leftLowerLegNormalMatrix = glm::transpose(glm::inverse(glm::mat3(leftLowerLegMatrix)));
    leftLowerLegMatrix = projectionMatrix * leftLowerLegMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leftLowerLegModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leftLowerLegNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leftLowerLegMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
    //left foot blue center(0, 0.25, -0.5)
    leftLegMatrix = glm::translate(leftLegMatrix, glm::vec3(0.0, -1.0, 0.0));
    glm::mat4 leftFootMatrix = glm::scale(leftLegMatrix, glm::vec3(0.5, 0.25, 0.25));
    leftFootMatrix = viewMatrix * leftFootMatrix;
    glm::mat4 leftFootModelviewMatrix = leftFootMatrix;
    glm::mat3 leftFootNormalMatrix = glm::transpose(glm::inverse(glm::mat3(leftFootMatrix)));
    leftFootMatrix = projectionMatrix * leftFootMatrix;
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &leftFootModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &leftFootNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &leftFootMatrix[0][0]);
    glDrawElements(GL_TRIANGLE_STRIP, sphereIndexTotal, GL_UNSIGNED_INT, (char*) NULL + 0);
    
}
void displayBunny() {
    glBindBuffer(GL_ARRAY_BUFFER, vboHandle[4]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO[4]);
    
    glVertexAttribPointer(positionShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+0);
    glVertexAttribPointer(normalShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+16);
    glVertexAttribPointer(colorShader, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+32);
    glVertexAttribPointer(textureCoordShader, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (char*) NULL+ 80);
    
    glm::mat4 bunnyMatrix = glm::mat4(1.0f);
    
    bunnyMatrix = glm::rotate(bunnyMatrix, x_angle, glm::vec3(0.0, 1.0, 0.0));
    bunnyMatrix = glm::rotate(bunnyMatrix, y_angle, glm::vec3(1.0, 0.0, 0.0));
    bunnyMatrix = glm::scale(bunnyMatrix, glm::vec3(scale_size, scale_size, scale_size));
    
    bunnyMatrix = glm::translate(bunnyMatrix, glm::vec3(8.0, -0.5, -2.0));
    bunnyMatrix = glm::rotate(bunnyMatrix, 90.0f, glm::vec3(0.0, 1.0, 0.0));
    bunnyMatrix = glm::scale(bunnyMatrix, glm::vec3(3.0, 3.0, 3.0));
    bunnyMatrix = viewMatrix * bunnyMatrix;
    glm::mat4 bunnyModelviewMatrix = bunnyMatrix;
    glm::mat3 bunnyNormalMatrix = glm::transpose(glm::inverse(glm::mat3(bunnyMatrix)));
    bunnyMatrix = projectionMatrix * bunnyMatrix;
    
    glUniformMatrix4fv(modelviewMatrixShader, 1, GL_FALSE, &bunnyModelviewMatrix[0][0]);
    glUniformMatrix3fv(normalMatrixShader, 1, GL_FALSE, &bunnyNormalMatrix[0][0]);
    glUniformMatrix4fv(local2clipShader, 1, GL_FALSE, &bunnyMatrix[0][0]);
    glDrawElements(GL_TRIANGLES, bunnyIndexTotal, GL_UNSIGNED_INT, (char*) NULL+0);
}
void display() {
    if (flagWireMode == true) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
    else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    
    if (flagEdgeDetection == false) {
        
        glUseProgram(programObject);

        glBindFramebuffer(GL_FRAMEBUFFER, shadowFBO);
        
        glPointSize(3.0f);
//        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_POINT_SPRITE);
        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
        glViewport(0, 0, 1024, 1024);
        glEnable(GL_DEPTH_TEST);    // need depth test to correctly draw 3D objects
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glUniform1i(textureShader, 0); // associate Tex1 to texture unit 0
        glUniform1i(cubeMapTextureShader, 2);
        glUniform1i(shadowMapShader, 5);
        glUniform1i(flagTextureShader, false);
        glUniform1i(flagCubemapShader, false);
        glUniform1i(flagShadowMapShader, false);
        
        glUniform1i(flagNormalMapShader, false);
        glUniform1i(normalMapShader, 1);
        glUniform1i(flagCartoonShader, cartoon);
        glUniform1f(timeShader, ntime);
        glUniform1i(flagEdgeShader, false);
        glUniform1i(flagAnimationShader, false);
        
        
        viewMatrix = viewMatrixLight;
        
    //    if (trackballMove)
    //	{
    //        viewMatrix = viewMatrix * trackball;
    //	}
        
        if (flagScene == true) {
            displayScene();
    //        displayWalls();
        }
        displayRobot();
        displayBunny();
        
        
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        
        glViewport(0, 0, 600, 600);
        glEnable(GL_DEPTH_TEST);    // need depth test to correctly draw 3D objects
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glUniform1i(textureShader, 0); // associate Tex1 to texture unit 0
        glUniform1i(cubeMapTextureShader, 2);
        glUniform1i(shadowMapShader, 5);
        glUniform1i(flagTextureShader, false);
        glUniform1i(flagCubemapShader, false);
        
    //    viewMatrix = viewMatrixCamera;
        viewMatrix = viewMatrixCamera * trackball;

    //    if (recordView == true) {
    //        viewMatrixCamera = viewMatrix;
    //        recordView = false;
    //    }
    //    
        if (flagScene == true) {
            displayScene1();
            displayWalls1();
        }
        displayRobot();
        displayBunny();
        
        glutSwapBuffers();
        
        glUseProgram(0);
    }
    else {
        glUseProgram(programObject);
        
        glBindFramebuffer(GL_FRAMEBUFFER, shadowFBO);
        
        glViewport(0, 0, 1024, 1024);
        glEnable(GL_DEPTH_TEST);    // need depth test to correctly draw 3D objects
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glUniform1i(textureShader, 0); // associate Tex1 to texture unit 0
        glUniform1i(cubeMapTextureShader, 2);
        glUniform1i(shadowMapShader, 5);
        glUniform1i(flagTextureShader, false);
        glUniform1i(flagCubemapShader, false);
        glUniform1i(flagShadowMapShader, false);
        
        glUniform1i(flagNormalMapShader, false);
        glUniform1i(normalMapShader, 1);
        glUniform1i(flagCartoonShader, cartoon);
        glUniform1f(timeShader, ntime);
        ntime += 1.0;
        glUniform1i(flagEdgeShader, true);
        viewMatrix = viewMatrixCamera;
        
        if (flagScene == true) {
            displayScene();
        }
        displayRobot();
        displayBunny();
        
        
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        
        glViewport(0, 0, 600, 600);
        glEnable(GL_DEPTH_TEST);    // need depth test to correctly draw 3D objects
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glUniform1i(textureShader, 0); // associate Tex1 to texture unit 0
        glUniform1i(cubeMapTextureShader, 2);
        glUniform1i(shadowMapShader, 5);
        glUniform1i(flagTextureShader, false);
        glUniform1i(flagCubemapShader, false);
                
        viewMatrix = viewMatrixCamera * trackball;
        
//        if (flagScene == true) {
//            displayScene1();
//            displayWalls1();
//        }
//        displayRobot();
//        displayBunny();
        glUniform1i(flagEdgeShader, true);
        displayEdge();
        glUniform1i(flagEdgeShader, true);
        glutSwapBuffers();
        
        glUseProgram(0);
    }
}

void store_ply(PlyFile* input, Vertex ***vertices, Face ***faces,
               unsigned int* vertexcount, unsigned int* facecount,
               int* vertexnormals, int* facenormals) {
    int i, j;
    
    // go through the element types
    for(i = 0; i < input->num_elem_types; i = i + 1) {
        int count;
        
        // setup the element for reading and get the element count
        char* element = setup_element_read_ply(input, i, &count);
        
        // get vertices
        if(strcmp("vertex", element) == 0) {
            *vertices = (Vertex**)malloc(sizeof(Vertex) * count);
            *vertexcount = count;
            
            // run through the properties and store them
            for(j = 0; j < input->elems[i]->nprops; j = j + 1) {
                PlyProperty* property = input->elems[i]->props[j];
                PlyProperty setup;
                
                if(strcmp("x", property->name) == 0 &&
                   property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[0];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Vertex, location[0]);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                }
                else if(strcmp("y", property->name) == 0 &&
                        property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[1];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Vertex, location[1]);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                }
                else if(strcmp("z", property->name) == 0 &&
                        property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[2];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Vertex, location[2]);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                }
                else if(strcmp("nx", property->name) == 0 &&
                        property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[3];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Vertex, normal[0]);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                    *vertexnormals = 1;
                }
                else if(strcmp("ny", property->name) == 0 &&
                        property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[4];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Vertex, normal[1]);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                    *vertexnormals = 1;
                }
                else if(strcmp("nz", property->name) == 0 &&
                        property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[5];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Vertex, normal[2]);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                    *vertexnormals = 1;
                }
                // dunno what it is
                else {
                    fprintf(stderr, "unknown property type found in %s: %s\n",
                            element, property->name);
                }
            }
            
            // do this if you want to grab the other data
            // list_pointer = get_other_properties_ply
            //                (input, offsetof(Vertex, struct_pointer));
            
            // copy the data
            for(j = 0; j < count; j = j + 1) {
                (*vertices)[j] = (Vertex*)malloc(sizeof(Vertex));
                
                get_element_ply(input, (void*)((*vertices)[j]));
            }
        }
        // get faces
        else if(strcmp("face", element) == 0) {
            *faces = (Face**)malloc(sizeof(Face) * count);
            *facecount = count;
            
            // run through the properties and store them
            for(j = 0; j < input->elems[i]->nprops; j = j + 1) {
                PlyProperty* property = input->elems[i]->props[j];
                PlyProperty setup;
                
                if(strcmp("vertex_indices", property->name) == 0 &&
                   property->is_list == PLY_LIST) {
                    
                    setup.name = string_list[6];
                    setup.internal_type = Uint32;
                    setup.offset = offsetof(Face, vertices);
                    setup.count_internal = Uint32;
                    setup.count_offset = offsetof(Face, count);
                    
                    setup_property_ply(input, &setup);
                }
                else if(strcmp("nx", property->name) == 0 &&
                        property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[3];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Face, nx);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                    *facenormals = 1;
                }
                else if(strcmp("ny", property->name) == 0 &&
                        property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[4];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Face, ny);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                    *facenormals = 1;
                }
                else if(strcmp("nz", property->name) == 0 &&
                        property->is_list == PLY_SCALAR) {
                    
                    setup.name = string_list[5];
                    setup.internal_type = Float32;
                    setup.offset = offsetof(Face, nz);
                    setup.count_internal = 0;
                    setup.count_offset = 0;
                    
                    setup_property_ply(input, &setup);
                    *facenormals = 1;
                }
                // dunno what it is
                else {
                    fprintf(stderr, "unknown property type found in %s: %s\n",
                            element, property->name);
                }
            }
            
            // do this if you want to grab the other data
            // list_pointer = get_other_properties_ply
            //                (input, offsetof(Face, struct_pointer));
            
            // copy the data
            for(j = 0; j < count; j = j + 1) {
                (*faces)[j] = (Face*)malloc(sizeof(Face));
                
                get_element_ply(input, (void*)((*faces)[j]));
            }
        }
        // who knows?
        else {
            fprintf(stderr, "unknown element type found: %s\n", element);
        }
    }
    
    // if you want to grab the other data do this
    // get_other_element_ply(input);
}

void find_center (float& cx, float& cy, float& cz,
                  float& minx, float& maxx, float&miny,
                  float &maxy, float &minz, float & maxz) {
    float x, y, z;
    float min_x=9999, max_x=-9999, min_y=9999, max_y=-9999;
    float min_z=9999, max_z=-9999;
    
    x = y= z = 0;
    for(int i = 0; i < vertexcount; i++) {
        x += vertices[i]->location[0];
        y += vertices[i]->location[1];
        z += vertices[i]->location[2];
        if (min_x >vertices[i]->location[0]) min_x = vertices[i]->location[0];
        if (max_x <vertices[i]->location[0]) max_x = vertices[i]->location[0];
        if (min_y >vertices[i]->location[1]) min_y = vertices[i]->location[1];
        if (max_y <vertices[i]->location[1]) max_y = vertices[i]->location[1];
        if (min_z >vertices[i]->location[2]) min_z = vertices[i]->location[2];
        if (max_z <vertices[i]->location[2]) max_z = vertices[i]->location[2];
    }
    cx = x / (float) vertexcount;
    cy = y / (float) vertexcount;
    cz = z / (float) vertexcount;
    minx = min_x; maxx = max_x;
    miny = min_y; maxy = max_y;
    minz = min_z; maxz = max_z;
}

void startMotion(int x, int y)
{

	startX = x;
	startY = y;
	curx = x;
	cury = y;

    lastPos.x = (float)x / winWidth * 2.0 - 1.0;
    lastPos.y = (float)y / winHeight * 2.0 - 1.0;

    if (lastPos.x * lastPos.x + lastPos.y * lastPos.y > 1.0) {
        lastPos = glm::normalize(lastPos);
    }
    
    lastPos.z = sqrt(1.0 - lastPos.x * lastPos. x - lastPos.y * lastPos.y);
    
    lastPos = glm::normalize(lastPos);
    
}

void stopMotion(int x, int y)
{
    recordView = true;
    viewMatrixCamera = viewMatrix;

//	/* check if position has changed */
//	if (startX != x || startY != y) {
//		redrawContinue = true;
//	}
//	else {
//		angle = 0.0;
//		redrawContinue = true;
//		trackballMove = true;
//	}
}

void mymouse(int button, int state, int x, int y) {
//    original mouse
//    if (state == GLUT_DOWN) {
//        press_x = x; press_y = y;
//        if (button == GLUT_LEFT_BUTTON)
//            xform_mode = XFORM_ROTATE;
//        else if (button == GLUT_RIGHT_BUTTON)
//            xform_mode = XFORM_SCALE;
//    }
//    else if (state == GLUT_UP) {
//        xform_mode = XFORM_NONE;
//    }
    
    if (state == GLUT_DOWN) {
        press_x = x; press_y = y;
        if (button == GLUT_RIGHT_BUTTON)
            xform_mode = XFORM_SCALE;
    }
    else if (state == GLUT_UP) {
        xform_mode = XFORM_NONE;
    }
    
    /**** trackball ****/
//    if(button==GLUT_RIGHT_BUTTON) {
//		exit(0);
//	}
	/* holding down left button allows user to rotate cube */
	if(button==GLUT_LEFT_BUTTON) {
        y = winHeight - y;

		switch(state)
		{
			case GLUT_DOWN:
                printf("mouse donw: (%d, %d);\n", x, y);
                startMotion(x, y);
                break;
			case GLUT_UP:
                stopMotion(x, y);
                printf("mouse up: (%d, %d);\n", x, y);
                printf("angle: %f, vec: (%f, %f, %f)\n", angle, axis.x, axis.y, axis.z);
                break; 
		}
	}
}
void mymotion(int x, int y) {
//    original motion
//    if (xform_mode==XFORM_ROTATE) {
//        x_angle += (x - press_x)/5.0;
//        if (x_angle > 180) x_angle -= 360;
//        else if (x_angle <-180) x_angle += 360;
//        press_x = x;
//        
//        y_angle += (y - press_y)/5.0;
//        if (y_angle > 180) y_angle -= 360;
//        else if (y_angle <-180) y_angle += 360;
//        press_y = y;
//    }
//	else if (xform_mode == XFORM_SCALE){
//        float old_size = scale_size;
//        scale_size *= (1+ (y - press_y)/60.0);
//        if (scale_size <0) scale_size = old_size;
//        press_y = y;
//    }
	if (xform_mode == XFORM_SCALE){
        float old_size = scale_size;
        scale_size *= (1+ (y - press_y)/60.0);
        if (scale_size <0) scale_size = old_size;
        press_y = y;
    }
    /**** trackball ****/
    

    
        y = winHeight - y;

        glm::vec3 curPos = glm::vec3(0.0);
        //	float dx, dy, dz;
        /* compute position on hemisphere */
        
        curPos.x = (float)x / winWidth * 2.0 - 1.0;
        curPos.y = (float)y / winHeight * 2.0 - 1.0;
        
//        cout << "x, y: " << x << ", " << y << endl;
//        cout << "curPos (before): " << curPos.x << ", " << curPos.y << ", " << curPos.z << endl;
    
        if (curPos.x * curPos.x + curPos.y * curPos.y > 1.0) {
            curPos = glm::normalize(curPos);
        }
        
        curPos.z = sqrt(1.0 - curPos.x * curPos. x - curPos.y * curPos.y);
        
        curPos = glm::normalize(curPos);
        
//        cout << "curPos (after): " << curPos.x << ", " << curPos.y << ", " << curPos.z << endl;
//        cout << "lastPos: " << lastPos.x << lastPos.y << lastPos.z << endl << endl;
        axis = glm::cross(lastPos, curPos);
        axis = glm::normalize(axis);
        angle = glm::acos(glm::dot(lastPos, curPos));
//        cout << "angle (before): " << angle << endl;
        angle = 360.0 * angle / (2.0 * M_PI);
//        cout << "angle (after): " << angle << endl;
//        lastPos = curPos;
//        cout << "start tracking mouse: " << angle << axis.x << axis.y << axis.z << endl;
        trackball = glm::rotate(glm::mat4(1.0), angle, axis);

//	if(trackingMouse)
//	{
//		/* compute the change in position on the hemisphere */
//		dx = curPos[0] - lastPos[0];
//		dy = curPos[1] - lastPos[1];
//		dz = curPos[2] - lastPos[2];
//		if (dx || dy || dz)
//		{
//			/* compute theta and cross product */
//			angle = 90.0 * sqrt(dx*dx + dy*dy + dz*dz);
//			axis[0] = lastPos[1]*curPos[2] - lastPos[2]*curPos[1];
//			axis[1] = lastPos[2]*curPos[0] - lastPos[0]*curPos[2];
//			axis[2] = lastPos[0]*curPos[1] - lastPos[1]*curPos[0];
//			/* update position */
//			lastPos[0] = curPos[0];
//			lastPos[1] = curPos[1];
//			lastPos[2] = curPos[2];
//		}
//        
//	}
    

    glutPostRedisplay();
}
void mykey(unsigned char key, int x, int y) {
	if (key == 'q' || key == 'Q') {
        exit(1);
    }
    if (key == 's' || key == 'S') {
        if (flagScene == false) {
            flagScene = true;
        }
        else {
            flagScene = false;
        }
    }
    if (key == 'w' || key == 'W') {
        if (flagWireMode == false) {
            flagWireMode = true;
        }
        else {
            flagWireMode = false;
        }
    }
    if (key == 'f' || key == 'F') {
        transformMatrix = glm::translate(transformMatrix, glm::vec3(0.5, 0.0, 0.0));
        if (robotStep != 1) {
            robotStep = 1;
        }
        else {
            robotStep = 2;
        }
    }
    if (key == 'b' || key == 'B') {
        transformMatrix = glm::translate(transformMatrix, glm::vec3(-0.5, 0.0, 0.0));
        if (robotStep != 1) {
            robotStep = 1;
        }
        else {
            robotStep = 2;
        }
    }
    if (key == 't') {
        transformMatrix = glm::rotate(transformMatrix, 15.0f, glm::vec3(0.0, 1.0, 0.0));
    }
    if (key == 'T') {
        transformMatrix = glm::rotate(transformMatrix, -15.0f, glm::vec3(0.0, 1.0, 0.0));
    }
    if (key == 'r' || key == 'R') {
        liftRightArmAngle += 15;
    }
    if (key == 'l' || key == 'L') {
        liftLeftArmAngle -= 15;
    }
    if (key == 'c' || key == 'C') {
        if (cartoon == true) {
            cartoon = false;
        }
        else {
            cartoon = true;
        }
    }
    if (key == 'e' || key == 'E') {
        if (flagEdgeDetection == true) {
            flagEdgeDetection = false;
        }
        else {
            flagEdgeDetection = true;
        }
    }
    if (key == 'a' || key == 'A') {
        if (flagAnimation == true) {
            flagAnimation = false;
        }
        else {
            flagAnimation = true;
            ntime = 0;
        }
    }
	glutPostRedisplay();
}
void myidle(int val) {
    int ani = 0;
    if (ani == 0) glutTimerFunc(15, myidle, 0);
    ntime += 0.05;
    if (ntime > 10) {
        ntime = 0;
    }
    glutPostRedisplay();
}
void myreadimage(char *filename, int i) {
    //    const char* filename = "transformers.png";
    
    //load and decode
    std::vector<unsigned char> image;
    unsigned width, height;
    unsigned error = lodepng::decode(image, width, height, filename);
    
    //if there's an error, display it
    if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    
    //    cout << "image Size: (" << width << ", " << height << ")" << endl;
    //    cout << "image total: " << image.size() << endl;
    
    switch (i) {
        case 0:     //desktop surface
            imageWidth = width;
            imageHeight = height;
            
            for (int i = 0; i < image.size(); i ++) {
                imageData[i] = image[i];
            }
            break;
        case 1:     //px
            imageWidth1 = width;
            imageHeight1 = height;
            
            for (int i = 0; i < image.size(); i ++) {
                imageData1[i] = image[i];
            }
            break;
        case 2:     //py
            imageWidth2 = width;
            imageHeight2 = height;
            
            for (int i = 0; i < image.size(); i ++) {
                imageData2[i] = image[i];
            }
            break;
        case 3:     //pz
            imageWidth3 = width;
            imageHeight3 = height;
            
            for (int i = 0; i < image.size(); i ++) {
                imageData3[i] = image[i];
            }
            break;
        case 4:     //nx
            imageWidth4 = width;
            imageHeight4 = height;
            
            for (int i = 0; i < image.size(); i ++) {
                imageData4[i] = image[i];
            }
            break;
        case 5:     //ny
            imageWidth5 = width;
            imageHeight5 = height;
            
            for (int i = 0; i < image.size(); i ++) {
                imageData5[i] = image[i];
            }
            break;
        case 6:     //nz
            imageWidth6 = width;
            imageHeight6 = height;
            
            for (int i = 0; i < image.size(); i ++) {
                imageData6[i] = image[i];
            }
            break;
        case 7:
            imageWidthNormalMap = width;
            imageHeightNormalMap = height;
            
            for (int i = 0; i < image.size(); i ++) {
                imageDataNormalMap[i] = image[i];
            }
            
            
        default:
            break;
    }
    
}

int main(int argc, char** argv) {
    
    cout << "hello! " << endl;
//    cout << "Hello! " << argv[1] << endl;
    
    PlyFile* input;
    
    // get the ply structure and open the file
    input = read_ply(stdin);
    
    // read in the data
    store_ply(input,
              &vertices, &faces,
              &vertexcount, &facecount,
              &vertexnormals, &facenormals);
    
    // close the file
    close_ply(input);
    find_center(cx, cy, cz, x_min, x_max,
                y_min, y_max, z_min, z_max);
    printf("geometry center = %f %f %f \n", cx, cy, cz);
    printf("geometry bound = x: %f %f y: %f %f z: %f %f\n",
           x_min, x_max, y_min, y_max, z_min, z_max);
    
    myreadimage((char *)"landscape.png", 0);
    myreadimage((char *)"px.png", 1);
    myreadimage((char *)"py.png", 2);
    myreadimage((char *)"pz.png", 3);
    myreadimage((char *)"nx.png", 4);
    myreadimage((char *)"ny.png", 5);
    myreadimage((char *)"nz.png", 6);
    
    myreadimage((char *)"normalMap.png", 7);
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH);
    glutInitWindowSize(600,600);
    
    glutCreateWindow("CSE 5542 lab 5 - Yuxiang Kou");
    glutDisplayFunc(display);
    glutMouseFunc(mymouse);
    glutMotionFunc(mymotion);
    glutKeyboardFunc(mykey);
    
    // initialize GLEW
    GLenum err = glewInit();
    if ( err != GLEW_OK)  printf(" Error initializing GLEW! \n");
    else printf("Initializing GLEW succeeded!\n"); 
    
    programObject = SetupGLSL((char *)"lab5Shader");  //create shaders - assume the shaders are lab3Shader.vert and lab3Shader.frag        

    initShader();
    initVBO();
    initTexture();
    initFBO();
    
    
    glutTimerFunc(15, myidle, 0);

    glutMainLoop();
    
}
