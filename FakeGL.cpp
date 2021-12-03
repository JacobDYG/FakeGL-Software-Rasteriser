//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  ------------------------
//  FakeGL.cpp
//  ------------------------
//  
//  A unit for implementing OpenGL workalike calls
//  
///////////////////////////////////////////////////

#include "FakeGL.h"
#include <math.h>

//-------------------------------------------------//
//                                                 //
// CONSTRUCTOR / DESTRUCTOR                        //
//                                                 //
//-------------------------------------------------//

// constructor
FakeGL::FakeGL()
// Initialisation list
    :
    // Buffers
    depthVal(0.0f, 0.0f, 0.0f, 0.0f),
    clearColorVal(0.0f, 0.0f, 0.0f, 0.0f),
    // Default colour should be white
    currentColor(255.0f, 255.0f, 255.0f, 255.0f),
    // Initialise current normal to something safe
    currentNormal(0.0f, 0.0f, 0.0f, 1.0f),
    // Lighting values
    lightPos(0.0f, 0.0f, 0.0f, 1.0f),
    ambientColor(1.0f, 1.0f, 1.0f, 1.0f),
    diffuseColor(1.0f, 1.0f, 1.0f, 1.0f),
    specularColor(1.0f, 1.0f, 1.0f, 1.0f),
    // Default material properties, from OpenGL spec
    // Reflectance
    currentAmbientReflectance(0.2f, 0.2f, 0.2f, 1.0f),
    currentDiffuseReflectance(0.8f, 0.8f, 0.8f, 1.0f),
    currentSpecularReflectance(0.0f, 0.0f, 0.0f, 1.0f),
    // Additives
    currentEmission(0.0f, 0.0f, 0.0f, 1.0f),
    currentShininess(0.0f)
{ // constructor
    // Set default matrix values
    modelViewMat.SetIdentity();
    projectionMat.SetIdentity();
} // constructor

// destructor
FakeGL::~FakeGL()
    { // destructor
    } // destructor

//-------------------------------------------------//
//                                                 //
// GEOMETRIC PRIMITIVE ROUTINES                    //
//                                                 //
//-------------------------------------------------//

// starts a sequence of geometric primitives
void FakeGL::Begin(unsigned int PrimitiveType)
{ // Begin()
    // Check the previous sequence has been closed with FakeGL::End()
    if (!primitiveAssembly)
    {
        // Check the primitive is valid
        if (PrimitiveType == FAKEGL_POINTS || PrimitiveType == FAKEGL_LINES || PrimitiveType == FAKEGL_TRIANGLES)
        {
            // Primitive is valid, set it
            primitiveMode = PrimitiveType;
            // Indicate that a sequence of vertices is being specified
            primitiveAssembly = true;
        }
        else
        {
            // The requested primitive mode is wrong. Set it to -1 as in (some) cases this will make it more obvious something is wrong
            primitiveMode = -1;
        }
    }
} // Begin()

// ends a sequence of geometric primitives
void FakeGL::End()
{ // End()
    // Calling later pipeline stages
    while (RasterisePrimitive())
    {
        while (fragmentQueue.size() > 0)
        {
            ProcessFragment();
        }
    }
    // Set the primitive value to -1
    primitiveMode = -1;
    // Indicate that primitive specification is done
    primitiveAssembly = false;
} // End()

// sets the size of a point for drawing
void FakeGL::PointSize(float size)
{ // PointSize()
    // Set the point size if it is larger than 0
    if (size > 0.0)
    {
        pointSize = size;
    }
} // PointSize()

// sets the width of a line for drawing purposes
void FakeGL::LineWidth(float width)
{ // LineWidth()
    // Set the line width if it is larger than 0
    if (width > 0.0)
    {
        lineWidth = width;
    }
} // LineWidth()

//-------------------------------------------------//
//                                                 //
// MATRIX MANIPULATION ROUTINES                    //
//                                                 //
//-------------------------------------------------//

// set the matrix mode (i.e. which one we change)   
void FakeGL::MatrixMode(unsigned int whichMatrix)
{ // MatrixMode()
    // Check this is a valid matrix mode
    if (whichMatrix == FAKEGL_MODELVIEW || whichMatrix == FAKEGL_PROJECTION)
    {
        // Set the matrix mode
        matrixMode = whichMatrix;
    }
} // MatrixMode()

// pushes a matrix on the stack
void FakeGL::PushMatrix()
{ // PushMatrix()
    // Check the matrix mode to determine which to push & which stack to change
    if (matrixMode == FAKEGL_MODELVIEW)
    {
        modelViewMatStack.push(modelViewMat);
    }
    else if (matrixMode == FAKEGL_PROJECTION)
    {
        projectionMatStack.push(projectionMat);
    }
} // PushMatrix()

// pops a matrix off the stack
void FakeGL::PopMatrix()
    { // PopMatrix()
    // Check the matrix mode to determine which stack to pop from & into which matrix
    if (matrixMode == FAKEGL_MODELVIEW)
    {
        // Read the matrix at the top of the stack, then pop it
        modelViewMat = modelViewMatStack.top();
        modelViewMatStack.pop();
    }
    else if (matrixMode == FAKEGL_PROJECTION)
    {
        projectionMat = projectionMatStack.top();
        projectionMatStack.pop();
    }
    } // PopMatrix()

// load the identity matrix
void FakeGL::LoadIdentity()
{ // LoadIdentity()
    // Check the matrix mode to determine which to modify
    if (matrixMode == FAKEGL_MODELVIEW)
    {
        // Load the identity matrix
        modelViewMat.SetIdentity();
    }
    else if (matrixMode == FAKEGL_PROJECTION)
    {
        projectionMat.SetIdentity();
    }
} // LoadIdentity()

// multiply by a known matrix in column-major format
void FakeGL::MultMatrixf(const float *columnMajorCoordinates)
{ // MultMatrixf()
    // Initialise a new matrix, store specified matrix as row major
    Matrix4 rowMajor;
    for (size_t col = 0; col < 4; col++)
    {
        for (size_t row = 0; row < 4; row++)
        {
            rowMajor.coordinates[row][col] = columnMajorCoordinates[(col * 4) + row];
        }
    }
    // Check the matrix mode
    if (matrixMode == FAKEGL_MODELVIEW)
    {
        // Multiply the current matrix with the new one
        modelViewMat = modelViewMat * rowMajor;
    }
    else if (matrixMode == FAKEGL_PROJECTION)
    {
        projectionMat = projectionMat * rowMajor;
    }
} // MultMatrixf()

// sets up a perspective projection matrix
void FakeGL::Frustum(float left, float right, float bottom, float top, float zNear, float zFar)
{ // Frustum()
    // Describes a perspective matrix
    // Applied to CURRENT matrix, though generally used for projection
    Matrix4 frustumMat;
    frustumMat.SetZero();
    // Set components
    frustumMat[0][0] = (2 * zNear) / (right - left);
    frustumMat[1][1] = (2 * zNear) / (top - bottom);
    frustumMat[2][0] = -(right + left) / (right - left);
    frustumMat[2][1] = (right + left) / (right - left);
    frustumMat[2][2] = (top + bottom) / (top - bottom);
    frustumMat[2][3] = -1;
    frustumMat[3][2] = -(2 * zFar * zNear) / (zFar - zNear);
    // Multiply with current matrix
    MultMatrixf(frustumMat.columnMajor().coordinates);
    // Update near and far planes
    FakeGL::zNear = zNear;
    FakeGL::zFar = zFar;
} // Frustum()

// sets an orthographic projection matrix
void FakeGL::Ortho(float left, float right, float bottom, float top, float zNear, float zFar)
{ // Ortho()
    // Describes an orthographic matrix (parallel projection)
    // Applied to CURRENT matrix, though generally used for projection
    Matrix4 orthoMat;
    orthoMat.SetIdentity();
    // Set components
    orthoMat[0][0] = 2 / (right - left);
    orthoMat[1][1] = 2 / (top - bottom);
    orthoMat[2][2] = -2 / (zFar - zNear);
    orthoMat[3][0] = -(right + left) / (right - left);
    orthoMat[3][1] = -(top + bottom) / (top - bottom);
    orthoMat[3][2] = -(zFar + zNear) / (zFar - zNear);
    // Multiply with current matrix
    MultMatrixf(orthoMat.columnMajor().coordinates);
    // Update near and far planes
    FakeGL::zNear = zNear;
    FakeGL::zFar = zFar;
} // Ortho()

// rotate the matrix
void FakeGL::Rotatef(float angle, float axisX, float axisY, float axisZ)
{ // Rotatef()
    // Produces a rotation of angle degrees around the axes x, y, z
    // Applied to CURRENT matrix, though generally used for modelview
    Matrix4 rotationMatrix;
    rotationMatrix.SetRotation(Cartesian3(axisX, axisY, axisZ), angle);
    // Multiply by the current matrix
    MultMatrixf(rotationMatrix.columnMajor().coordinates);
} // Rotatef()

// scale the matrix
void FakeGL::Scalef(float xScale, float yScale, float zScale)
{ // Scalef()
    // Scales along the x, y and z axes
    // Applied to CURRENT matrix, though generally used for modelview
    Matrix4 scalingMatrix;
    scalingMatrix.SetScale(xScale, yScale, zScale);
    // Multiply by the current matrix
    MultMatrixf(scalingMatrix.columnMajor().coordinates);
} // Scalef()

// translate the matrix
void FakeGL::Translatef(float xTranslate, float yTranslate, float zTranslate)
{ // Translatef()
    // Declare a new matrix and set it to the identity
    Matrix4 translation;
    translation.SetIdentity();
    // Change the translation part (col 4, rows 1-3)
    translation[0][3] = xTranslate;
    translation[1][3] = yTranslate;
    translation[2][3] = zTranslate;
    // Multiply with the active matrix
    MultMatrixf(translation.columnMajor().coordinates);
} // Translatef()

// sets the viewport
void FakeGL::Viewport(int x, int y, int width, int height)
{ // Viewport()
    // OpenGL allows negative values for x and y, but not width and height
    // I will not change anything if any value is incorrect
    if (width > 0 && height > 0)
    {
        // Set viewport values
        viewportX = x;
        viewportY = y;
        viewportWidth = width;
        viewportHeight = height;
    }
} // Viewport()

//-------------------------------------------------//
//                                                 //
// VERTEX ATTRIBUTE ROUTINES                       //
//                                                 //
//-------------------------------------------------//

// sets colour with floating point
void FakeGL::Color3f(float red, float green, float blue)
{ // Color3f()
    // Our colours are stored in the range 0-255, but specified as 0-1, so multiply
    float newRed = red * 255.0, newGreen = green * 255.0, newBlue = blue * 255.0;
    // Clip each colour
    for (auto& colour : {&newRed, &newGreen, &newBlue})
    {
        if (*colour < 0.0f)
        {
            *colour = 0.0f;
        }
        else if (*colour > 255.0f)
        {
            *colour = 255.0f;
        }
    }
    // Set the colour
    currentColor = RGBAValue(newRed, newGreen, newBlue, 255.0f);
} // Color3f()

// sets material properties
void FakeGL::Materialf(unsigned int parameterName, const float parameterValue)
{ // Materialf()
    // As OpenGL spec, this may only update the shininess component
    if (parameterName == FAKEGL_SHININESS)
    {
        // Also clamp as per spec
        if (parameterValue >= 0.0f && parameterValue <= 128.0f)
        {
            currentShininess = parameterValue;
        }
        // OpenGL spec implies putting a wrong value spits out an error and nothing else, so I will do nothing else
    }
} // Materialf()

void FakeGL::Materialfv(unsigned int parameterName, const float *parameterValues)
{ // Materialfv()
    // Check the requested param
    switch (parameterName)
    {
        // Set the appropriate property
    case FAKEGL_AMBIENT:
        currentAmbientReflectance = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
        break;
    case FAKEGL_DIFFUSE:
        currentDiffuseReflectance = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
        break;
    case FAKEGL_AMBIENT_AND_DIFFUSE:
        // Equivalent to calling this function twice
        Materialfv(FAKEGL_AMBIENT, parameterValues);
        Materialfv(FAKEGL_DIFFUSE, parameterValues);
        break;
    case FAKEGL_SPECULAR:
        currentSpecularReflectance = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
        break;
    case FAKEGL_EMISSION:
        currentEmission = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
        break;
    case FAKEGL_SHININESS:
        // Call the other material function, it does the same
        Materialf(FAKEGL_SHININESS, *parameterValues);
        break;
    }
} // Materialfv()

// sets the normal vector
void FakeGL::Normal3f(float x, float y, float z)
{ // Normal3f()
    currentNormal.x = x;
    currentNormal.y = y;
    currentNormal.z = z;
    currentNormal.w = 1.0f;
} // Normal3f()

// sets the texture coordinates
void FakeGL::TexCoord2f(float u, float v)
{ // TexCoord2f()
    // OpenGL texture behaviour can be changed, default behaviour however is to repeat
    // As there is no way to change this, I will hardcode the default here
    if (u > 1.0f)
    {
        u = u - (int)u;
    }
    else if (u < 0.0f)
    {
        u = 1.0f - (u - (int)u);
    }
    if (v > 1.0f)
    {
        v = v - (int)v;
    }
    else if (v < 0.0f)
    {
        v = 1.0f - (v - (int)v);
    }
    // Now set current uvs
    currentU = u;
    currentV = v;
} // TexCoord2f()

// sets the inVertex & launches it down the pipeline
void FakeGL::Vertex3f(float x, float y, float z)
    { // Vertex3f()
    // Declare a temporary inVertex, use for vertexWithAttributes and add to the back of the inVertex queue
    Homogeneous4 myVertex;
    myVertex.x = x;
    myVertex.y = y;
    myVertex.z = z;
    myVertex.w = 1.0f;

    // Use position of temporary inVertex for attributed
    vertexWithAttributes inVertex;
    inVertex.position = myVertex;

    // Set to current colour
    inVertex.colour = currentColor;

    // Set to current normals
    inVertex.normal = currentNormal;

    // Set the material properties
    inVertex.ambientReflectance = currentAmbientReflectance;
    inVertex.diffuseReflectance = currentDiffuseReflectance;
    inVertex.specularReflectance = currentSpecularReflectance;
    inVertex.emission = currentEmission;
    inVertex.shininess = currentShininess;

    // Set the texture coordinates
    inVertex.u = currentU;
    inVertex.v = currentV;

    // Push to the queue, and for the moment, call the next pipeline stage
    vertexQueue.push_back(inVertex);

    FakeGL::TransformVertex();
    } // Vertex3f()

//-------------------------------------------------//
//                                                 //
// STATE VARIABLE ROUTINES                         //
//                                                 //
//-------------------------------------------------//

// disables a specific flag in the library
void FakeGL::Disable(unsigned int property)
    { // Disable()
    // Check the property against the flags, and change to false if they match
    switch (property)
    {
    case FAKEGL_LIGHTING:
        lighting = false;
        break;
    case FAKEGL_TEXTURE_2D:
        textureEnabled = false;
        break;
    case FAKEGL_DEPTH_TEST:
        depthTest = false;
        break;
    case FAKEGL_PHONG_SHADING:
        phongShading = false;
        break;
    }
    } // Disable()

// enables a specific flag in the library
void FakeGL::Enable(unsigned int property)
{ // Enable()
    // Check the property against the flags, and change to true if they match
    switch (property)
    {
    case FAKEGL_LIGHTING:
        lighting = true;
        break;
    case FAKEGL_TEXTURE_2D:
        textureEnabled = true;
        break;
    case FAKEGL_DEPTH_TEST:
        depthTest = true;
        break;
    case FAKEGL_PHONG_SHADING:
        phongShading = true;
        break;
    }
} // Enable()

//-------------------------------------------------//
//                                                 //
// LIGHTING STATE ROUTINES                         //
//                                                 //
//-------------------------------------------------//

// sets properties for the one and only light
void FakeGL::Light(int parameterName, const float *parameterValues)
{ // Light()
    // Check which property to set
    if ((parameterName & FAKEGL_AMBIENT_AND_DIFFUSE) == FAKEGL_AMBIENT_AND_DIFFUSE)
    {
        // Set ambient and diffuse to the same
        ambientColor = diffuseColor = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
    }
    else if ((parameterName & FAKEGL_AMBIENT) == FAKEGL_AMBIENT)
    {
        // Set ambient colour
        ambientColor = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
    }
    else if ((parameterName & FAKEGL_DIFFUSE) == FAKEGL_DIFFUSE)
    {
        // Set diffuse colour
        diffuseColor = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
    }
    else if ((parameterName & FAKEGL_SPECULAR) == FAKEGL_SPECULAR)
    {
        // Set specular light
        specularColor = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
    }
    else if ((parameterName & FAKEGL_POSITION) == FAKEGL_POSITION)
    {
        // Set the light position
        Homogeneous4 lightPosTemp(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
        // Multiply by the current marix and set pos
        lightPos = modelViewMat * lightPosTemp;
    }
} // Light()

//-------------------------------------------------//
//                                                 //
// TEXTURE PROCESSING ROUTINES                     //
//                                                 //
// Note that we only allow one texture             //
// so glGenTexture & glBindTexture aren't needed   //
//                                                 //
//-------------------------------------------------//

// sets whether textures replace or modulate
void FakeGL::TexEnvMode(unsigned int textureMode)
{ // TexEnvMode()
    // Check the texture mode, and change the flag appropriately
    if (textureMode == FAKEGL_REPLACE)
    {
        FakeGL::textureMode = FAKEGL_REPLACE;
    }
    else if (textureMode == FAKEGL_MODULATE)
    {
        FakeGL::textureMode = FAKEGL_MODULATE;
    }
} // TexEnvMode()

// sets the texture image that corresponds to a given ID
void FakeGL::TexImage2D(const RGBAImage &textureImage)
{ // TexImage2D()
    // Check the size is greater than 0 in both dimensions (avoiding divisions by 0 later)
    if (textureImage.height > 0 && textureImage.width > 0)
    {
        currentTexture = &textureImage;
    }
} // TexImage2D()

//-------------------------------------------------//
//                                                 //
// FRAME BUFFER ROUTINES                           //
//                                                 //
//-------------------------------------------------//

// clears the frame buffer
void FakeGL::Clear(unsigned int mask)
    { // Clear()
      // Check mask, indicating if colour and/or depth buffer should be cleared
    if (mask & FAKEGL_COLOR_BUFFER_BIT)
    {
        // Set all pixels to previously set color value
        for (size_t row = 0; row < frameBuffer.height; row++)
        {
            for (size_t col = 0; col < frameBuffer.width; col++)
            {
                frameBuffer[row][col] = clearColorVal;
            }
        }
    }
    if (mask & FAKEGL_DEPTH_BUFFER_BIT)
    {
        // Set all per pixel depth values to previously set depth value
        for (size_t row = 0; row < frameBuffer.height; row++)
        {
            for (size_t col = 0; col < frameBuffer.width; col++)
            {
                depthBuffer[row][col] = depthVal;
            }
        }
    }
    } // Clear()

// sets the clear colour for the frame buffer
void FakeGL::ClearColor(float red, float green, float blue, float alpha)
    { // ClearColor()
    // Multiply by 255 to convert from OpenGL format (float, [0-1]) to FakeGL (float, 0-255)
    clearColorVal = RGBAValue(red * 255, green * 255, blue * 255, alpha * 255);
    } // ClearColor()

//-------------------------------------------------//
//                                                 //
// MAJOR PROCESSING ROUTINES                       //
//                                                 //
//-------------------------------------------------//

// transform one inVertex & shift to the raster queue
void FakeGL::TransformVertex()
{ // TransformVertex()
    // Loop until the queue is empty
    while (!vertexQueue.empty())
    {
        // DCS vertex to use
        screenVertexWithAttributes myScreenVertex;

        // Get a inVertex from the queue
        vertexWithAttributes inVertex = vertexQueue.back();
        vertexQueue.pop_back();

        // Apply transformations: ModelView first, then projection
        inVertex.position = modelViewMat * inVertex.position;
        inVertex.position = projectionMat * inVertex.position;
        // Transform normals
        inVertex.normal = modelViewMat * inVertex.normal;
        inVertex.normal = projectionMat * inVertex.normal;
        
        // Perspective divison (by w): Convert to cartesian then divide as scalar
        Cartesian3 cartesian(inVertex.position.x, inVertex.position.y, inVertex.position.z);
        cartesian = cartesian / inVertex.position.w;

        // Lighting calculations (if requested). Don't waste time calculating them if phong shading.
        if (lighting && !phongShading)
        {
            // Tempoary Homogenous4 for light floats
            Homogeneous4 light(0.0f, 0.0f, 0.0f, 0.0f);
            // Ambient can be applied straight away
            light = Homogeneous4(ambientColor.x * inVertex.ambientReflectance.x, ambientColor.y * inVertex.ambientReflectance.y, ambientColor.z * inVertex.ambientReflectance.z, 0.0f);
            
            // Diffuse based on light pos (incidence angle)
            Cartesian3 lightPosTmp = lightPos.Vector().unit();
            Cartesian3 lightDirection(lightPosTmp.x, lightPosTmp.y, -lightPosTmp.z);
            Cartesian3 surfaceNormal = inVertex.normal.Vector().unit();
            // Set diffuse amount based on incidence angle
            float diffuseAmount = surfaceNormal.dot(lightDirection);
            if (diffuseAmount > 0.0f)
            {
                light = light + (Homogeneous4(diffuseColor.x * inVertex.diffuseReflectance.x, diffuseColor.y * inVertex.diffuseReflectance.y, diffuseColor.z * inVertex.diffuseReflectance.z, diffuseColor.w * inVertex.diffuseReflectance.w) * diffuseAmount);
            }
            
            // Set specular based on incidence angle, viewing angle and shininess exponent
            // 'Camera' is at 0,0,0
            Cartesian3 eyeVec = Cartesian3(0.0f, 0.0f, 0.0f) - Cartesian3(0.0f, 0.0f, 0.0f);
            //Cartesian3 specularLightDirection = (lightPos.Vector() - inVertex.position.Vector()).unit();
            Cartesian3 bisector = ((eyeVec + lightDirection) / 2.0f ).unit();
            // Calculate specular
            // Check if the dot product is negative before raising to exponent, to avoid negatives becoming positives
            float dotProduct = surfaceNormal.dot(bisector);
            if (dotProduct < 0.0f)
            {
                dotProduct = 0.0f;
            }
            float specularAmount = pow(dotProduct, inVertex.shininess);
            // If there is any specular, add it to the light
            if (specularAmount > 0.0f)
            {
                light = light + (Homogeneous4(specularColor.x * inVertex.specularReflectance.x, specularColor.y * inVertex.specularReflectance.y, specularColor.z * inVertex.specularReflectance.z, 0.0f) * specularAmount);
            }
            
            // Add emission
            light = light + inVertex.emission;

            // Clip values
            for (auto& component : { &light.x, &light.y, &light.z, &light.w })
            {
                if (*component > 1.0f)
                {
                    *component = 1.0f;
                }
                else if (*component < 0.0f)
                {
                    *component = 0.0f;
                }
            }

            // Convert to colour
            myScreenVertex.colour = RGBAValue(255.0f * light.x, 255.0f * light.y, 255.0f * light.z, 255.0f * light.w);
        }
        else
        {
            // Pass the colour through
            myScreenVertex.colour = inVertex.colour;
        }
        // Convert NDC to DC
        myScreenVertex.position.x = viewportX + ((cartesian.x + 1.0) * 0.5 * viewportWidth);
        myScreenVertex.position.y = viewportY + ((cartesian.y + 1.0) * 0.5 * viewportHeight);
        // Change depth to range 0-255, as our depth buffer is an 8 bit int
        // This does not prevent against values outside this range, as clipping is performed in raster
        myScreenVertex.position.z = std::round(((cartesian.z - zFar) / (zNear - zFar)) * 255.0);

        // Pass through remaining attributes
        myScreenVertex.normal = inVertex.normal;
        myScreenVertex.ambientReflectance = inVertex.ambientReflectance;
        myScreenVertex.diffuseReflectance = inVertex.diffuseReflectance;
        myScreenVertex.specularReflectance = inVertex.specularReflectance;
        myScreenVertex.emission = inVertex.emission;
        myScreenVertex.shininess = inVertex.shininess;
        myScreenVertex.u = inVertex.u;
        myScreenVertex.v = inVertex.v;

        rasterQueue.push_back(myScreenVertex);
    }
} // TransformVertex()

// rasterise a single primitive if there are enough vertices on the queue
bool FakeGL::RasterisePrimitive()
{ // RasterisePrimitive(
    // Check the primitive
    // Results where a primitive hasn't been set are undefined in OpenGL- I will refuse 
    //      to rasterise anything without a primitive being set.
    if (primitiveMode == FAKEGL_POINTS)
    {
        // Check the queue size
        if (rasterQueue.size() >= 1)
        {
            // Primitive can be drawn, pop a inVertex
            screenVertexWithAttributes inVertex = rasterQueue.front();
            rasterQueue.pop_front();
            // Now call the appropriate drawing function and return true
            RasterisePoint(inVertex);
            return true;
        }
        else
        {
            // There are not enough vertices on the queue, return false
            return false;
        }
    }
    else if (primitiveMode == FAKEGL_LINES)
    {
        // Check the queue size
        if (rasterQueue.size() >= 2)
        {
            // Primitive can be drawn, pop two vertices
            screenVertexWithAttributes vertex0 = rasterQueue.front();
            rasterQueue.pop_front();
            screenVertexWithAttributes vertex1 = rasterQueue.front();
            rasterQueue.pop_front();
            // Now call the appropriate drawing function and return true
            RasteriseLineSegment(vertex0, vertex1);
            return true;
        }
        else
        {
            // There are not enough vertices on the queue, return false
            return false;
        }
    }
    else if (primitiveMode == FAKEGL_TRIANGLES)
    {
        // Check the queue size
        if (rasterQueue.size() >= 3)
        {
            // Primitive can be drawn, pop three vertices
            screenVertexWithAttributes vertex0 = rasterQueue.front();
            rasterQueue.pop_front();
            screenVertexWithAttributes vertex1 = rasterQueue.front();
            rasterQueue.pop_front();
            screenVertexWithAttributes vertex2 = rasterQueue.front();
            rasterQueue.pop_front();
            // Now call the appropriate drawing function and return true
            RasteriseTriangle(vertex0, vertex1, vertex2);
            return true;
        }
        else
        {
            // There are not enough vertices on the queue, return false
            return false;
        }
    }
    // Return false as a primitive was not specified
    return false;
} // RasterisePrimitive()

// rasterises a single point
void FakeGL::RasterisePoint(screenVertexWithAttributes &vertex0)
{ // RasterisePoint()
    // Find which pixels intersecting the point- depends on point size, default 1
    // Calculate the min and max columns from the point size. Round to ints.
    int rowMin = std::round(vertex0.position.y - pointSize / 2.0);
    int rowMax = std::round(vertex0.position.y + pointSize / 2.0);
    int colMin = std::round(vertex0.position.x - pointSize / 2.0);
    int colMax = std::round(vertex0.position.x + pointSize / 2.0);

    // Create a fragment for reuse
    fragmentWithAttributes rasterFragment;
    // Set the colour now, as it'll be the same for every pixel on a point
    rasterFragment.colour = vertex0.colour;
    // Depth
    RGBAValue depth(0.0f, 0.0f, 0.0f, 0.0f);

    // Colour fragments and send to fragment queue
    for (rasterFragment.row = rowMin; rasterFragment.row < rowMax; rasterFragment.row++)
    {
        // Clipping in rows
        if (rasterFragment.row < 0 || rasterFragment.row > frameBuffer.height)
        {
            continue;
        }
        // Inner loop
        for (rasterFragment.col = colMin; rasterFragment.col < colMax; rasterFragment.col++)
        {
            // Clipping in columns
            if (rasterFragment.col < 0 || rasterFragment.col > frameBuffer.width)
            {
                continue;
            }
            // No interpolation is neccasary for points, so everything is done, push the fragment to the fragment queue
            fragmentQueue.push_back(rasterFragment);
            // Store depth in alpha, clipping first
            if (vertex0.position.z < 0.0 || vertex0.position.z > 255.0)
            {
                depth.alpha = 255.0;
            }
            else
            {
                depth.alpha = vertex0.position.z;
            }
            depthBuffer[rasterFragment.row][rasterFragment.col] = depth;

            // Might want to move this :O
            //ProcessFragment();
        }
    }
} // RasterisePoint()

// rasterises a single line segment
void FakeGL::RasteriseLineSegment(screenVertexWithAttributes &vertex0, screenVertexWithAttributes &vertex1)
{ // RasteriseLineSegment()
    // Draw as quads
    // Get vector and unit normal
    Cartesian3 vec = vertex1.position - vertex0.position;
    Cartesian3 normal = vec;
    normal.x = -vec.y;
    normal.y = vec.x;
    normal = normal.unit();
    // Rasterise quad as two triangles
    // Generate new vertices
    screenVertexWithAttributes q0 = vertex0, q1 = vertex0, p0 = vertex1, p1 = vertex1;
    q0.position = vertex0.position - ((lineWidth / 2) * normal);
    q1.position = vertex0.position + ((lineWidth / 2) * normal);
    p0.position = vertex1.position - ((lineWidth / 2) * normal);
    p1.position = vertex1.position + ((lineWidth / 2) * normal);
    // Rasterise as two triangles
    RasteriseTriangle(q0, p0, p1);
    RasteriseTriangle(p1, q1, q0);
} // RasteriseLineSegment()

// rasterises a single triangle
void FakeGL::RasteriseTriangle(screenVertexWithAttributes &vertex0, screenVertexWithAttributes &vertex1, screenVertexWithAttributes &vertex2)
    { // RasteriseTriangle()
    // compute a bounding box that starts inverted to frame size
    // clipping will happen in the raster loop proper
    float minX = frameBuffer.width, maxX = 0.0;
    float minY = frameBuffer.height, maxY = 0.0;
    
    // test against all vertices
    if (vertex0.position.x < minX) minX = vertex0.position.x;
    if (vertex0.position.x > maxX) maxX = vertex0.position.x;
    if (vertex0.position.y < minY) minY = vertex0.position.y;
    if (vertex0.position.y > maxY) maxY = vertex0.position.y;
    
    if (vertex1.position.x < minX) minX = vertex1.position.x;
    if (vertex1.position.x > maxX) maxX = vertex1.position.x;
    if (vertex1.position.y < minY) minY = vertex1.position.y;
    if (vertex1.position.y > maxY) maxY = vertex1.position.y;
    
    if (vertex2.position.x < minX) minX = vertex2.position.x;
    if (vertex2.position.x > maxX) maxX = vertex2.position.x;
    if (vertex2.position.y < minY) minY = vertex2.position.y;
    if (vertex2.position.y > maxY) maxY = vertex2.position.y;

    // now for each side of the triangle, compute the line vectors
    Cartesian3 vector01 = vertex1.position - vertex0.position;
    Cartesian3 vector12 = vertex2.position - vertex1.position;
    Cartesian3 vector20 = vertex0.position - vertex2.position;

    // now compute the line normal vectors
    Cartesian3 normal01(-vector01.y, vector01.x, 0.0);  
    Cartesian3 normal12(-vector12.y, vector12.x, 0.0);  
    Cartesian3 normal20(-vector20.y, vector20.x, 0.0);  

    // we don't need to normalise them, because the square roots will cancel out in the barycentric coordinates
    float lineConstant01 = normal01.dot(vertex0.position);
    float lineConstant12 = normal12.dot(vertex1.position);
    float lineConstant20 = normal20.dot(vertex2.position);

    // and compute the distance of each inVertex from the opposing side
    float distance0 = normal12.dot(vertex0.position) - lineConstant12;
    float distance1 = normal20.dot(vertex1.position) - lineConstant20;
    float distance2 = normal01.dot(vertex2.position) - lineConstant01;

    // if any of these are zero, we will have a divide by zero error
    // but notice that if they are zero, the vertices are collinear in projection and the triangle is edge on
    // we can render that as a line, but the better solution is to render nothing.  In a surface, the adjacent
    // triangles will eventually take care of it
    if ((distance0 == 0) || (distance1 == 0) || (distance2 == 0))
        return; 

    // create a fragment for reuse
    fragmentWithAttributes rasterFragment;

    // Copy colours
    RGBAValue v0Color = vertex0.colour, v1Color = vertex1.colour, v2Color = vertex2.colour;

    // loop through the pixels in the bounding box
    for (rasterFragment.row = minY; rasterFragment.row <= maxY; rasterFragment.row++)
        { // per row
        // this is here so that clipping works correctly
        if (rasterFragment.row < 0) continue;
        if (rasterFragment.row >= frameBuffer.height) continue;
        for (rasterFragment.col = minX; rasterFragment.col <= maxX; rasterFragment.col++)
            { // per pixel
            // this is also for correct clipping
            if (rasterFragment.col < 0) continue;
            if (rasterFragment.col >= frameBuffer.width) continue;
            
            // the pixel in cartesian format
            Cartesian3 pixel(rasterFragment.col, rasterFragment.row, 0.0);
            
            // right - we have a pixel inside the frame buffer AND the bounding box
            // note we *COULD* compute gamma = 1.0 - alpha - beta instead
            float alpha = (normal12.dot(pixel) - lineConstant12) / distance0;           
            float beta = (normal20.dot(pixel) - lineConstant20) / distance1;            
            float gamma = (normal01.dot(pixel) - lineConstant01) / distance2;           

            // now perform the half-plane test
            if ((alpha < 0.0) || (beta < 0.0) || (gamma < 0.0))
                continue;

            
            // Compute color by specified method
            if (!phongShading)
            {
                // No phong shading, lighting (if applicable) was already computed and stored in colour, so pass it through
                rasterFragment.colour = alpha * v0Color + beta * v1Color + gamma * v2Color;
            }
            else if (phongShading && lighting)
            {
                // Phong shading- must interpolate all values
                // Ambient Reflectance
                float ambientReflectanceRed = alpha * vertex0.ambientReflectance.x + beta * vertex1.ambientReflectance.x + gamma * vertex2.ambientReflectance.x;
                float ambientReflectanceGreen = alpha * vertex0.ambientReflectance.y + beta * vertex1.ambientReflectance.y + gamma * vertex2.ambientReflectance.y;
                float ambientReflectanceBlue = alpha * vertex0.ambientReflectance.z + beta * vertex1.ambientReflectance.z + gamma * vertex2.ambientReflectance.z;
                // Diffuse Reflectance
                float diffuseReflectanceRed = alpha * vertex0.diffuseReflectance.x + beta * vertex1.diffuseReflectance.x + gamma * vertex2.diffuseReflectance.x;
                float diffuseReflectanceGreen = alpha * vertex0.diffuseReflectance.y + beta * vertex1.diffuseReflectance.y + gamma * vertex2.diffuseReflectance.y;
                float diffuseReflectanceBlue = alpha * vertex0.diffuseReflectance.z + beta * vertex1.diffuseReflectance.z + gamma * vertex2.diffuseReflectance.z;
                float diffuseReflectanceAlpha = alpha * vertex0.diffuseReflectance.w + beta * vertex1.diffuseReflectance.w + gamma * vertex2.diffuseReflectance.w;
                // Specular Reflectance
                float specularReflectanceRed = alpha * vertex0.specularReflectance.x + beta * vertex1.specularReflectance.x + gamma * vertex2.specularReflectance.x;
                float specularReflectanceGreen = alpha * vertex0.specularReflectance.y + beta * vertex1.specularReflectance.y + gamma * vertex2.specularReflectance.y;
                float specularReflectanceBlue = alpha * vertex0.specularReflectance.z + beta * vertex1.specularReflectance.z + gamma * vertex2.specularReflectance.z;
                // Emission
                float emissionRed = alpha * vertex0.emission.x + beta * vertex1.emission.x + gamma * vertex2.emission.x;
                float emissionGreen = alpha * vertex0.emission.y + beta * vertex1.emission.y + gamma * vertex2.emission.y;
                float emissionBlue = alpha * vertex0.emission.z + beta * vertex1.emission.z + gamma * vertex2.emission.z;
                // Shininess
                float shininess = alpha * vertex0.shininess + beta * vertex1.shininess + gamma * vertex2.shininess;
                // Normal
                Cartesian3 fragmentNormal = alpha * vertex0.normal.Vector().unit() + beta * vertex1.normal.Vector().unit() + gamma * vertex2.normal.Vector().unit();
                
                // Tempoary Homogenous4 for light floats
                Homogeneous4 light(0.0f, 0.0f, 0.0f, 0.0f);
                // Ambient can be applied straight away
                light = Homogeneous4(ambientColor.x * ambientReflectanceRed, ambientColor.y * ambientReflectanceGreen, ambientColor.z * ambientReflectanceBlue, 0.0f);

                // Diffuse based on light pos (incidence angle)
                Cartesian3 lightPosTmp = lightPos.Vector().unit();
                Cartesian3 lightDirection(lightPosTmp.x, lightPosTmp.y, -lightPosTmp.z);
                // Set diffuse amount based on incidence angle
                float diffuseAmount = fragmentNormal.dot(lightDirection);
                if (diffuseAmount > 0.0f)
                {
                    light = light + (Homogeneous4(diffuseColor.x * diffuseReflectanceRed, diffuseColor.y * diffuseReflectanceGreen, diffuseColor.z * diffuseReflectanceBlue, diffuseColor.w * diffuseReflectanceAlpha) * diffuseAmount);
                }

                // Set specular based on incidence angle, viewing angle and shininess exponent
                // 'Camera' is at 0,0,0
                Cartesian3 eyeVec = Cartesian3(0.0f, 0.0f, 0.0f);
                //Cartesian3 specularLightDirection = (lightPos.Vector() - inVertex.position.Vector()).unit();
                Cartesian3 bisector = ((eyeVec + lightDirection) / 2.0f).unit();
                // Calculate specular
                // Check if the dot product is negative before raising to exponent, to avoid negatives becoming positives
                float dotProduct = fragmentNormal.dot(bisector);
                if (dotProduct < 0.0f)
                {
                    dotProduct = 0.0f;
                }
                float specularAmount = pow(dotProduct, shininess);
                // If there is any specular, add it to the light
                if (specularAmount > 0.0f)
                {
                    light = light + (Homogeneous4(specularColor.x * specularReflectanceRed, specularColor.y * specularReflectanceGreen, specularColor.z * specularReflectanceBlue, 0.0f) * specularAmount);
                }

                // Add emission
                light = light + Homogeneous4(emissionRed, emissionGreen, emissionBlue, 0.0f);

                // Clip values
                for (auto& component : { &light.x, &light.y, &light.z, &light.w })
                {
                    if (*component > 1.0f)
                    {
                        *component = 1.0f;
                    }
                    else if (*component < 0.0f)
                    {
                        *component = 0.0f;
                    }
                }

                // Convert to colour
                rasterFragment.colour = RGBAValue(255.0f * light.x, 255.0f * light.y, 255.0f * light.z, 255.0f * light.w);
            }
            // compute textures
            if (textureEnabled)
            {
                // Interpolate UVs
                float fragU = alpha * vertex0.u + beta * vertex1.u + gamma * vertex2.u;
                float fragV = alpha * vertex0.v + beta * vertex1.v + gamma * vertex2.v;
                // Convert to discrete texture coords. Flip u/col as OpenGL interpolates from bottom left, while texture is stored from top left.
                int texCol = std::round(fragU * currentTexture->width);
                int texRow = std::round(fragV * currentTexture->height);
                // Check interpolation mode
                if (textureMode == FAKEGL_MODULATE)
                {
                    // Modulate texture: multiply by current colours to apply lighting
                    rasterFragment.colour.red = std::round(((float)rasterFragment.colour.red / 255.0f) * (float)((*currentTexture)[texRow][texCol]).red);
                    rasterFragment.colour.green = std::round(((float)rasterFragment.colour.green / 255.0f) * (float)((*currentTexture)[texRow][texCol]).green);
                    rasterFragment.colour.blue = std::round(((float)rasterFragment.colour.blue / 255.0f) * (float)((*currentTexture)[texRow][texCol]).blue);
                    rasterFragment.colour.alpha = std::round(((float)rasterFragment.colour.alpha / 255.0f) * (float)((*currentTexture)[texRow][texCol]).alpha);
                }
                else
                {
                    // Replace texture (ignore lighting)
                    rasterFragment.colour = (*currentTexture)[texRow][texCol];
                }
            }
            // compute depth
            float depth = alpha * vertex0.position.z + beta * vertex1.position.z + gamma * vertex2.position.z;
            if (depth > 255.0)
            {
                depth = 255.0;
            }
            else if (depth < 0.0)
            {
                depth = 0.0;
            }
            rasterFragment.depth.alpha = (char)std::round(depth);

            // now we add it to the queue for fragment processing
            fragmentQueue.push_back(rasterFragment);
            } // per pixel
        } // per row
    } // RasteriseTriangle()

// process a single fragment
void FakeGL::ProcessFragment()
{ // ProcessFragment()
    // Check a fragment exists, return if not
    if (!(fragmentQueue.size() >= 1))
    {
        return;
    }
    // Get the fragment from the queue and pop it
    fragmentWithAttributes fragment = fragmentQueue.front();
    fragmentQueue.pop_front();

    // Depth testing ( if enabled )
    if (depthTest)
    {
        // Check if this fragment is in front
        if (depthBuffer[fragment.row][fragment.col].alpha < fragment.depth.alpha)
        {
            // This fragment is in front, update depth buffer and set in framebuffer
            depthBuffer[fragment.row][fragment.col].alpha = fragment.depth.alpha;
            // Set the value in the framebuffer- relying on clipping from earlier
            frameBuffer[fragment.row][fragment.col] = fragment.colour;
        }
    }
    else
    {
        // Set the value in the framebuffer- relying on clipping from earlier
        frameBuffer[fragment.row][fragment.col] = fragment.colour;
    }
} // ProcessFragment()

// standard routine for dumping the entire FakeGL context (except for texture / image)
std::ostream &operator << (std::ostream &outStream, FakeGL &fakeGL)
    { // operator <<
    outStream << "=========================" << std::endl;
    outStream << "Dumping FakeGL Context   " << std::endl;
    outStream << "=========================" << std::endl;


    outStream << "-------------------------" << std::endl;
    outStream << "Vertex Queue:            " << std::endl;
    outStream << "-------------------------" << std::endl;
    for (auto inVertex = fakeGL.vertexQueue.begin(); inVertex < fakeGL.vertexQueue.end(); inVertex++)
        { // per matrix
        outStream << "Vertex " << inVertex - fakeGL.vertexQueue.begin() << std::endl;
        outStream << *inVertex;
        } // per matrix


    outStream << "-------------------------" << std::endl;
    outStream << "Raster Queue:            " << std::endl;
    outStream << "-------------------------" << std::endl;
    for (auto inVertex = fakeGL.rasterQueue.begin(); inVertex < fakeGL.rasterQueue.end(); inVertex++)
        { // per matrix
        outStream << "Vertex " << inVertex - fakeGL.rasterQueue.begin() << std::endl;
        outStream << *inVertex;
        } // per matrix


    outStream << "-------------------------" << std::endl;
    outStream << "Fragment Queue:          " << std::endl;
    outStream << "-------------------------" << std::endl;
    for (auto fragment = fakeGL.fragmentQueue.begin(); fragment < fakeGL.fragmentQueue.end(); fragment++)
        { // per matrix
        outStream << "Fragment " << fragment - fakeGL.fragmentQueue.begin() << std::endl;
        outStream << *fragment;
        } // per matrix


    return outStream;
    } // operator <<

// subroutines for other classes
std::ostream &operator << (std::ostream &outStream, vertexWithAttributes &inVertex)
    { // operator <<
    std::cout << "Vertex With Attributes" << std::endl;
    std::cout << "Position:   " << inVertex.position << std::endl;
    std::cout << "Colour:     " << inVertex.colour << std::endl;

    // you

    return outStream;
    } // operator <<

std::ostream &operator << (std::ostream &outStream, screenVertexWithAttributes &inVertex) 
    { // operator <<
    std::cout << "Screen Vertex With Attributes" << std::endl;
    std::cout << "Position:   " << inVertex.position << std::endl;
    std::cout << "Colour:     " << inVertex.colour << std::endl;

    return outStream;
    } // operator <<

std::ostream &operator << (std::ostream &outStream, fragmentWithAttributes &fragment)
    { // operator <<
    std::cout << "Fragment With Attributes" << std::endl;
    std::cout << "Row:        " << fragment.row << std::endl;
    std::cout << "Col:        " << fragment.col << std::endl;
    std::cout << "Colour:     " << fragment.colour << std::endl;

    return outStream;
    } // operator <<


    
