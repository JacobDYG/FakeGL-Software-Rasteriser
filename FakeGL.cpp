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
    // Lighting values
    lightPos(0.0f, 0.0f, 0.0f, 1.0f),
    ambientColor(1.0f, 1.0f, 1.0f, 1.0f),
    diffuseColor(1.0f, 1.0f, 1.0f, 1.0f),
    specularColor(1.0f, 1.0f, 1.0f, 1.0f)
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
    // Set the primitive value to -1
    primitiveMode = -1;
    // Indicate that primitive specification is done
    primitiveAssembly = false;
    // This could potentially be used to call shader stages
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
    } // Frustum()

// sets an orthographic projection matrix
void FakeGL::Ortho(float left, float right, float bottom, float top, float zNear, float zFar)
    { // Ortho()
    } // Ortho()

// rotate the matrix
void FakeGL::Rotatef(float angle, float axisX, float axisY, float axisZ)
    { // Rotatef()
    } // Rotatef()

// scale the matrix
void FakeGL::Scalef(float xScale, float yScale, float zScale)
    { // Scalef()
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
    } // Materialf()

void FakeGL::Materialfv(unsigned int parameterName, const float *parameterValues)
    { // Materialfv()
    } // Materialfv()

// sets the normal vector
void FakeGL::Normal3f(float x, float y, float z)
    { // Normal3f()
    } // Normal3f()

// sets the texture coordinates
void FakeGL::TexCoord2f(float u, float v)
    { // TexCoord2f()
    } // TexCoord2f()

// sets the vertex & launches it down the pipeline
void FakeGL::Vertex3f(float x, float y, float z)
    { // Vertex3f()
    // Declare a tempoary vertex, use for vertexWithAttributes and add to the back of the vertex queue
    Homogeneous4 myVertex;
    myVertex.x = x;
    myVertex.y = y;
    myVertex.z = z;
    myVertex.w = 1.0f;

    vertexWithAttributes myVertexWithAttributes;
    myVertexWithAttributes.position = myVertex;

    // Set to current colour
    myVertexWithAttributes.colour = currentColor;

    vertexQueue.push_back(myVertexWithAttributes);

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
        textures = false;
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
        textures = true;
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
        lightPos = Homogeneous4(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3]);
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
    } // TexEnvMode()

// sets the texture image that corresponds to a given ID
void FakeGL::TexImage2D(const RGBAImage &textureImage)
    { // TexImage2D()
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

// transform one vertex & shift to the raster queue
void FakeGL::TransformVertex()
    { // TransformVertex()
    // Loop until the queue is empty
    while (!vertexQueue.empty())
    {
        // Get a vertex from the queue
        vertexWithAttributes myVertexWithAttributes = vertexQueue.back();
        vertexQueue.pop_back();

        // Apply transformations: ModelView first, then projection
        myVertexWithAttributes.position = modelViewMat * myVertexWithAttributes.position;
        myVertexWithAttributes.position = projectionMat * myVertexWithAttributes.position;

        // Perspective divison (by w)
        screenVertexWithAttributes myScreenVertex;
        // Convert to cartesian then divide as scalar
        Cartesian3 cartesian(myVertexWithAttributes.position.x, myVertexWithAttributes.position.y, myVertexWithAttributes.position.z);
        cartesian = cartesian / myVertexWithAttributes.position.w;
        
        // Convert NDC to DC
        myScreenVertex.position.x = viewportX + ((cartesian.x + 1.0) * 0.5 * viewportHeight);
        myScreenVertex.position.y = viewportY + ((cartesian.y + 1.0) * 0.5 * viewportWidth);
        // Change depth to range 0-255, as our depth buffer is an 8 bit int
        // This does not prevent against values outside this range, as clipping is performed in raster
        myScreenVertex.position.z = std::round(((cartesian.z - zFar) / (zNear - zFar)) * 255.0);

        // Lighting calculations (if requested)
        if (lighting)
        {
            // Ambient can be applied straight away
            myScreenVertex.colour.red = (char)((int)std::round((float)myVertexWithAttributes.colour.red * (float)ambientColor.x));
            myScreenVertex.colour.green = (char)((int)std::round((float)myVertexWithAttributes.colour.green * (float)ambientColor.y));
            myScreenVertex.colour.blue = (char)((int)std::round((float)myVertexWithAttributes.colour.blue * (float)ambientColor.z));
            myScreenVertex.colour.alpha = myVertexWithAttributes.colour.alpha;
        }
        else
        {
            // Pass the colour through
            myScreenVertex.colour = myVertexWithAttributes.colour;
        }
        

        rasterQueue.push_back(myScreenVertex);
    }

    FakeGL::RasterisePrimitive();
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
            // Primitive can be drawn, pop a vertex
            screenVertexWithAttributes vertex = rasterQueue.front();
            rasterQueue.pop_front();
            // Now call the appropriate drawing function and return true
            RasterisePoint(vertex);
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
            ProcessFragment();
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

    // and compute the distance of each vertex from the opposing side
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

            // compute colour
            rasterFragment.colour = alpha * vertex0.colour + beta * vertex1.colour + gamma * vertex2.colour; 
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

            // delet this
            ProcessFragment();
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
    // Check if this fragment is in front
    if (depthBuffer[fragment.row][fragment.col].alpha < fragment.depth.alpha)
    {
        // This fragment is in front, update depth buffer and set in framebuffer
        depthBuffer[fragment.row][fragment.col].alpha = fragment.depth.alpha;
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
    for (auto vertex = fakeGL.vertexQueue.begin(); vertex < fakeGL.vertexQueue.end(); vertex++)
        { // per matrix
        outStream << "Vertex " << vertex - fakeGL.vertexQueue.begin() << std::endl;
        outStream << *vertex;
        } // per matrix


    outStream << "-------------------------" << std::endl;
    outStream << "Raster Queue:            " << std::endl;
    outStream << "-------------------------" << std::endl;
    for (auto vertex = fakeGL.rasterQueue.begin(); vertex < fakeGL.rasterQueue.end(); vertex++)
        { // per matrix
        outStream << "Vertex " << vertex - fakeGL.rasterQueue.begin() << std::endl;
        outStream << *vertex;
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
std::ostream &operator << (std::ostream &outStream, vertexWithAttributes &vertex)
    { // operator <<
    std::cout << "Vertex With Attributes" << std::endl;
    std::cout << "Position:   " << vertex.position << std::endl;
    std::cout << "Colour:     " << vertex.colour << std::endl;

    // you

    return outStream;
    } // operator <<

std::ostream &operator << (std::ostream &outStream, screenVertexWithAttributes &vertex) 
    { // operator <<
    std::cout << "Screen Vertex With Attributes" << std::endl;
    std::cout << "Position:   " << vertex.position << std::endl;
    std::cout << "Colour:     " << vertex.colour << std::endl;

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


    
