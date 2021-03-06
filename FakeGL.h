//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  ------------------------
//  FakeGL.h
//  ------------------------
//  
//  A unit for implementing OpenGL workalike calls
//  
///////////////////////////////////////////////////

// include guard
#ifndef FAKEGL_H
#define FAKEGL_H

#include "Cartesian3.h"
#include "Homogeneous4.h"
#include "Matrix4.h"
#include "RGBAImage.h"
#include <vector>
#include <deque>
#include <stack>

// we will store all of the FakeGL context in a class object
// this is similar to the real OpenGL which handles multiple windows
// by having separate variables for each stored in separate graphics
// contexts.   
//
// We have simplified the calls, so the code is not *identical* to
// the OpenGL calls, but it's pretty close

// class constants
// constants for Begin()
const unsigned int FAKEGL_POINTS = 0;
const unsigned int FAKEGL_LINES = 1;
const unsigned int FAKEGL_TRIANGLES = 2;
// bitflag constants for Clear()
const unsigned int FAKEGL_COLOR_BUFFER_BIT = 1;
const unsigned int FAKEGL_DEPTH_BUFFER_BIT = 2;
// constants for Enable()/Disable()
const unsigned int FAKEGL_LIGHTING = 1;
const unsigned int FAKEGL_TEXTURE_2D = 2;
const unsigned int FAKEGL_DEPTH_TEST = 3;
const unsigned int FAKEGL_PHONG_SHADING = 4;
// constants for Light() - actually bit flags
const unsigned int FAKEGL_POSITION = 1;
const unsigned int FAKEGL_AMBIENT = 2;
const unsigned int FAKEGL_DIFFUSE = 4;
const unsigned int FAKEGL_AMBIENT_AND_DIFFUSE = 6;
const unsigned int FAKEGL_SPECULAR = 8;
// additional constants for Material()
const unsigned int FAKEGL_EMISSION = 16;
const unsigned int FAKEGL_SHININESS = 32;
// constants for matrix operations
const unsigned int FAKEGL_MODELVIEW = 1;
const unsigned int FAKEGL_PROJECTION = 2;
// constants for texture operations
const unsigned int FAKEGL_MODULATE = 1;
const unsigned int FAKEGL_REPLACE = 2;

// class with vertex attributes
class vertexWithAttributes
{ // class vertexWithAttributes
    public:
    // Position in OCS
    Homogeneous4 position;
    
    // Colour
    RGBAValue colour;
    
    // Normal
    Homogeneous4 normal;
    
    // Material properties
    // Reflectance
    Homogeneous4 ambientReflectance;
    Homogeneous4 diffuseReflectance;
    Homogeneous4 specularReflectance;
    // Additives
    Homogeneous4 emission;
    float shininess = 0.0f;

    // Texture coordinates
    float u = 0.0f;
    float v = 0.0f;

}; // class vertexWithAttributes

// class for a vertex after transformation to screen space
class screenVertexWithAttributes
{ // class screenVertexWithAttributes
    public:
    // Position in DCS
    Cartesian3 position;
    
    // Colour
    RGBAValue colour;
    
    // Normal
    Homogeneous4 normal;

    // Material properties
    // Reflectance
    Homogeneous4 ambientReflectance;
    Homogeneous4 diffuseReflectance;
    Homogeneous4 specularReflectance;
    // Additives
    Homogeneous4 emission;
    float shininess = 0.0f;

    // Texture coordinates
    float u = 0.0f;
    float v = 0.0f;

}; // class screenVertexWithAttributes

// class for a fragment with attributes
class fragmentWithAttributes
{ // class fragmentWithAttributes
    public:
    // the row & column address in the framebuffer
    int row, col;
    // the RGBA colour of the fragment
    RGBAValue colour;
    // the depth value of the fragment
    RGBAValue depth;

}; // class fragmentWithAttributes

// the class storing the FakeGL context
class FakeGL
{ // class FakeGL
    // for the sake of simplicity & clarity, we make everything public
    public:
    //-----------------------------
    // MATRIX STATE
    //-----------------------------
    
    // Matrix mode- initial value in OpenGL is GL_MODELVIEW
    unsigned int matrixMode = FAKEGL_MODELVIEW;
    // Matrix stacks
    std::stack<Matrix4> modelViewMatStack;
    std::stack<Matrix4> projectionMatStack;
    // The current matrices -default to identity
    Matrix4 modelViewMat;
    Matrix4 projectionMat;
    
    //-----------------------------
    // ATTRIBUTE STATE
    //-----------------------------

    // Store the current colour as set by Color3f (defaults to white, alpha = 1, := 255 in this case)
    RGBAValue currentColor;
    // The current normals as set by Normal3f
    Homogeneous4 currentNormal;
    // The current material properties, as set by Materialfv (or Materialf in case of shininess)
    // Reflectance
    Homogeneous4 currentAmbientReflectance;
    Homogeneous4 currentDiffuseReflectance;
    Homogeneous4 currentSpecularReflectance;
    // Additives
    Homogeneous4 currentEmission;
    float currentShininess;
    //Texture coords
    float currentU = 0.0;
    float currentV = 0.0;

    //-----------------------------
    // OUTPUT FROM INPUT STAGE
    // INPUT TO TRANSFORM STAGE
    //-----------------------------

    // we want a queue of vertices with attributes for passing to the rasteriser
    std::deque<vertexWithAttributes> vertexQueue;

    //-----------------------------
    // TRANSFORM/LIGHTING STATE
    //-----------------------------

    // Light position
    Homogeneous4 lightPos;
    // Light colour. Using Homogenous4 so they can be stored as floats
    Homogeneous4 ambientColor;
    Homogeneous4 diffuseColor;
    Homogeneous4 specularColor;

    // Projection
    float zNear = -1.0f;
    float zFar = 1.0f;

    //-----------------------------
    // OUTPUT FROM TRANSFORM STAGE
    // INPUT TO RASTER STAGE
    //-----------------------------
    std::deque<screenVertexWithAttributes> rasterQueue;

    //-----------------------------
    // RASTERISE STATE
    //-----------------------------

    // Primitive
    // Results where a primitive hasn't been set are undefined in OpenGL. I will refuse to rasterise
    //      anything without a primitive being set. Default is -1 as 0/null is actually FAKEGL_POINTS
    int primitiveMode = -1;

    // In OpenGL, only certian calls are permitted between glBegin and glEnd, to be completely accurate,
    //      blocking some functions (could) be implemented here, so track the state
    bool primitiveAssembly = false;

    // Primitive attributes- defaults are 1
    float pointSize = 1.0;
    float lineWidth = 1.0;

    // Lighting flags
    bool lighting = false;
    bool phongShading = false;

    //-----------------------------
    // OUTPUT FROM RASTER STAGE
    // INPUT TO FRAGMENT STAGE
    //-----------------------------
    std::deque<fragmentWithAttributes> fragmentQueue;

    //-----------------------------
    // TEXTURE STATE
    //-----------------------------

    // Enabled/disabled flag
    bool textureEnabled = false;
    // Pointer to the current texture
    const RGBAImage* currentTexture;
    // Texture mode
    int textureMode = FAKEGL_REPLACE;

    //-----------------------------
    // FRAMEBUFFER STATE
    //-----------------------------
    
    // Values to use for clearing colour and depth buffers
    RGBAValue clearColorVal;
    RGBAValue depthVal;
    // Depth test on/off
    bool depthTest = false;

    // Viewport state- x, y specify the lower left corner of the viewport, and default to 0. 
    //      Width/height are set by the window system, and have no default.
    int viewportX = 0;
    int viewportY = 0;
    unsigned int viewportWidth;
    unsigned int viewportHeight;

    //-----------------------------
    // OUTPUT FROM FRAGMENT STAGE
    //-----------------------------
    
    // the frame buffer itself
    RGBAImage frameBuffer;
     
    // rather than define an extra class, we will cheat and use 
    // a second RGBAImage in which the alpha stores the depth buffer
    RGBAImage depthBuffer;
    
    //-------------------------------------------------//
    //                                                 //
    // CONSTRUCTOR / DESTRUCTOR                        //
    //                                                 //
    //-------------------------------------------------//
    
    // constructor
    FakeGL();
    
    // destructor
    ~FakeGL();
    
    //-------------------------------------------------//
    //                                                 //
    // GEOMETRIC PRIMITIVE ROUTINES                    //
    //                                                 //
    //-------------------------------------------------//
    
    // starts a sequence of geometric primitives
    void Begin(unsigned int PrimitiveType);
    
    // ends a sequence of geometric primitives
    void End();

    // sets the size of a point for drawing
    void PointSize(float size);

    // sets the width of a line for drawing purposes
    void LineWidth(float width);
    
    //-------------------------------------------------//
    //                                                 //
    // MATRIX MANIPULATION ROUTINES                    //
    //                                                 //
    //-------------------------------------------------//

    // set the matrix mode (i.e. which one we change)   
    void MatrixMode(unsigned int whichMatrix);

    // pushes a matrix on the stack
    void PushMatrix();

    // pops a matrix off the stack
    void PopMatrix();

    // load the identity matrix
    void LoadIdentity();
    
    // multiply by a known matrix in column-major format
    void MultMatrixf(const float *columnMajorCoordinates);

    // sets a perspective projection matrix
    void Frustum(float left, float right, float bottom, float top, float zNear, float zFar);

    // sets an orthographic projection matrix
    void Ortho(float left, float right, float bottom, float top, float zNear, float zFar);

    // rotate the matrix
    void Rotatef(float angle, float axisX, float axisY, float axisZ);

    // scale the matrix
    void Scalef(float xScale, float yScale, float zScale);
    
    // translate the matrix
    void Translatef(float xTranslate, float yTranslate, float zTranslate);
    
    // sets the viewport
    void Viewport(int x, int y, int width, int height);

    //-------------------------------------------------//
    //                                                 //
    // VERTEX ATTRIBUTE ROUTINES                       //
    //                                                 //
    //-------------------------------------------------//

    // sets colour with floating point
    void Color3f(float red, float green, float blue);
    
    // sets material properties
    void Materialf(unsigned int parameterName, const float parameterValue);
    void Materialfv(unsigned int parameterName, const float *parameterValues);

    // sets the normal vector
    void Normal3f(float x, float y, float z);

    // sets the texture coordinates
    void TexCoord2f(float u, float v);

    // sets the vertex & launches it down the pipeline
    void Vertex3f(float x, float y, float z);

    //-------------------------------------------------//
    //                                                 //
    // STATE VARIABLE ROUTINES                         //
    //                                                 //
    //-------------------------------------------------//

    // disables a specific flag in the library
    void Disable(unsigned int property);
    
    // enables a specific flag in the library
    void Enable(unsigned int property);
    
    //-------------------------------------------------//
    //                                                 //
    // LIGHTING STATE ROUTINES                         //
    //                                                 //
    //-------------------------------------------------//

    // sets properties for the one and only light
    void Light(int parameterName, const float *parameterValues);

    //-------------------------------------------------//
    //                                                 //
    // TEXTURE PROCESSING ROUTINES                     //
    //                                                 //
    // Note that we only allow one texture             //
    // so glGenTexture & glBindTexture aren't needed   //
    //                                                 //
    //-------------------------------------------------//

    // sets whether textures replace or modulate
    void TexEnvMode(unsigned int textureMode);

    // sets the texture image that corresponds to a given ID
    void TexImage2D(const RGBAImage &textureImage);

    //-------------------------------------------------//
    //                                                 //
    // FRAME BUFFER ROUTINES                           //
    //                                                 //
    //-------------------------------------------------//

    // clears the frame buffer
    void Clear(unsigned int mask);
    
    // sets the clear colour for the frame buffer
    void ClearColor(float red, float green, float blue, float alpha);
    
    //-------------------------------------------------//
    //                                                 //
    // ROUTINE TO FLUSH THE PIPELINE                   //
    //                                                 //
    //-------------------------------------------------//
    
    // flushes the pipeline
    void Flush();

    //-------------------------------------------------//
    //                                                 //
    // MAJOR PROCESSING ROUTINES                       //
    //                                                 //
    //-------------------------------------------------//

    // transform one vertex & shift to the transformed queue
    void TransformVertex();

    // rasterise a single primitive if there are enough vertices on the queue
    bool RasterisePrimitive();

    // rasterises a single point
    void RasterisePoint(screenVertexWithAttributes &vertex0);

    // rasterises a single line segment
    void RasteriseLineSegment(screenVertexWithAttributes &vertex0, screenVertexWithAttributes &vertex1);
    
    // rasterises a single triangle
    void RasteriseTriangle(screenVertexWithAttributes &vertex0, screenVertexWithAttributes &vertex1, screenVertexWithAttributes &vertex2);
    
    // process a single fragment
    void ProcessFragment();
    
}; // class FakeGL

// standard routine for dumping the entire FakeGL context (except for texture / image)
std::ostream &operator << (std::ostream &outStream, FakeGL &fakeGL); 

// subroutines for other classes
std::ostream &operator << (std::ostream &outStream, vertexWithAttributes &vertex); 
std::ostream &operator << (std::ostream &outStream, screenVertexWithAttributes &vertex); 
std::ostream &operator << (std::ostream &outStream, fragmentWithAttributes &fragment); 

// include guard        
#endif