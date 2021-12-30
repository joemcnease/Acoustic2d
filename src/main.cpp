#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdlib.h>

#include <glm/ext/matrix_clip_space.hpp>
#include <glm/ext/matrix_float4x3.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/ext/vector_float3.hpp>
#include <glm/geometric.hpp>
#include <glm/trigonometric.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "stb_image.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "shader.h"
#include "camera.h"
#include "ui.h"
#include "vao.h"
#include "vbo.h"
#include "ebo.h"


// Function declarations
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

void generate_axes(std::vector<float>& vertices, float range);
void generate_surface(std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int size_x, unsigned int size_z, float dx, float dz, float& min, float& max);
void generate_surface_from_func(float (*func)(float, float), std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int size_x, unsigned int size_z, float dx, float dz);
float surf(float x, float z, float t);
float zero_surf(float x, float z);
void time_step(float (*func)(float, float, float), std::vector<float>& vertices, unsigned int size_x, unsigned int size_z, float dx, float dz, float& min, float& max, float t, float wave_scale);
void wave_time_step(std::vector<double>& p, std::vector<double>& pOld, std::vector<double>& pNew, std::vector<double>& d2px, std::vector<double>& d2pz, std::vector<double>& c, std::vector<double>& src, unsigned int nx, unsigned int nz, unsigned int xs, unsigned int zs, double dx, double dz, double dt, unsigned int it);
void update_wave_surf(std::vector<double>& p, std::vector<float>& vertices, unsigned int size_x, unsigned int size_z, float dx, float dz, float& min, float& max, float wave_scale);
void generate_sphere(std::vector<float>* vertices, std::vector<float>* indices, float r, int stackCount, int sectorCount);
void printVerticesAndIndices(std::vector<float>& vertices, std::vector<unsigned int>& indices);
void zeros(std::vector<float>& vertices);
void source_func(std::vector<double>& source, const std::vector<double>& time, const double t0, const double dt, const double T);
void diff(std::vector<double>& vec, const unsigned int order, bool retainSize = false);
void reset_wave_surf(std::vector<double>& p, std::vector<double>& pOld, std::vector<double>& pNew, std::vector<double>& d2px, std::vector<double>& d2pz, std::vector<double>& src);

// Window Settings
const unsigned int SCR_WIDTH = 1200;
const unsigned int SCR_HEIGHT = 800;

// ImGui edit mode
bool editMode = false;

// Reset wave surface
bool reset_wave = false;
bool frame_update = false;

// Camera
Camera camera(glm::vec3(0.0f, 650.0f, 650.0f), glm::vec3(0.0f, 1.0f, 0.0f), -90.0f, -45.0f);
// Camera camera(glm::vec3(0.0f, 800.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), -90.0f, -89.0f);
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// Global timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    // glfw window creation
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // setup Dear ImGui context
    UI ui;
    ui.SetStyle(CLASSIC);
    ui.SetBackend(window, "#version 330");

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // enable global opengl state
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);

    // build and compile our shader program
    // ------------------------------------
    Shader axesShader("src/shader.vs", "src/shader.fs");
    Shader surfaceShader("src/surfaceShader.vs", "src/surfaceShader.fs");

    // init graph axes and surface
    // -------------------------------------------------------------------------------------
    std::vector<float> axesVertices;
    generate_axes(axesVertices, 1000.0f);
    VAO axesVAO;
    VBO axesVBO(axesVertices.data(), axesVertices.size() * sizeof(float), GL_STATIC_DRAW);
    axesVAO.VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    
    // grids allocated for simulation
    unsigned int size_x = 400;
    unsigned int size_z = 400;
    unsigned int nt = 2000;
    double dx = 10.0;
    double dz = 10.0;
    double dt = 0.001;
    double c0 = 3000.0;
    double f0 = 100.0;
    double t0 = 100.0*dt;
    double T = 1.0 / f0;

    // source location
    unsigned int xs = 10; // std::floor(size_x/2) - 1;
    unsigned int zs = std::floor(size_z/2) - 1;

    std::vector<double> p(size_z * size_x, 0);
    std::vector<double> pOld = p;
    std::vector<double> pNew = p;
    std::vector<double> d2dx = p;
    std::vector<double> d2dz = p;
    std::vector<double> c(p.size(), c0);

    // add some layers to the model
    for (unsigned int i = 0; i < size_z; ++i) {
        for (unsigned int j = 0; j < size_x; ++j) {
            unsigned int index = i*size_x + j;
            if ((std::pow((i-std::floor(size_x/2)), 2) + std::pow((j-std::floor(size_z/2)), 2)) < 1000) {
                c[index] = 2000.0;
            }
        }
    }

    double eps = c0 * dt / dx;
    std::cout << "Stability criterion: " << eps << std::endl;
    std::cout << "Source frequency: " << f0 << " Hz" << std::endl;

    // number time steps and time vector
    std::vector<double> t(nt);
    for (unsigned int i=0; i<nt; ++i) {
        t[i] = i * dt;
    }
    
    // source time function
    std::vector<double> src(t.size());
    source_func(src, t, t0, dt, T);
    diff(src, 1);
    // for (unsigned int i = 0; i < src.size(); ++i) {
    //     src[i] = src[i] / dt;
    //     std::cout << "src[i] = " << src[i] << "\n";
    // }
    // std::exit(0);

    // get min and max of velocity
    // const auto minmax = std::minmax_element(std::begin(c), std::end(c));
    // const float cmin = *minmax.first;
    // const float cmax = *minmax.second;

    std::vector<float> surfaceVertices;
    std::vector<unsigned int> surfaceIndices;
    unsigned int vsize_x = size_x - 1;
    unsigned int vsize_z = size_z - 1;
    float vdx = 1.0f;
    float vdz = 1.0f;
    float min;
    float max;
    generate_surface_from_func(&zero_surf, surfaceVertices, surfaceIndices, vsize_x, vsize_z, vdx, vdz);
    VAO surfaceVAO;
    VBO surfaceVBO(surfaceVertices.data(), surfaceVertices.size() * sizeof(float), GL_STATIC_DRAW);
    surfaceVAO.VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    EBO surfaceEBO(surfaceIndices.data(), surfaceIndices.size() * sizeof(unsigned int), GL_STATIC_DRAW);



    bool show_window = false;
    float clear_color[4] = {0.45f, 0.55f, 0.60f, 1.00f};
    int frames_per_update = 60;
    int frames_since_update = 0;
    float wave_scale = 1.0f;
    float smooth_time = 0.0f;
    int time_scale = 1;

    unsigned int time_step = 0;

    // render loop
    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // poll events events
        glfwPollEvents();

        // start the dear imgui frame
        ui.NewFrame();

        // Keyboard input
        processInput(window);

        // Rendering commands
        glClearColor(clear_color[0], clear_color[1], clear_color[2], clear_color[3]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        // Initialize MVP matrices
        glm::mat4 model = glm::mat4(1.0f);
        glm::mat4 view = camera.GetViewMatrix();
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_HEIGHT / (float)SCR_WIDTH, 0.1f, 2000.0f);

        // activate shader
        axesShader.use();
        axesVAO.Bind();
        model = glm::translate(model, glm::vec3(0.0f, 0.0f, 0.0f));
        axesShader.setMat4("model", model);
        axesShader.setMat4("projection", projection);
        axesShader.setMat4("view", view);
        axesShader.setVec3("ourColor", glm::vec3(1.0f, 1.0f, 1.0f));
        glDrawArrays(GL_LINES, 0, 6);

        // axesShader.use();
        // gridVAO.Bind();
        // model = glm::mat4(1.0f);
        // model = glm::translate(model, glm::vec3(0.0f, 0.0f, 0.0f));
        // axesShader.setMat4("model", model);
        // axesShader.setMat4("projection", projection);
        // axesShader.setMat4("view", view);
        // axesShader.setVec3("ourColor", glm::vec3(1.0f, 1.0f, 1.0f));
        // glDrawElements(GL_TRIANGLES, gridIndices.size(), GL_UNSIGNED_INT, 0);

        surfaceShader.use();
        surfaceVAO.Bind(); surfaceEBO.Bind();
        model = glm::mat4(1.0f);
        model = glm::translate(model, glm::vec3(0.0f, 0.0f, 0.0f));
        surfaceShader.setMat4("model", model);
        surfaceShader.setMat4("projection", projection);
        surfaceShader.setMat4("view", view);
        surfaceShader.setFloat("min", min);
        surfaceShader.setFloat("max", max);
        glDrawElements(GL_TRIANGLES, surfaceIndices.size(), GL_UNSIGNED_INT, 0);

        frames_since_update++;
        smooth_time += 0.000001*time_scale;
        if (frame_update && frames_since_update > frames_per_update)
        {
            wave_time_step(p, pOld, pNew, d2dx, d2dz, c, src, size_x, size_z, xs, zs, dx, dz, dt, time_step); 
            update_wave_surf(p, surfaceVertices, vsize_x, vsize_z, vdx, vdz, min, max, wave_scale);
            surfaceVBO.SubBuffer(surfaceVertices.data(), surfaceVertices.size() * sizeof(float), 0);
            frames_since_update = 0;
            if (time_step < nt) {
                time_step++;
            }
            else {time_step = 0;}
        }

        if (reset_wave) {
            reset_wave_surf(p, pOld, pNew, d2dx, d2dz, src);
            time_step = 0;
            reset_wave = false;
        }

        // User interface layout
        if (editMode)
            ui.DefaultLayout(show_window, clear_color, size_x, size_z, vdx, vdz, wave_scale, frames_per_update, smooth_time, time_scale, frame_update);

        // Rendering
        ui.Render();

        glfwSwapBuffers(window);
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}


void processInput(GLFWwindow *window)
{
    // escape key
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, true);
    }

    // camera motion with wasd
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    if (editMode == false && glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        editMode = true;
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
    if (editMode == true && glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        editMode = false;
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    }
    if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) {
        reset_wave = true;
    }

}


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}


void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;

    lastX = xpos;
    lastY = ypos;

    if (!editMode) {
        camera.ProcessMouseMovement(xoffset, yoffset);
    }
}


void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}


void generate_axes(std::vector<float>& vertices, float range)
{
        vertices.push_back(-range);
        vertices.push_back(0);
        vertices.push_back(0);

        vertices.push_back(range);
        vertices.push_back(0);
        vertices.push_back(0);

        vertices.push_back(0);
        vertices.push_back(-range);
        vertices.push_back(0);

        vertices.push_back(0);
        vertices.push_back(range);
        vertices.push_back(0);

        vertices.push_back(0);
        vertices.push_back(0);
        vertices.push_back(-range);

        vertices.push_back(0);
        vertices.push_back(0);
        vertices.push_back(range);
}


void generate_surface(std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int size_x, unsigned int size_z, float dx, float dz, float& min, float& max)
{
    float xmax = std::floor(size_x/2)*dx;
    float zmax = std::floor(size_z/2)*dz;

    float x, y, z;

    for (unsigned int i = 0; i < (size_x + 1); ++i)
    {
        for (unsigned int j = 0; j < (size_z + 1); ++j)
        {
            x = -xmax + i*dx;
            z = -zmax + j*dz;
            y = std::sin(std::sqrt(x*x + z*z));

            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);

            if (i < size_x && j < size_z)
            {
                indices.push_back(i*(size_z+1) + j);
                indices.push_back(i*(size_z+1) + (size_z+1) + j);
                indices.push_back(i*(size_z+1) + (size_z+1) + j + 1);
                indices.push_back(i*(size_z+1) + j);
                indices.push_back(i*(size_z+1) + j + 1);
                indices.push_back(i*(size_z+1) + (size_z+1) + j + 1);
            }

            if (i == 0 && j == 0) {min = y; max = y;}
            if (y < min) {min = y;}
            if (y > max) {max = y;}
        }
    }
}


float surf(float x, float z, float t)
{
    return std::sin(x * t) + std::cos(z * t);
}


float zero_surf(float x, float z)
{
    return 0;
}


void zeros(std::vector<float>& vertices)
{
    for (unsigned int i=0; i<vertices.size()/3; i+=3)
    {
        vertices[i] = 0;
        vertices[i+1] = 0;
        vertices[i+2] = 0;
    }
}


void source_func(std::vector<double>& source, const std::vector<double>& time, const double t0, const double dt, const double T) {
    for (unsigned int i=0; i<source.size(); ++i) {
        source[i] = std::exp(-1.0f / (T * T) * ((time[i] - t0) * (time[i] - t0)));
    }
}


void diff(std::vector<double>& vec, const unsigned int order, bool retainSize) {
    for (unsigned int i=0; i<order; ++i) {
        for (unsigned int i=0; i<vec.size()-1; ++i) {
            vec[i] = vec[i+1] - vec[i];
        }
    }

    if (retainSize) {
        for (unsigned int i = 0; i < order; ++i) {
            vec.push_back(vec[vec.size() - 1]);
        }
    }
}


void wave_time_step(std::vector<double>& p, std::vector<double>& pOld, std::vector<double>& pNew, std::vector<double>& d2px, std::vector<double>& d2pz, std::vector<double>& c, std::vector<double>& src, unsigned int nx, unsigned int nz, unsigned int xs, unsigned int zs, double dx, double dz, double dt, unsigned int it) {
    // update pressure field based on 2d finite difference solver
    for (unsigned int j = 0; j < nz - 1; ++j) {
        for (unsigned int k = 1; k < nx - 2; ++k) {
            unsigned int index = j*nx + k;
            d2px[index] = (p[index - 1] - 2*p[index] + p[index + 1]) / (dx * dx);
        }
    }
    for (unsigned int j = 1; j < nz - 2; ++j) {
        for (unsigned int k = 0; k < nx - 1; ++k) {
            unsigned int index = j*nx + k;
            d2pz[index] = (p[index - nx] - 2*p[index] + p[index + nx]) / (dz * dz);
        }
    }

    for (unsigned int i = 0; i < pNew.size() - 1; ++i) {
        pNew[i] = 2*p[i] - pOld[i] + (dt * dt) * (c[i]*c[i]) * (d2px[i] + d2pz[i]);
    }
    pNew[zs*nx + xs] = pNew[zs*nx + xs] + src[it];
    pOld = p;
    p = pNew;
}


void printVerticesAndIndices(std::vector<float>& vertices, std::vector<unsigned int>& indices)
{
    std::cout << "Vertices" << std::endl;
    for (unsigned int i = 0; i < vertices.size()/3; i += 3)
    {
        std::printf("(x, y, z) = (%.3f, %.3f, %.3f)\n", vertices[i], vertices[i+1], vertices[i+2]);
    }
    std::cout << "Indices" << std::endl;
    for (unsigned int i = 0; i < indices.size()/6; i += 6)
    {
        std::printf("(index1, index2, index3, index4, index5, index6) = (%3i, %3i, %3i, %3i, %3i, %3i)\n", indices[i], indices[i+1], indices[i+2], indices[i+3], indices[i+4], indices[i+5]);
    }

}


void generate_surface_from_func(float (*func)(float, float), std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int size_x, unsigned int size_z, float dx, float dz)
{
    float xmax = std::floor(size_x/2)*dx;
    float zmax = std::floor(size_z/2)*dz;

    float x, y, z;

    for (unsigned int i = 0; i < (size_x + 1); ++i)
    {
        for (unsigned int j = 0; j < (size_z + 1); ++j)
        {
            x = -xmax + i*dx;
            z = -zmax + j*dz;
            y = func(x, z);

            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);

            if (i < size_x && j < size_z)
            {
                indices.push_back(i*(size_z+1) + j);
                indices.push_back(i*(size_z+1) + (size_z+1) + j);
                indices.push_back(i*(size_z+1) + (size_z+1) + j + 1);
                indices.push_back(i*(size_z+1) + j);
                indices.push_back(i*(size_z+1) + j + 1);
                indices.push_back(i*(size_z+1) + (size_z+1) + j + 1);
            }
        }
    }
}

void update_wave_surf(std::vector<double>& p, std::vector<float>& vertices, unsigned int size_x, unsigned int size_z, float dx, float dz, float& min, float& max, float wave_scale)
{
    float xmax = std::floor(size_x/2)*dx;
    float zmax = std::floor(size_z/2)*dz;

    float x, y, z;

    int index = 0;
    int wave_index = 0;
    for (unsigned int i = 0; i < (size_x + 1); ++i)
    {
        for (unsigned int j = 0; j < (size_z + 1); ++j)
        {
            x = -xmax + i*dx;
            z = -zmax + j*dz;
            y = wave_scale * (float)p[wave_index];

            vertices[index] = x;
            vertices[index + 1] = y;
            vertices[index + 2] = z;

            if (i == 0 && j == 0) {min = y; max = y;}
            if (y < min) {min = y;}
            if (y > max) {max = y;}

            index += 3;
            wave_index += 1;
        }
    }
}

void reset_wave_surf(std::vector<double>& p, std::vector<double>& pOld, std::vector<double>& pNew, std::vector<double>& d2px, std::vector<double>& d2pz, std::vector<double>& src) {
    std::fill(p.begin(), p.end(), 0.0);
    std::fill(pOld.begin(), pOld.end(), 0.0);
    std::fill(pNew.begin(), pNew.end(), 0.0);
    std::fill(d2px.begin(), d2px.end(), 0.0);
    std::fill(d2pz.begin(), d2pz.end(), 0.0);
}

void time_step(float (*func)(float, float, float), std::vector<float>& vertices, unsigned int size_x, unsigned int size_z, float dx, float dz, float& min, float& max, float t, float wave_scale)
{
    float xmax = std::floor(size_x/2)*dx;
    float zmax = std::floor(size_z/2)*dz;

    float x, y, z;

    int index = 0;
    for (unsigned int i = 0; i < (size_x + 1); ++i)
    {
        for (unsigned int j = 0; j < (size_z + 1); ++j)
        {
            x = -xmax + i*dx;
            z = -zmax + j*dz;
            y = wave_scale * func(x, z, t);

            vertices[index] = x;
            vertices[index + 1] = y;
            vertices[index + 2] = z;

            if (i == 0 && j == 0) {min = y; max = y;}
            if (y < min) {min = y;}
            if (y > max) {max = y;}

            index += 3;
        }
    }
}
