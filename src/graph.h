#pragma once

#include <map>
#include <glad/glad.h>
#include <glm/glm.hpp>

#include "shader.h"
#include "vao.h"
#include "vbo.h"
#include "ebo.h"


// The graph class encapsulates all the data needed to plot a line or surface
// plot with axes, tick marks, and a grid.
class Graph
{
    private:
        Shader m_shader;

        VAO axesVAO;
        VBO axesVBO;

        VAO axesTicksVAO;
        VBO axesTicksVBO;

        VAO gridVAO;
        VBO gridVBO;
        EBO gridEBO;

        VAO surfaceVAO;
        VBO surfaceVBO;
        EBO surfaceEBO;

    public:
        Graph();
        ~Graph();

        void line();
        void surface();
        void getColorSchemes();
};
