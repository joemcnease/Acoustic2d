#pragma once

#include <iostream>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// eunumerator to add typical ImGui components to list and then create
//

enum UI_THEME {
    CLASSIC,
    DARK,
    LIGHT
};

class UI {
    public:
        ImGuiIO context;

        // Constructor for User Interface class
        UI();
        // Destructor for deleting User Interface context
        ~UI();

        void SetStyle(UI_THEME theme);

        void SetBackend(GLFWwindow* window, const char* glsl_version);

        void LoadFont(ImGuiIO context);

        void NewFrame();

        void DefaultLayout(bool& show_window, float clear_color[], unsigned int& size_x, unsigned int& size_z, float& dx, float& dz, float& wave_scale, int& frames_per_update, float& smooth_time, int& time_scale, bool& frame_update);

        void Render();
};
