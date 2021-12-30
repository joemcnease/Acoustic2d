#include <iostream>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "ui.h"

// Constructor for User Interface class
UI::UI()
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    context = io;
    LoadFont(context);
}

// Destructor for deleting User Interface context
UI::~UI()
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void UI::SetStyle(UI_THEME theme)
{
    switch (theme) {
        case CLASSIC:
            ImGui::StyleColorsClassic();
            break;
        case DARK:
            ImGui::StyleColorsDark();
            break;
        case LIGHT:
            ImGui::StyleColorsLight();
            break;
    }
}

void UI::SetBackend(GLFWwindow* window, const char* glsl_version)
{
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);
}

void UI::LoadFont(ImGuiIO io)
{
    io.Fonts->AddFontDefault();
}

void UI::NewFrame()
{
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void UI::DefaultLayout(bool& show_window, float clear_color[], unsigned int& size_x, unsigned int& size_z, float& dx, float& dz, float& wave_scale, int& frames_per_update, float& smooth_time, int& time_scale, bool& frame_update)
{
    // Basic window
    {
        // const ImU32 grid_x_min = 0, grid_x_max = 1000;
        // const ImU32 grid_z_min = 0, grid_z_max = 1000;

        ImGui::Begin("OpenViz Editor");
        ImGui::Checkbox("Another Window", &show_window);
        ImGui::ColorEdit3("Background Color", clear_color);
        ImGui::SliderInt("Frames Per Wave Update", &frames_per_update, 1, 60);
        ImGui::Checkbox("Update Wave (time-step)", &frame_update);
        ImGui::SliderFloat("dx", &dx, 0, 10);
        ImGui::SliderFloat("dz", &dz, 0, 10);
        ImGui::SliderFloat("Wave Scale", &wave_scale, -1000, 1000);
        ImGui::SliderInt("Time Scale", &time_scale, 1, 1000);
        ImGui::Text("MS/FRAME: %.3f -- FPS: %.1f FPS", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        ImGui::End();
    }

    if (show_window) {
        ImGui::Begin("Another Window", &show_window);
        ImGui::Text("Hello from another window!");
        if (ImGui::Button("Close Me"))
            show_window = false;
        ImGui::End();
    }
}

void UI::Render() {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
