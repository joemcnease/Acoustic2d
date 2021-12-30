#version 330 core

layout (location = 0) in vec3 aPos;

varying vec4 surfaceColor;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float max;
uniform float min;

vec4 RED = vec4(1.0f, 0.0f, 0.0f, 1.0f);
vec4 BLUE = vec4(0.0f, 0.0f, 1.0f, 1.0f);

void main() {
    gl_Position = projection * view * model * vec4(aPos, 1.0);

    float r = (aPos.y - min)/(max - min);
    // if      (r <= 1.0 && r > 0.8) {
    //     surfaceColor = vec4(1.0f, 0.0f, 0.0f, 1.0f);
    // }
    // else if (r <= 0.8 && r > 0.6) {
    //     surfaceColor = vec4(1.0f, 1.0f, 0.0f, 1.0f);
    // }
    // else if (r <= 0.6 && r > 0.4) {
    //     surfaceColor = vec4(0.0f, 1.0f, 0.0f, 1.0f);
    // }
    // else if (r <= 0.4 && r > 0.2) {
    //     surfaceColor = vec4(0.0f, 1.0f, 1.0f, 1.0f);
    // }
    // else if (r <= 0.2 && r >= 0.0) {
    //     surfaceColor = vec4(0.0f, 0.0f, 1.0f, 1.0f);
    // }
    // else {
    //     surfaceColor = vec4(r, 1.0, r, 1.0f);
    // }
    surfaceColor = mix(BLUE, RED, r);
};
