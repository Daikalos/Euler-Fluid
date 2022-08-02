#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/OpenGL.hpp>

#include <execution>

#include "Config.h"

struct Vertex
{
	GLfloat x, y;
};

struct Color
{
	GLfloat r, g, b;
};

class Fluid
{
public:
	Fluid(Config* config, const size_t& width, const size_t& height, const float& diff, const float& visc);

	inline int IX(int x, int y)
	{
		return x + y * W;
	}

	inline int IX(int i)
	{
		return (i % W) + i;
	}

	inline int safe_IX(int x, int y)
	{
		if (x < 0) x = 0;
		else if (x > W - 1) x = W - 1;

		if (y < 0) y = 0;
		else if (y > H - 1) y = H - 1;
		
		return IX(x, y);
	}

	void add_density(int x, int y, float amount)
	{
		density[safe_IX(x, y)] += amount;
	}

	void add_velocity(int x, int y, float vx, float vy)
	{
		int index = safe_IX(x, y);

		this->vx[index] += vx;
		this->vy[index] += vy;
	}

	void step_line(int x0, int y0, int x1, int y1, int dx, int dy, float a);

	void update(const float& dt);
	void draw();

private:
	void lin_solve(int b, float* x, float* x0, float a, float c);

	void set_bnd(int b, float* x);

	void diffuse(int b, float* x, float* x0, float dt);
	void project(float* u, float* v, float* p, float* div);
	void advect(int b, float* d, float* d0, float* u, float* v, float dt);

	void fade_density();

	inline float map_to_range(float val, float minIn, float maxIn, float minOut, float maxOut)
	{
		float x = (val - minIn) / (maxIn - minIn);
		return minOut + (maxOut - minOut) * x;
	}

private:
	Config* config;

	size_t W, H, N, V;
	float diff, visc;

	std::vector<float> vx;
	std::vector<float> vy;
	std::vector<float> vx_prev;
	std::vector<float> vy_prev;

	std::vector<float> density;
	std::vector<float> density_prev;

	std::vector<int> range;

	std::vector<Vertex> vertices;
	std::vector<Color> colors;
};

