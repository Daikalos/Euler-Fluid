#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/OpenGL.hpp>

#include <thread>
#include <execution>

#include "Config.h"
#include "ThreadPool.h"
#include "Utilities.h"

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
	Fluid(Config* config, const size_t& width, const size_t& height, const float diff, const float visc);

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

	void update(const float dt);
	void draw();

private:
	[[nodiscard]] constexpr int IX(const int x, const int y) const noexcept { return x + y * W; }
	[[nodiscard]] constexpr int IX(const int i) const noexcept { return (i % W) + i; }

	[[nodiscard]] constexpr int safe_IX(int x, int y) const noexcept
	{
		if (x < 0) x = 0;
		else if (x > W - 1) x = W - 1;

		if (y < 0) y = 0;
		else if (y > H - 1) y = H - 1;

		return IX(x, y);
	}

private:
	void lin_solve(float* x, const float* x0, const float a, const int b, const float c);

	void set_bnd(float* x, const int b);

	void diffuse(float* x, const float* x0, const int b, const float dt);
	void advect(float* d, const float* d0, const float* vx, const float* vy, const int b, const float dt);
	void project(float* u, float* v, float* p, float* div);

	void fade_density();

private:
	Config* config;

	const size_t W, H, N, V;
	const float diff, visc;

	float* vx;
	float* vy;
	float* vx_prev;
	float* vy_prev;

	float* density;
	float* density_prev;

	Vertex* vertices;
	Color* colors;

	ThreadPool threadpool; // using threadpooling improves the performance, but currently looking if the loops can be parallized
};

