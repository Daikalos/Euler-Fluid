#include "Fluid.h"

Fluid::Fluid(Config* config, const size_t& width, const size_t& height, const float& diff, const float& visc)
	: config(config), W(width), H(height), N(width * height), V(N * 4), diff(diff), visc(visc)
{
	vx.resize(N);
	vy.resize(N);

	vx_prev.resize(N);
	vy_prev.resize(N);

	density.resize(N);
	density_prev.resize(N);

	range.reserve(N);
	for (int i = 0; i < N; ++i)
		range.push_back(i);

	vertices.resize(V);
	colors.resize(V);

	for (int y = 0; y < H; ++y)
		for (int x = 0; x < W; ++x)
		{
			int i = IX(x, y) * 4;

			int x0 = x * config->scale;
			int x1 = (x + 1) * config->scale;
			int y0 = y * config->scale;
			int y1 = (y + 1) * config->scale;

			vertices[i + 0] = Vertex(x0, y0);
			vertices[i + 1] = Vertex(x0, y1);
			vertices[i + 2] = Vertex(x1, y1);
			vertices[i + 3] = Vertex(x1, y0);
		}

	glVertexPointer(2, GL_FLOAT, 0, vertices.data());
	glColorPointer(3, GL_FLOAT, 0, colors.data());
}

void Fluid::lin_solve(int b, float* x, float* x0, float a, float c)
{
	for (int k = 0; k < 2; ++k)
	{
		for (int i = 1; i < H - 1; ++i)
		{
			for (int j = 1; j < W - 1; ++j)
			{
				x[IX(j, i)] = (x0[IX(j, i)] + a *
					(x[IX(j - 1, i)] +
						x[IX(j + 1, i)] +
						x[IX(j, i - 1)] +
						x[IX(j, i + 1)])) / c;
			}
		}
		set_bnd(b, x);
	}
}

void Fluid::set_bnd(int b, float* x)
{
	int i;
	for (i = 1; i < H - 1; ++i)
	{
		x[IX(0,		i)] = b == 1 ? -x[IX(1, i	 )] : x[IX(1, i	   )];
		x[IX(W - 1, i)] = b == 1 ? -x[IX(W - 2, i)] : x[IX(W - 2, i)];
	}

	for (i = 1; i < W - 1; ++i)
	{
		x[IX(i, 0	 )] = b == 2 ? -x[IX(i, 1	 )] : x[IX(i, 1	   )];
		x[IX(i, H - 1)] = b == 2 ? -x[IX(i, H - 2)] : x[IX(i, H - 2)];
	}

	x[IX(0,		0	 )] = 0.33f * (x[IX(1,	   0	)] + x[IX(0,	 1	  )]);
	x[IX(0,		H - 1)] = 0.33f * (x[IX(1,	   H - 1)] + x[IX(0,	 H - 2)]);
	x[IX(W - 1, 0	 )] = 0.33f * (x[IX(W - 2, 0	)] + x[IX(W - 1, 1	  )]);
	x[IX(W - 1, H - 1)] = 0.33f * (x[IX(W - 2, H - 1)] + x[IX(W - 1, H - 2)]);

}

void Fluid::diffuse(int b, float* x, float* x0, float dt)
{
	float a = dt * diff * (W - 2) * (H - 2);
	lin_solve(b, x, x0, a, 1 + 6 * a);
}

void Fluid::advect(int b, float* d, float* d0, float* vx, float* vy, float dt)
{
	float i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1;

	float dtx = dt * float(W - 2);
	float dty = dt * float(H - 2);

	float wf = float(W);
	float hf = float(H);

	for (int i = 1; i < H - 1; ++i) 
	{
		for (int j = 1; j < W - 1; ++j) 
		{
			x = float(j) - dtx * vx[IX(j, i)]; 
			y = float(i) - dty * vy[IX(j, i)];

			if (x < 0.5f) x = 0.5f; 
			if (x > wf + 0.5f) x = wf + 0.5f;

			i0 = std::floorf(x); 
			i1 = i0 + 1.0f;

			if (y < 0.5f) y = 0.5f;
			if (y > hf + 0.5f) y = hf + 0.5f;
			
			j0 = std::floorf(y);
			j1 = j0 + 1.0f;

			s1 = x - i0; 
			s0 = 1.0f - s1; 
			t1 = y - j0; 
			t0 = 1.0f - t1;

			int i0i = int(i0);
			int i1i = int(i1);
			int j0i = int(j0);
			int j1i = int(j1);

			d[IX(j, i)] =
				s0 * (t0 * d0[safe_IX(i0i, j0i)] + t1 * d0[safe_IX(i0i, j1i)]) +
				s1 * (t0 * d0[safe_IX(i1i, j0i)] + t1 * d0[safe_IX(i1i, j1i)]);
		}
	}

	set_bnd(b, d);
}

void Fluid::project(float* vx, float* vy, float* p, float* div)
{
	for (int y = 1; y < H - 1; ++y)
	{
		for (int x = 1; x < W - 1; ++x) 
		{
			div[IX(x, y)] = -0.5f * (
				vx[IX(x + 1, y	  )] - 
				vx[IX(x - 1, y	  )] +
				vy[IX(x,	 y + 1)] - 
				vy[IX(x,	 y - 1)]) / ((W + H) * 0.5f);

			p[IX(x, y)] = 0;
		}
	}

	set_bnd(0, div); 
	set_bnd(0, p);

	lin_solve(0, p, div, 1, 6);

	for (int y = 1; y < H - 1; ++y)
	{
		for (int x = 1; x < W - 1; ++x) 
		{
			vx[IX(x, y)] -= 0.5f * (p[IX(x + 1, y	 )] - p[IX(x - 1, y	   )]) * W;
			vy[IX(x, y)] -= 0.5f * (p[IX(x,	    y + 1)] - p[IX(x,	  y - 1)]) * H;
		}
	}

	set_bnd(1, vx); 
	set_bnd(2, vy);
}

void Fluid::fade_density()
{
	std::for_each(std::execution::par_unseq,
		density.begin(), density.end(),
		[](float& d)
		{
			d = (d - 0.05f < 0) ? 0 : d - 0.05f;
		});
}

void Fluid::step_line(int x0, int y0, int x1, int y1, int dx, int dy, float a)
{
	if (x0 == x1 && y0 == y1)
		return;

	if (dx > dy)
	{
		add_velocity(x0, y0, (x1 - x0) * a, (y1 - y0) * a);

		int pk = 2 * dy - dx;
		for (int i = 0; i < dx; ++i)
		{
			x0 < x1 ? ++x0 : --x0;
			if (pk < 0)
				pk += 2 * dy;
			else
			{
				y0 < y1 ? ++y0 : --y0;
				pk += 2 * dy - 2 * dx;
			}
			add_velocity(x0, y0, (x1 - x0) * a, (y1 - y0) * a);
		}
	}
	else
	{
		add_velocity(x0, y0, (x1 - x0) * a, (y1 - y0) * a);

		int pk = 2 * dx - dy;
		for (int i = 0; i < dy; ++i)
		{
			y0 < y1 ? ++y0 : --y0;
			if (pk < 0)
				pk += 2 * dx;
			else
			{
				x0 < x1 ? ++x0 : --x0;
				pk += 2 * dx - 2 * dy;
			}
			add_velocity(x0, y0, (x1 - x0) * a, (y1 - y0) * a);
		}
	}
}

void Fluid::update(const float& dt)
{
	diffuse(1, vx_prev.data(), vx.data(), dt);
	diffuse(2, vy_prev.data(), vy.data(), dt);

	project(vx_prev.data(), vy_prev.data(), vx.data(), vy.data());
	
	advect(1, vx.data(), vx_prev.data(), vx_prev.data(), vy_prev.data(), dt);
	advect(2, vy.data(), vy_prev.data(), vx_prev.data(), vy_prev.data(), dt);

	project(vx.data(), vy.data(), vx_prev.data(), vy_prev.data());

	diffuse(0, density_prev.data(), density.data(), dt);
	advect(0, density.data(), density_prev.data(), vx.data(), vy.data(), dt);

	fade_density();
}

void Fluid::draw()
{
	std::for_each(std::execution::par_unseq,
		range.begin(), range.end(),
		[&](const int& i)
		{
			float r = map_to_range(vx[i], -0.05f, 0.05f, 0.0f, 1.0f);
			float g = map_to_range(vy[i], -0.05f, 0.05f, 0.0f, 1.0f);

			int vy = i * 4;

			colors[vy + 0] = Color(r, g, 0.0f);
			colors[vy + 1] = Color(r, g, 0.0f);
			colors[vy + 2] = Color(r, g, 0.0f);
			colors[vy + 3] = Color(r, g, 0.0f);
		});

	glDrawArrays(GL_QUADS, 0, V);
}
