#include "Fluid.h"

Fluid::Fluid(const size_t& width, const size_t& height, const float& diff, const float& visc)
	: W(width), H(height), N(width * height), diff(diff), visc(visc)
{
	u.resize(N, 0.0f);
	v.resize(N, 0.0f);

	u_prev.resize(N, 0.0f);
	v_prev.resize(N, 0.0f);

	density.resize(N, 0.0f);
	density_prev.resize(N, 0.0f);

	rectangles.reserve(N);
	for (int y = 0; y < H; ++y)
		for (int x = 0; x < W; ++x)
		{
			int index = IX(x, y);

			sf::RectangleShape rect;

			rect.setSize(sf::Vector2f(10.0f, 10.0f));
			rect.setPosition(sf::Vector2f(x * 10.0f, y * 10.0f));

			rectangles.push_back(rect);
		}
}

void Fluid::diffuse(int b, float* x, float* x0, float dt)
{
	float a = dt * diff * (W - 2) * (H - 2);
	lin_solve(b, x, x0, a, 1 + 4 * a);
}

void Fluid::lin_solve(int b, float* x, float* x0, float a, float c)
{
	for (int k = 0; k < 20; ++k)
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

void Fluid::advect(int b, float* d, float* d0, float* u, float* v, float dt)
{
	float i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1;

	float dtx = dt * (W - 2);
	float dty = dt * (H - 2);

	for (int i = 1; i < H - 1; ++i) 
	{
		for (int j = 1; j < W - 1; ++j) 
		{
			x = float(i) - dtx * u[IX(i, j)]; 
			y = float(j) - dty * v[IX(i, j)];

			if (x < 0.5f) x = 0.5f; 
			if (x > W + 0.5f) x = W + 0.5f; 

			i0 = std::floorf(x); 
			i1 = i0 + 1.0f;

			if (y < 0.5f) y = 0.5f;
			if (y > H + 0.5f) y = H + 0.5f; 
			
			j0 = std::floorf(y);
			j1 = j0 + 1.0f;

			s1 = x - i0; 
			s0 = 1.0f - s1; 
			t1 = y - j0; 
			t0 = 1.0f - t1;

			d[IX(i, j)] = 
				s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
				s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
		}
	}

	set_bnd(b, d);
}

void Fluid::project(float* u, float* v, float* p, float* div)
{
	for (int y = 1; y < H - 1; ++y)
	{
		for (int x = 1; x < W - 1; ++x) 
		{
			div[IX(x, y)] = -0.5f * (
				u[IX(x + 1, y	 )] - 
				u[IX(x - 1, y	 )] +
				v[IX(x,		y + 1)] - 
				v[IX(x,		y - 1)]) / W;

			p[IX(x, y)] = 0;
		}
	}

	set_bnd(0, div); 
	set_bnd(0, p);

	lin_solve(0, p, div, 1, 4);

	for (int y = 1; y < H - 1; ++y)
	{
		for (int x = 1; x < W - 1; ++x) 
		{
			u[IX(x, y)] -= 0.5f * (p[IX(x + 1, y	)] - p[IX(x - 1, y	  )]) * W;
			v[IX(x, y)] -= 0.5f * (p[IX(x,	   y + 1)] - p[IX(x,	 y - 1)]) * H;
		}
	}

	set_bnd(1, u); 
	set_bnd(2, v);
}

void Fluid::set_bnd(int b, float* x)
{
	int i;
	for (i = 1; i < W - 1; ++i) 
	{
		x[IX(i,		0)] = b == 2 ? -x[IX(i, 1	 )] : x[IX(i, 1	   )];
		x[IX(i, H - 1)] = b == 2 ? -x[IX(i, H - 2)] : x[IX(i, H - 2)];
	}

	for (i = 1; i < H - 1; ++i)
	{
		x[IX(0,		i)] = b == 1 ? -x[IX(1,		i)] : x[IX(1,	  i)];
		x[IX(W - 1, i)] = b == 1 ? -x[IX(W - 2, i)] : x[IX(W - 2, i)];
	}

	x[IX(0,		0	 )] = 0.5f * (x[IX(1, 0	   )] + x[IX(0,	    1)]);
	x[IX(0,		H - 1)] = 0.5f * (x[IX(1, H + 1)] + x[IX(0,	    H)]);
	x[IX(W - 1, 0	 )] = 0.5f * (x[IX(W, 0	   )] + x[IX(W + 1, 1)]);
	x[IX(W - 1, H - 1)] = 0.5f * (x[IX(W, H + 1)] + x[IX(W + 1, H)]);

}

void Fluid::update(const float& dt)
{
	diffuse(1, u_prev.data(), u.data(), dt);
	diffuse(2, v_prev.data(), v.data(), dt);

	project(u_prev.data(), v_prev.data(), u.data(), v.data());
	
	advect(1, u.data(), u_prev.data(), u_prev.data(), v_prev.data(), dt);
	advect(2, v.data(), v_prev.data(), u_prev.data(), v_prev.data(), dt);

	project(u.data(), v.data(), u_prev.data(), v_prev.data());

	diffuse(0, density_prev.data(), density.data(), dt);
	advect(0, density.data(), density_prev.data(), u.data(), v.data(), dt);
}

void Fluid::draw(sf::RenderWindow& window)
{
	std::for_each(rectangles.begin(), rectangles.end(),
		[&](sf::RectangleShape& rect)
		{
			int i = (&rect - rectangles.data());

			int r = (int)map_to_range(u[i], -0.05, 0.05f, 0, 255);
			int g = (int)map_to_range(v[i], -0.05, 0.05f, 0, 255);

			rect.setFillColor(sf::Color(r, g, 255));

			window.draw(rect);
		});
}
