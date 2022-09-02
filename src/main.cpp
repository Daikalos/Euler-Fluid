#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>

#include "Camera.h"
#include "InputHandler.h"
#include "Fluid.h"
#include "Config.h"

int main()
{
	srand(time(NULL));

	std::vector<sf::VideoMode> fullscreen_modes = sf::VideoMode::getFullscreenModes();

	if (fullscreen_modes.size() == 0)
		return -1;

	sf::VideoMode video_mode = fullscreen_modes.front();

	if (!video_mode.isValid())
		return -1;

	sf::RenderWindow window(video_mode, "Euler Fluid", sf::Style::Fullscreen);

	if (!window.setActive(true))
		return -1;

	Config config;
	config.load();

	window.setVerticalSyncEnabled(true);

	Camera camera(window);
	InputHandler input_handler;
	
	sf::Vector2i mouse_pos = camera.get_mouse_world_position() / config.scale;
	sf::Vector2i mouse_pos_prev = mouse_pos;

	sf::Clock clock;
	float dt = FLT_EPSILON;

	Fluid fluid(&config, 2 + video_mode.size.x / config.scale, 2 + video_mode.size.y / config.scale, 0.0f, 0.0000001f);

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glScalef(1.0f, -1.0f, 1.0f);
	glOrtho(0, window.getSize().x, 0, window.getSize().y, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	int i = 0;
	while (window.isOpen())
	{
		dt = std::fminf(clock.restart().asSeconds(), 0.075f);

		input_handler.update();

		sf::Event event;
		while (window.pollEvent(event))
		{
			switch (event.type)
			{
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::Resized:
				glViewport(0, 0, window.getSize().x, window.getSize().y);
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				glScalef(1.0f, -1.0f, 1.0f);
				glOrtho(0, window.getSize().x, 0, window.getSize().y, -1.0, 1.0);
				glMatrixMode(GL_MODELVIEW);

				camera.set_position((sf::Vector2f)window.getSize() / 2.0f);

				break;
			case sf::Event::MouseWheelScrolled:
				input_handler.set_scrollDelta(event.mouseWheelScroll.delta);
				break;
			}
		}

		camera.update(input_handler);

		mouse_pos = camera.get_mouse_world_position() / config.scale;
		sf::Vector2i amount = mouse_pos - mouse_pos_prev;

		if (input_handler.get_left_held())
			fluid.add_density(mouse_pos.x, mouse_pos.y, 200.0f);

		fluid.step_line(
			mouse_pos_prev.x, mouse_pos_prev.y,
			mouse_pos.x, mouse_pos.y, 
			std::fabsf(amount.x), std::fabsf(amount.y), 1.0f);

		mouse_pos_prev = mouse_pos;

		fluid.update(dt);
		
		glClear(GL_COLOR_BUFFER_BIT);
		glPushMatrix();
		glLoadMatrixf(camera.get_world_matrix());
		fluid.draw();
		glPopMatrix();
		window.display();
	}

	return 0;
}