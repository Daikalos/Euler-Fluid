#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>

#include "Camera.h"
#include "InputHandler.h"
#include "Fluid.h"
#include "Cell.h"

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

	window.setVerticalSyncEnabled(true);

	sf::CircleShape circle(100.0f);
	circle.setFillColor(sf::Color::Red);

	Camera camera(window);
	InputHandler input_handler;
	
	sf::Vector2f mouse_pos = sf::Vector2f(camera.get_mouse_world_position());;
	sf::Vector2f mouse_pos_prev = mouse_pos;

	sf::Clock clock;
	float dt = FLT_EPSILON;

	Fluid fluid(video_mode.size.x / 10.0f, video_mode.size.y / 10.0f, 0.0f, 0.0000001f);

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

		mouse_pos = sf::Vector2f(camera.get_mouse_world_position());
		sf::Vector2f amount = mouse_pos - mouse_pos_prev;
		mouse_pos_prev = mouse_pos;

		if (input_handler.get_left_held())
			fluid.add_density(mouse_pos.x / 10.0f, mouse_pos.y / 10.0f, 200.0f);

		fluid.add_velocity(mouse_pos.x / 10.0f, mouse_pos.y / 10.0f, amount.x / 10.0f, amount.y / 10.0f);
		fluid.update(dt);
		
		glClear(GL_COLOR_BUFFER_BIT);

		window.pushGLStates();



		window.popGLStates();

		glPushMatrix();

		glLoadMatrixf(camera.get_world_matrix());

		fluid.draw(window);

		glPopMatrix();

		window.display();
	}

	return 0;
}