#include "Config.h"

Config::Config()
{
	scale = 10;
}

Config::~Config()
{

}

void Config::load()
{
	std::ifstream project_file(FILE_NAME + ".json", std::ifstream::binary);
	if (project_file.good())
	{
		try
		{
			nlohmann::json json = nlohmann::json::parse(project_file);
			load_var(json);
		}
		catch (nlohmann::json::parse_error) {}
		catch (nlohmann::detail::type_error e) {}
	}
}

void Config::load_var(nlohmann::json& json)
{
	nlohmann::basic_json<>::value_type config = json[FILE_NAME];

	scale = config["scale"];
}

std::vector<Reconstruct> Config::refresh(Config& prev)
{
	std::vector<Reconstruct> result;
	result.reserve(8);

	std::ifstream project_file(FILE_NAME + ".json", std::ifstream::binary);
	if (project_file.good())
	{
		try
		{
			nlohmann::json json = nlohmann::json::parse(project_file);
			load_var(json);
		}
		catch (nlohmann::json::parse_error e) {}
		catch (nlohmann::detail::type_error e) {}
	}

	if (prev.scale != scale)
		result.push_back(Reconstruct::RFluid);

	return result;
}
