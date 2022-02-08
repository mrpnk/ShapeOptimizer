#include <vector>
#include <string>
#include <variant>


const double pi = 3.14159265359;
int divround(int n, int d) { return (n + d - 1) / d; }
inline double minmod(double a, double b)
{
	return a * b <= 0 ? 0 : abs(a) < abs(b) ? a : b;
}

std::vector<std::variant<std::monostate,double, std::string>> getArgs(int argc, char* argv[], std::vector<std::string> flags)
{
	std::vector<std::variant<std::monostate,double, std::string>> args(flags.size());
	for (int i = 0; i < argc - 1; i++)
		for (size_t j = 0; j < flags.size(); j++)
		{
			std::string& f = flags[j];
			if (std::string(argv[i]) == f)
			{
				char* p;
				double converted = strtof(argv[i + 1], &p);
				if (*p) args[j] = std::string(argv[i + 1]);
				else args[j] = converted;
				i++;
			}
		}
	return args;
}

