#include "ShapeOptimizer.h"

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <random>

struct vec2 {
	float x, y;
	friend std::ostream& operator << (std::ostream& os, vec2 const& ve) {
		return os << ve.x << " " << ve.y;
	}
};
struct edge {
	vec2 a, b;
	friend std::ostream& operator << (std::ostream& os, edge const& ed) {
		return os << ed.a << " " << ed.b;
	}
};
struct boundary {
	std::vector<edge> edges;
	friend std::ostream& operator << (std::ostream& os, boundary const& bd) {
		for (auto const& e : bd.edges)
			os << e << std::endl;
		return os << std::endl;
	}
};
class scene {
	int w, h;// width and height
	std::vector<std::vector<bool>> dataset;
public:
	void init(int w, int h) {
		this->w = w; this->h = h;

		// initialize randomly
		std::random_device rd;
		std::mt19937 gen(214); // seed
		std::uniform_int_distribution<> distrib(0, 1);

		dataset.resize(h);
		for (int y = 0; y < h; ++y) {
			dataset[y].resize(w);
			for (int x = 0; x < w; ++x) {
				dataset[y][x] = (bool)distrib(gen);
				if (x * y == 0 || x == w - 1 || y == h - 1)
					dataset[y][x] = false;
			}
		}
	}
	boundary generateBoundary() {
		// Use marching squares
		boundary bd;
		const std::array<std::vector<std::array<int, 2>>, 16> lu = {{
		{ }, // nothing set
		{ {1, 3}}, // 3 set
		{ {3, 2}}, // 2 set
		{ {1, 2}}, // 3 and 2
		{ {0, 1}}, // 1
		{ {0, 3}}, // 1 and 3
		{ {0, 1}, { 3, 2 }}, // 1 and 2 set
		{ {0, 2}}, // 1,2,3
		{ {2, 0}}, // 0
		{ {2, 0}, { 1, 3 }}, // 0 and 3
		{ {0, 3}},
		{ {1, 0}},
		{ {2, 1}},
		{ {2, 3}},
		{ {3, 1}},
		{ } // all set
		}};

		for (int y = 0; y < h-1; ++y) {
			for (int x = 0; x < w - 1; ++x) {
				vec2 centerNodes[4] = { {x + 0.5f,y},{x + 1,y + 0.5f},{x,y + 0.5f},{x + 0.5f,y + 1} };
				size_t idx = (dataset[y][x] << 3) + (dataset[y][x + 1] << 2) + (dataset[y + 1][x] << 1) + (dataset[y + 1][x + 1]);
				auto const& edgeidx = lu[idx];
				for (auto const& a : edgeidx) {
					edge e = { centerNodes[a[0]],centerNodes[a[1]] };
					bd.edges.push_back(e);
				}
			}
		}

	
		return bd;
	}
	friend std::ostream& operator << (std::ostream& os, scene const& sc) {
		for (int y = 0; y < sc.h; ++y) {
			for (int x = 0; x < sc.w; ++x) {
				os << sc.dataset[y][x] << " ";
			}
			os << std::endl;
		}
		return os;
	}
};

int main(){
	scene sc;
	sc.init(20, 20);
	std::ofstream file;

	file.open("output-raw.txt");
	file << sc;
	file.close();


	boundary b = sc.generateBoundary();

	file.open("output-boundary.txt");
	file << b;
	file.close();

	
	return 0;
}
