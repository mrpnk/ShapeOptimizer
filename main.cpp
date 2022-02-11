#include "ShapeOptimizer.h"

#include "hydro/sph.h"

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <random>
#include <set>

std::vector<vec<2>> myshape;
double myscale = 0.05;
double solidMass = 0.2;

double simulateDrag(std::vector<vec<2>> const& shape)
{
	srand(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	srand(3245);

	using namespace Fluid;
	SPH<2> sph;
	sph.init(1.4, 3.0, 16, false);

	sph.setDomain({ {{0,1},{0,0.6}} },
		{ { {BoundaryType::inflow, BoundaryType::outflow},
			{BoundaryType::wall, BoundaryType::wall} } });

	sph.createParticles(100, 1.0);
	sph.createSolid(shape, solidMass, myscale);

	StreamFileWriter pfw;
	StreamFileWriter sfw;

	int logCounter = 0;
	double avgDrag = 0;
	sph.simulate<false, true>(pfw, sfw, 4, 2.0, 0.8, 0.01, [&](logInfo const& li) {
		if (li.time > 1) {
			logCounter++; avgDrag += li.drag.x;
		}});
	avgDrag /= logCounter;

	std::cout << "drag = " << avgDrag << std::endl;

	return avgDrag;
}

struct oi {
	double energy;
	double nFilled;
};
struct optimizationInfo {
	
	std::vector<oi> infos;
	friend std::ostream& operator << (std::ostream& os, optimizationInfo const& odi) {
		for (oi const& i : odi.infos)
			os << i.energy << " " << i.nFilled << std::endl;
		return os;
	}
};

struct edge {
	vec<2> a, b;
	friend std::ostream& operator << (std::ostream& os, edge const& ed) {
		return os << ed.a << " " << ed.b;
	}
};
struct boundary {
	std::vector<edge> edges;
	double getDistance() {
		double dist=0;
		for (auto const& e:edges)
			dist += sqrt(pow(e.a.x - e.b.x, 2) + pow(e.a.y - e.b.y, 2));
		return dist;
	}
	friend std::ostream& operator << (std::ostream& os, boundary const& bd) {
		for (auto const& e : bd.edges)
			os << e << std::endl;
		return os << std::endl;
	}
};
class scene {
	int w, h;// width and height
	std::vector<std::vector<bool>> dataset;
	std::set<std::pair<int, int>> filleds, emptys;

	mutable std::random_device rd;

	bool isBoundary(int x, int y) { // between filled and empty
		return (dataset[y][x] != dataset[y][x - 1])
			|| (dataset[y][x] != dataset[y][x + 1])
			|| (dataset[y][x] != dataset[y - 1][x])
			|| (dataset[y][x] != dataset[y + 1][x]);
	}

	std::mt19937 gen;
	std::uniform_int_distribution<> uid, uid2;
	std::uniform_real_distribution<> urd;

	bool isConnected() {
		std::vector<std::vector<bool>> visited;
		visited.resize(h);
		for (int y = 0; y < h; ++y) {
			visited[y].resize(w,false);
		}
		
		std::queue<std::pair<int, int>> toSpread;

		// find some filled cell
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				if (dataset[y][x]) {
					visited[y][x] = true;
					toSpread.push({x,y});
					goto breakout;
				}
			}
		}
	breakout:

		// propagate the "visited" status to neighbouring cells
		while (!toSpread.empty()) {
			auto[x,y] = toSpread.front();
			toSpread.pop();

			if (dataset[y - 1][x] && !visited[y - 1][x]) {
				visited[y - 1][x] = true;
				toSpread.push({ x,y-1 });
			}
			if (dataset[y + 1][x] && !visited[y + 1][x]) {
				visited[y + 1][x] = true;
				toSpread.push({ x,y + 1 });
			}
			if (dataset[y][x - 1] && !visited[y][x - 1]) {
				visited[y][x - 1] = true;
				toSpread.push({ x - 1,y });
			}
			if (dataset[y][x + 1] && !visited[y][x + 1]) {
				visited[y][x + 1] = true;
				toSpread.push({ x + 1,y });
			}
		}


		// check if there are filled cells that are not visited
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				if (dataset[y][x] && !visited[y][x]) {
					return false;
				}	
			}
		}
		return true;
	}

	void fillCell(int x, int y) {
		dataset[y][x] = true;
		filleds.insert({ x,y });
		emptys.erase({ x,y });
	}
	void emptyCell(int x, int y) {
		dataset[y][x] = false;
		emptys.insert({ x,y });
		filleds.erase({ x,y });
	}

	void log(StreamFileWriter& sfw, StreamFileWriter& sfw2, int i) {
		// Write to files
		auto sh = getShape(false);
		sfw2.writeBlocks(sh);
		sfw2.finishSection(i);

		boundary b = generateBoundary();
		sfw.writeBlocks(b.edges);
		sfw.finishSection(i);
	}

public:
	void init(int w, int h) {
		gen.seed(324);// rd()); // seed

		uid = std::uniform_int_distribution<>(0, (w - 2) * (h - 2) - 1);
		urd = std::uniform_real_distribution<>(0, 1);
		uid2 = std::uniform_int_distribution<>(0, 999999);

		this->w = w; this->h = h;

		// initialize randomly
		std::mt19937 gen(rd()); // seed
		std::uniform_int_distribution<> distrib(0, 1);

		dataset.resize(h);
		for (int y = 0; y < h; ++y) {
			dataset[y].resize(w);
			for (int x = 0; x < w; ++x) {
				
				// place a square in the center
				if (abs(x - w / 2) < 3 && abs(y - h / 2) < 3)
					fillCell(x,y);
				else if (x * y == 0 || x == w - 1 || y == h - 1)
					dataset[y][x] = false; // dont write into empties
				else emptyCell(x, y);
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
				vec<2> centerNodes[4] = { {x + 0.5,(double)y},{x + 1.,y + 0.5},{(double)x,y + 0.5},{(double)x + 0.5,y + 1.} };
				size_t idx = (dataset[y][x] << 3) + (dataset[y][x + 1] << 2) + (dataset[y + 1][x] << 1) + (dataset[y + 1][x + 1]);
				for (auto const& a : lu[idx]) {
					edge e;
					e.a = centerNodes[a[0]];
					bd.edges.push_back(edge{ centerNodes[a[0]],centerNodes[a[1]] });
				}
			}
		}

		return bd;
	}

	double getVolume() {
		double vol=0;
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				vol += dataset[y][x];
			}
		}
		return vol;
	}
	double getEnergy() { // The quantity to minimize
		double targetVolume = 25;
		double targetVolumeTolerance = 2;


		double deviation = std::max(0., abs(getVolume() - targetVolume) - targetVolumeTolerance);

		return simulateDrag(getShape()) * (1 + deviation * 0.1);
		//return generateBoundary().getDistance()/getVolume() *(1+abs(getVolume()-100)*0.01);
	}

	void findCell(int& x,int& y, bool filled, bool onlyBoundary) {
		
		if (filled) {
			auto it = std::begin(filleds);
			std::advance(it, uid2(gen) % filleds.size());
			std::tie(x, y) = *it;
		}
		else {
			auto it = std::begin(emptys);
			std::advance(it, uid2(gen) % emptys.size());
			std::tie(x, y) = *it;
		}

		if (onlyBoundary)
		while (!isBoundary(x, y)) {
			findCell(x, y, filled, false);
		}
	}

	void optimize(StreamFileWriter& sfw, StreamFileWriter& sfw2, int nIterations, double temperature, bool onlyBoundary, optimizationInfo& out_info) {
		// Use Monte-Carlo
		AutoTimer at(g_timer, "optimize"); 
		
		int x0, y0, x1, y1;
		int accCounter = 0;

		log(sfw, sfw2,0);

		double oldEnergy = getEnergy(), newEnergy;
		for (int i = 0; i < nIterations; ++i) {
		
			findCell(x0, y0, true, onlyBoundary);
			findCell(x1, y1, false, onlyBoundary);

			// move the filling
			emptyCell(x0, y0);
			fillCell(x1, y1);

			//std::cout << "Propose: " << x0 << " " << y0 << " -> " << x1 << " " << y1 << std::endl;

			bool accepted = false;
			if (isConnected()) {
				newEnergy = getEnergy(); // compute new energy
				double deltaEnergy = newEnergy - oldEnergy;
				double prop = std::min(1., exp(-deltaEnergy / temperature));
				accepted = urd(gen) < prop;
			}
			else {
			/*	std::cout << "Propose: " << x0 << " " << y0 << " -> " << x1 << " " << y1 << std::endl;

				std::cout << "would disconnect" << std::endl; */

				emptyCell(x1, y1);
				fillCell(x0, y0);

				i--;
				continue;
			}
			if (accepted) {
				oldEnergy = newEnergy;
				std::cout << "Accepted!" <<std::endl;
				accCounter++;
			}
			else { // rejected: move it back
				emptyCell(x1, y1);
				fillCell(x0, y0);
			}
			out_info.infos.push_back({ oldEnergy,getVolume() });


			// write the state to files
			log(sfw, sfw2, i+1);
			

			//if (accCounter == 15) break;
		}
	}

	std::vector<vec<2>> getShape(bool center = true) {
		std::vector<vec<2>> shape;
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				if(dataset[y][x])
					shape.push_back({ (double)x - (double)w / 2* center,(double)y - (double)h / 2* center });
			}
		}
		return shape;
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

// ---------------------------------------------------------------

void windchannel()
{
	using namespace Fluid;
	SPH<2> sph;
	sph.init(1.4, 3.0, 16, true); // dry air, 0 degree Celsius, normal pressure
	
	// No velocity field. Only 'pumping' the particles from the right boundary to the left,
	// creating a pressure gradient.

	sph.setDomain({ {{0,1},{0,.6}} },
		{ { {BoundaryType::inflow, BoundaryType::outflow},
			{BoundaryType::wall, BoundaryType::wall} } });

	sph.setExternalTime(0);
	sph.setExternalAcceleration({ 0,0 });
	sph.setExternalForce({ 0.01,0 });


	std::vector<vec<2>> dish = { {-7, -12}, { -2,-9 }, { 0,-5 }, { 0,0 }, { 0,5 }, { -2,9 }, { -7,12 } };
	std::vector<vec<2>> wing = { {-20, 3},{-15, 5},{-10, 5},{-5, 3},{0, 0},{5, -5} };
	std::vector<vec<2>> torpedo = { {-7, -12}, { -2,-9 }, { 0,-5 }, { 0,0 }, { 0,5 }, { -2,9 }, { -7,12 } };

	sph.createParticles(100, 1.0);
	sph.createSolid(myshape, solidMass, myscale);


	StreamFileWriter pfw;
	pfw.setBlockShape({ fileParticleData<2>::numElements()});
	pfw.open("particles.binary");

	StreamFileWriter sfw;
	sfw.setBlockShape({ fileSolidData<2>::numElements() });
	sfw.open("solids.binary");

	sph.simulate<false, true>(pfw, sfw, 4, 2.0, 0.8, 0.005);



	pfw.close();
	std::cout << "Particles written to file!" << std::endl;

	sfw.close();
	std::cout << "Solids written to file!" << std::endl;

	sph.getTree().toFile("kdtree.binary");
	std::cout << "Tree written to file!" << std::endl;
}

int main(){
	std::cout << "Max number of OMP threads: " << omp_get_max_threads() << std::endl;

	// --------------------------------------------
	size_t w = 20, h = 20;

	scene sc;
	sc.init(w,h);
	std::ofstream file;

	StreamFileWriter sfw;
	sfw.setBlockShape({ 2,2 });
	sfw.open("output-boundary.binary");

	StreamFileWriter sfw2;
	sfw2.setBlockShape({ 2 });
	sfw2.open("output-filling.binary");

	optimizationInfo oi;
	sc.optimize(sfw, sfw2, 1000, 0.004, true, oi);

	file.open("output-info.txt");
	file << oi;
	file.close();


	sfw.close();
	sfw2.close();

	myshape = sc.getShape();

	// --------------------------------------------
	
	//std::cout << "drag = " << simulateDrag(myshape);

	// --------------------------------------------

	windchannel();
	
	
	// --------------------------------------------
	
	g_timer.print();
	//std::cin.get();

	return 0;
}

