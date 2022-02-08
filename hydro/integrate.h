#pragma once
#include "filewriter.h"

struct butcher
{
	std::vector<std::vector<double>> a;
	std::vector<double> b;
	std::vector<double> c;

	inline static const butcher& euler()
	{
		static const auto a = butcher{ { { 0 } },{ 1 },{ 0 } };
		return a;
	}
	inline static const butcher& heun2()
	{
		static const auto a = butcher{ { { 0,0 },{1,0} },{ 0.5,0.5 },{ 0,1 } };
		return a;
	}
	inline static const butcher& simpson()
	{
		static const auto a = butcher{ { { 0,0,0 },{ 0.5,0,0 },{ -1,2,0 }}, { 1. / 6, 4. / 6, 1. / 6 }, {0,0.5,1} };
		return a;
	}
	inline static const butcher& rk4()
	{
		static const auto a = butcher{ { { 0,0,0,0 },{ 0.5,0,0,0 },{ 0,0.5,0,0 },{0,0,1,0} }, { 1. / 6, 1. / 3, 1. / 3, 1. / 6 },{ 0,0.5,0.5,1 } };
		return a;
	}
};

void integrate(const butcher& butch, std::string filename)
{
	auto f = [](double t, double y) {return sin(y) - 0.7 * tanh(y * 100) * pow(t, 1. / 5) - cos(y); };
	std::vector<double> ts;
	std::vector<double> ys;
	ts.push_back(1);
	ys.push_back(1);
	int s = butch.b.size();

	std::vector<double> k(s);

	while (ts.back() < 50)
	{
		double dt = 0.01;

		double t = ts.back();
		double y = ys.back();
		double temp = 0;
		for (int j = 0; j < s; j++)
		{
			double temp2 = 0;
			for (int l = 0; l < j; l++)
			{
				temp2 += butch.a[j][l] * k[l];
			}
			k[j] = f(t + dt * butch.c[j], y + dt * temp2);
			temp += butch.b[j] * k[j];
		}
		y += dt * temp;
		ys.push_back(y);
		t += dt;
		ts.push_back(t);
	}
	CubeFileWriter fw;
	fw.setShape({ (int)ys.size() });
	fw.open(filename);
	fw.new_row();
	fw.to_data(ts);
	fw.new_row();
	fw.to_data(ys);
	fw.close();

	std::cout << " integration done" << std::endl;
}