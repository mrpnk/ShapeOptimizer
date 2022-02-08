#pragma once
#include <fstream>
#include <string>

class CubeFileWriter
{
	struct i2
	{
		int x;
		int y;
	};

	std::vector<int> sizes;
	std::vector<std::vector<double>> header;
	std::fstream myfile;
public:

	template<typename T>
	class block_iter
	{
		std::vector<T>* m_container;

		i2 m_size;
		i2 m_step;
		i2 m_padding;

		int m_pos;

	public:
		block_iter(std::vector<T> &container, i2 size, int pos, i2 step = { 1,1 }, i2 padding = { 0,0 }) :
			m_container(&container),
			m_pos(pos),
			m_size(size),
			m_step(step),
			m_padding(padding) {}

		block_iter(const block_iter & other) :
			m_container(other.m_container),
			m_pos(other.m_pos),
			m_size(other.m_size),
			m_step(other.m_step),
			m_padding(other.m_padding) {}

		block_iter& operator= (const block_iter& other)
		{
			this->m_container = other.m_container;
			this->m_pos = other.m_pos;
			return *this;
		}

		block_iter& operator++ ()
		{
			int w = (m_size.x + m_padding.x * 2);
			if ((m_pos + m_step.x) % w >= m_size.x + m_padding.x)
				if (m_pos / w + m_step.y >= m_size.y + m_padding.y)
					m_pos = w * (m_size.y + 2 * m_padding.y);
				else
					m_pos += w * m_step.y - m_pos % w + m_padding.x;
			else
				m_pos += m_step.x;

			return *this;
		}
		bool operator!= (const block_iter& other)
		{
			return  m_pos != other.m_pos || m_container != other.m_container;
		}
		T&    operator* ()
		{
			return (*m_container)[m_pos];
		}
	};

	void setShape(std::vector<int> subelementDimensions)
	{
		sizes = subelementDimensions;
		sizes.insert(sizes.begin(), 0);
		header.resize(sizes.size());
		for (size_t i = 0; i < sizes.size(); i++)
			header[i].resize(sizes[i], 0);
	}
	void open(std::string filename)
	{
		int d = sizes.size();
		myfile = std::fstream(filename, std::ios::out | std::ios::binary);
		myfile.write((char*)&d, sizeof(int));
		myfile.write((char*)&sizes[0], sizeof(int)*d);
	}

	void new_row()
	{
		// starts a new row in the hypercube
		sizes[0]++;
		header[0].push_back(0);
	}

	template<typename T>
	void to_data(const std::vector<T>& data)
	{
		// appends data to the current row
		myfile.write((char*)&data[0], sizeof(T)*data.size());
	}

	template<typename T>
	void to_data(const block_iter<T>& begin, const block_iter<T>& end)
	{
		// appends data to the current row
		for (auto it = begin; it != end; ++it)
			myfile.write((char*)&(*it), sizeof(T));
	}

	void to_header(double value, int d, int i = -1)
	{
		// labels a header in dimension d with index i
		if (i < 0)
			i = header[d].size() - 1;
		header[d][i] = value;
	}

	void close()
	{
		if (myfile.is_open())
		{
			for (size_t i = 0; i < sizes.size(); i++)
				myfile.write((char*)&header[i][0], sizeof(double)*header[i].size());

			// Update the dimension in the file
			int d = sizes.size();
			myfile.seekp(sizeof(int) * 1);
			myfile.write((char*)&sizes[0], sizeof(int)*d);
			myfile.close();
		}
	}

	void reset()
	{
		sizes.clear();
		header.clear();
		myfile.clear();
	}
};

class NestFileWriter
{
	std::fstream myfile;
	std::vector<int> layout;

public:

	void open(std::string filename)
	{
		int pos = 0;
		myfile = std::fstream(filename, std::ios::out | std::ios::binary);
		myfile.write((char*)&pos, sizeof(int));
	}

	void beginList()
	{
		layout.push_back(0);
	}
	void endList()
	{
		layout.push_back(-1);
	}
	template<typename T>
	void to_data(const std::vector<T>& data)
	{
		myfile.write((char*)&data[0], sizeof(T)*data.size());
		layout.push_back(data.size()*(sizeof(T)/sizeof(double)));
	}

	void close()
	{
		if (myfile.is_open())
		{
			int pos = 0;
			for (int& i : layout)
				if(i>0)pos += i;

			myfile.write((char*)&layout[0], sizeof(int)*layout.size());

			myfile.seekp(sizeof(double) * 0);
			myfile.write((char*)&pos, sizeof(int));
			myfile.close();
		}
	}

	void reset()
	{
		layout.clear();
		myfile.clear();
	}
};
