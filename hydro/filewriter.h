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

	std::vector<size_t> dim;
	std::vector<std::vector<double>> header;
	std::fstream myfile;

	void writeDims() {
		size_t d = dim.size();
		myfile.seekp(0); 
		myfile.write((char*)&d, sizeof(size_t));
		myfile.write((char*)&dim[0], sizeof(size_t) * d);
	}

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
		T& operator* ()
		{
			return (*m_container)[m_pos];
		}
	};

	void setShape(std::vector<size_t> subelementDimensions)
	{
		dim = subelementDimensions;
		dim.insert(dim.begin(), 0);
		header.resize(dim.size());
		for (size_t i = 0; i < dim.size(); i++)
			header[i].resize(dim[i], 0);
	}
	void open(std::string filename)
	{
		myfile = std::fstream(filename, std::ios::out | std::ios::binary);
		writeDims();
	}

	void newLayer()
	{
		// starts a new layer of dimensions "dim" in the hypercube
		dim[0]++;
		header[0].push_back(0);
	}

	template<typename T>
	void writeData(const std::vector<T>& data)
	{
		// appends data to the current layer
		if(!data.empty())
		myfile.write((char*)&data[0], sizeof(T)*data.size());
	}

	template<typename T>
	void writeData(const block_iter<T>& begin, const block_iter<T>& end)
	{
		// appends data to the current layer
		for (auto it = begin; it != end; ++it)
			myfile.write((char*)&(*it), sizeof(T));
	}

	void writeHeader(double value, int d, int i = -1)
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
			for (size_t i = 0; i < dim.size(); i++)
				myfile.write((char*)&header[i][0], sizeof(double)*header[i].size());

			// Update the dimension in the file
			writeDims();
			myfile.close();
		}
		else {
			std::cout << "WARNING: The file has never been opened." << std::endl;
		}
	}

	void reset()
	{
		dim.clear();
		header.clear();
		myfile.clear();
	}
};

// Writes a stream of blocks of specified dimensions. The stream can be organized in sections
class StreamFileWriter
{
	std::vector<size_t> dim;
	std::vector<double> annotations;
	std::fstream myfile;

	std::vector<int> sectionSizes;
	int currentSectionSize;

	void writeHeader() {
		auto oldpos = myfile.tellp();

		size_t d = dim.size();
		myfile.seekp(0);
		myfile.write((char*)&d, sizeof(size_t));
		myfile.write((char*)&dim[0], sizeof(size_t) * d);

		size_t s = sectionSizes.size();
		myfile.write((char*)&s, sizeof(size_t));

		if(oldpos!=0)
			myfile.seekp(oldpos);
	}
	void writeFooter() {
		myfile.write((char*)&annotations[0], sizeof(double) * annotations.size());
		myfile.write((char*)&sectionSizes[0], sizeof(int) * sectionSizes.size());
	}

public:
	void setBlockShape(std::vector<size_t> blockDimensions){
		dim = blockDimensions;
		dim.insert(dim.begin(), 0);
	}
	void open(std::string filename){
		myfile = std::fstream(filename, std::ios::out | std::ios::binary);
		writeHeader();
		currentSectionSize = 0;
	}

	template<typename T>
	void writeBlocks(const std::vector<T>& blocks){
		// appends data to the current section
		if (!blocks.empty() && myfile.is_open()) {
			myfile.write((char*)&blocks[0], sizeof(T) * blocks.size());
			dim[0] += blocks.size();
			currentSectionSize += blocks.size();
		}
	}

	void finishSection(double annotation){
		// Binds all blocks written since the last call or the beginning into a section.
		sectionSizes.push_back(currentSectionSize);
		annotations.push_back(annotation);
		currentSectionSize = 0;
	}

	void close(){
		if (myfile.is_open()){
			// Update the dimension in the file
			writeHeader();
			writeFooter();
			myfile.close();
		}
		else {
			std::cout << "WARNING: The file has never been opened." << std::endl;
		}
	}

	void reset(){
		dim.clear();
		annotations.clear();
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
