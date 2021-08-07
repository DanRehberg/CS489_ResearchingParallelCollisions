/*
Author: Dan Rehberg
Date Modified: 8/5/2021

Purpose: Simple mesh class for shapes.
	Note, this is not a class that would be ideal to reuse with modifications.
	This uses manual memory allocation and limits the size of the
		mesh to the most complex mesh being tested in this program!
*/

#ifndef __SHAPE__
#define __SHAPE__

#include <cstdint>
#include "VectorMath.hpp"

struct Indices3
{
	int ind[3];
};

struct TriEdges
{
	vec3 edge[2];
};

class Shape
{
public:
	Shape()
	{
		vertices[0] = vec3(12.0f, 0.0f, 0.0f);
		if (edges == nullptr)
		{
			edges = new int* [242];
			for (unsigned int i = 0; i < 242; ++i)
			{
				edges[i] = new int[16];
			}
		}
		if (edgeAlt == nullptr)edgeAlt = new int[242 * 16];
		if (edgeCount == nullptr)edgeCount = new int[242];
	}
	Shape(vec3* positions, int N)
	{
		//vertices = new vec3[N];
		count = N;
		for (int i = 0; i < N; ++i)
		{
			vertices[i] = positions[i];
		}
		if (edges == nullptr)
		{
			edges = new int* [242];
			for (unsigned int i = 0; i < 242; ++i)
			{
				edges[i] = new int[16];
			}
		}
		if (edgeAlt == nullptr)edgeAlt = new int[242 * 16];
		if (edgeCount == nullptr)edgeCount = new int[242];
		faceCount = 2 * (N - 2);
		if (faces == nullptr)faces = new vec3[faceCount];
		if (faceVerts == nullptr)faceVerts = new Indices3[faceCount];
		if (faceEdges == nullptr)faceEdges = new TriEdges[faceCount];
	}
	~Shape()
	{
		//if (vertices != NULL)delete[] vertices;
		if (edges != nullptr)
		{
			for (unsigned int i = 0; i < 242; ++i)
			{
				delete[] edges[i];
			}
			delete[] edges;
		}
		if (edgeAlt != nullptr)
		{
			delete[] edgeAlt;
		}
		if (edgeCount != nullptr)
		{
			delete[] edgeCount;
		}
		if (faces != nullptr)
		{
			delete[] faces;
		}
		if (faceVerts != nullptr)
		{
			delete[] faceVerts;
		}
		if (faceEdges != nullptr)
		{
			delete[] faceEdges;
		}
	}
	//Host only on the copy constructor and copy assignment
	Shape(const Shape& cp)
	{
		if (cp.vertices != NULL)
		{
			this->count = cp.count;
			//Below should not happen, just adding it
			//if (this->vertices != NULL)delete[] this->vertices;
			//this->vertices = new vec3[cp.count];
			for (int i = 0; i < cp.count; ++i)
			{
				this->vertices[i] = cp.vertices[i];
			}
		}
		if (edges == nullptr)
		{
			edges = new int* [242];
			for (unsigned int i = 0; i < 242; ++i)
			{
				edges[i] = new int[16];
			}
		}
		if (edgeAlt == nullptr)edgeAlt = new int[242 * 16];
		if (edgeCount == nullptr)edgeCount = new int[242];
		faceCount = cp.faceCount;
		if (faceCount > 0)
		{
			faces = new vec3[faceCount];
			faceVerts = new Indices3[faceCount];
			faceEdges = new TriEdges[faceCount];
		}
	}
	Shape& operator=(const Shape& cp)
	{
		/*if (this->vertices == NULL)
		{
			this->count = cp.count;
			this->vertices = new vec3[cp.count];
		}
		else if (this->count != cp.count)
		{
			this->count = cp.count;
			delete[] this->vertices;
			this->vertices = new vec3[cp.count];
		}*/
		this->count = cp.count;
		for (int i = 0; i < cp.count; ++i)
		{
			this->vertices[i] = cp.vertices[i];
		}
		if (edges == nullptr)
		{
			edges = new int* [242];
			for (unsigned int i = 0; i < 242; ++i)
			{
				edges[i] = new int[16];
			}
		}
		if (edgeAlt == nullptr)edgeAlt = new int[242 * 16];
		if (edgeCount == nullptr)edgeCount = new int[242];
		faceCount = cp.faceCount;
		if (faceCount > 0)
		{
			faces = new vec3[faceCount];
			faceVerts = new Indices3[faceCount];
			faceEdges = new TriEdges[faceCount];
		}
		return *this;
	}
	vec3 getVertex(const uint32_t& index) const
	{
		return vertices[index];
	}
	/*
	__device__ vec3* getVertices()
	{
		return vertices;
	}
	*/
	uint32_t supportPoint(const vec3& direction) const
	{
		float magnitude = -999999;
		uint32_t tempID = 0;
		for (int i = 0; i < count; ++i)
		{
			float dR = dot(vertices[i], direction);
			if (dR > magnitude)
			{
				magnitude = dR;// dot(vertices[i], direction);
				tempID = i;
			}
		}
		return tempID;
	}
	uint32_t supportPoint(const vec3& direction, const vec3& t) const
	{
		float magnitude = -999999;
		uint32_t tempID = 0;
		for (int i = 0; i < count; ++i)
		{
			float dR = dot((vertices[i] + t), direction);
			if (dR > magnitude)
			{
				magnitude = dR;
				tempID = i;
			}
		}
		return tempID;
	}
	uint32_t supportPointHillClimb(const vec3& direction, const int& prevID)
	{
		float magnitude = -999999;
		float nMag = magnitude;
		uint32_t tempID = 0;
		uint32_t curID = prevID;
		uint32_t lastID = prevID;
		bool end = false;
		while (!end)
		{
			for (int i = 0; i < edgeCount[curID]; ++i)
			{
				float dR = dot(vertices[edges[curID][i]], direction);
				if (dR > nMag)
				{
					nMag = dR;
					tempID = edges[curID][i];
				}
			}
			if (nMag > magnitude)
			{
				lastID = curID;
				curID = tempID;
				magnitude = nMag;
			}
			else end = true;
		}
		return tempID;
	}
	vec3 vertices[242];
	int count = 0;
	int edgeTotal = 496;
	int faceCount = 0;
	//int edges[242][4];
	//int edgeCount[242];
	int** edges = nullptr;
	int* edgeAlt = nullptr;
	int* edgeCount = nullptr;
	vec3* faces = nullptr;
	Indices3* faceVerts = nullptr;
	TriEdges* faceEdges = nullptr;
};

#endif