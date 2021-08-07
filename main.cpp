/*
Author: Dan Rehberg
Advisor: Dr. Francisco Ortega
Class: CS498 Independent Research

Date: 8/4/2021 - Final Code and Trials

Purpose: This research is meant to be both an exercise in mathematics and parallel computing,
		while also considering the continued growth of consumer CPU with strong parallel capabilities.
	Notably, while parallelism is not new, the features of CPU design allow for different execution of
		parallelism compared to batching for large problems on a GPU.
	For moderately sized problems, a single application can utilize similar parallel dispatch 
		designs to minimize the execution time of a complex task.
	The task focused on is Collision Detection, with a goal to see why a parallel system could
		be useful for this type of problem.
		~The reason found, is to solve for the time of intersection exclusively for 
			continuous collision detection systems.
	This is being compared to an optimal serial algorithm - GJK - for both standard proximity testing
		and in an additionally iterative case for time of intersection solving.
	The nature of design in this collision detection is to use thread execution in place of iterations.
		However, for some of the problem sizes, threads will be iterating over work items.
			~When the problem size is greater than thread count.
*/

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <utility>
#include "ThreadPool.hpp"
#include "GJK.hpp"
#include "VectorMath.hpp"
#include "Meshes.hpp"
#include "Shape.hpp"

//Hello World style Parallel Task - Dot product
float vectorA[14] = { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 0.12f, 0.13f, 0.14f, 0.15f, 0.16f };
float vectorB[14] = { 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 0.12f, 0.13f, 0.14f, 0.15f, 0.16f, 1.0f };
float vectorResult = 0.0f;
void fourteenDotProduct(std::mutex& m, unsigned int index)
{
	float component = vectorA[index] * vectorB[index];
	{
		std::lock_guard<std::mutex> lock(m);
		vectorResult += component;
	}
}

std::vector<float> nVectorA;
std::vector<float> nVectorB;
std::atomic_uint32_t atomicResult;
void nDotProduct(std::mutex& m, unsigned int index)
{
	float component = nVectorA[index] * nVectorB[index];
	{
		//std::lock_guard<std::mutex> lock(m);
		//vectorResult += component;
	}
	//Alternative
	//uint32_t component = static_cast<uint32_t>(nVectorA[index] * nVectorB[index] * 100000.0f);
	//atomicResult.fetch_add(component, std::memory_order_relaxed);
}

vec3 translation, velocity;

std::atomic_uint32_t entry;
float projTimes[242];
void hyperplaneCubeAllFacesToI(std::mutex& m, unsigned int index)
{
	float min = 1000.0f;
	vec3 displaced = cube.vertices[index] + translation;
	for (int i = 0; i < cube.faceCount; ++i)
	{
		float normalVelocity = dot(cube.faces[i], normalize(velocity));
		if (normalVelocity < 0.000001f)continue;
		//Moller-Trumbore ray triangle intersection
		const vec3& A = cube.faceEdges[i].edge[0];
		const vec3& B = cube.faceEdges[i].edge[1];
		vec3 h = cross(velocity, B);
		float a = dot(A, h);
		if (a > -0.00001f && a < 0.00001f)//If a is essentially zero
			continue;
		float f = 1.0f / a;
		vec3 s = displaced - cube.vertices[cube.faceVerts[i].ind[0]];
		float u = f * dot(s, h);
		if (u < 0.0f || u > 1.0f) //Not within the first Barycentric coordinate
			continue;
		vec3 q = cross(s, A);
		float v = f * dot(velocity, q);
		if (v < 0.0f || (u + v) > 1.0f) //Not within Barycentric bounds
			continue;
		float t = f * dot(B, q);
		if (t > 0.000001f && t < 1.00001f)
		{
			//valid ToI
			if (t < min) min = t;
		}

	}
	projTimes[index] = min;
}

void hyperplaneSuzanneAllFacesToI(std::mutex& m, unsigned int index)
{
	float min = 1000.0f;
	vec3 displaced = suzanne.vertices[index] + translation;
	for (int i = 0; i < suzanne.faceCount; ++i)
	{
		float normalVelocity = dot(suzanne.faces[i], normalize(velocity));
		if (normalVelocity < 0.000001f)continue;
		//Moller-Trumbore ray triangle intersection
		const vec3& A = suzanne.faceEdges[i].edge[0];
		const vec3& B = suzanne.faceEdges[i].edge[1];
		vec3 h = cross(velocity, B);
		float a = dot(A, h);
		if (a > -0.00001f && a < 0.00001f)//If a is essentially zero
			continue;
		float f = 1.0f / a;
		vec3 s = displaced - suzanne.vertices[suzanne.faceVerts[i].ind[0]];
		float u = f * dot(s, h);
		if (u < 0.0f || u > 1.0f) //Not within the first Barycentric coordinate
			continue;
		vec3 q = cross(s, A);
		float v = f * dot(velocity, q);
		if (v < 0.0f || (u + v) > 1.0f) //Not within Barycentric bounds
			continue;
		float t = f * dot(B, q);
		if (t > 0.000001f && t < 1.00001f)
		{
			//valid ToI
			if (t < min) min = t;
		}

	}
	projTimes[index] = min;
}

void hyperplaneSphereAllFacesToI(std::mutex& m, unsigned int index)
{
	float min = 1000.0f;
	vec3 displaced = sphere.vertices[index] + translation;
	for (int i = 0; i < sphere.faceCount; ++i)
	{
		float normalVelocity = dot(sphere.faces[i], normalize(velocity));
		if (normalVelocity < 0.0001f)continue;
		//Moller-Trumbore ray triangle intersection
		const vec3& A = sphere.faceEdges[i].edge[0];
		const vec3& B = sphere.faceEdges[i].edge[1];
		vec3 h = cross(velocity, B);
		float a = dot(A, h);
		if (a > -0.00001f && a < 0.00001f)//If a is essentially zero
			continue;
		float f = 1.0f / a;
		vec3 s = displaced - sphere.vertices[sphere.faceVerts[i].ind[0]];
		float u = f * dot(s, h);
		if (u < 0.0f || u > 1.0f) //Not within the first Barycentric coordinate
			continue;
		vec3 q = cross(s, A);
		float v = f * dot(velocity, q);
		if (v < 0.0f || (u + v) > 1.0f) //Not within Barycentric bounds
			continue;
		float t = f * dot(B, q);
		if (t > 0.000001f && t < 1.00001f)
		{
			//valid ToI
			if (t < min) min = t;
		}
		
	}
	projTimes[index] = min;
}

bool validFaces[480];
void hyperplaneCullSuzanne(std::mutex& m, unsigned int index)
{
	for (unsigned int i = 0; i < 8; ++i)
	{
		float normalVelocity = dot(suzanne.faces[index], velocity);
		if (normalVelocity < 0.0001f)continue;
		validFaces[index] = true;
		break;
	}
}

void hyperplaneCullSphere(std::mutex& m, unsigned int index)
{
	for (unsigned int i = 0; i < 8; ++i)
	{
		float normalVelocity = dot(sphere.faces[index], velocity);
		if (normalVelocity < 0.0001f)continue;
		validFaces[index] = true;
		break;
	}
}

int testFaces[480];
int testCount = 0;

void hyperplaneReducedSuzanne(std::mutex& m, unsigned int index)
{
	float min = 1000.0f;
	vec3 displaced = suzanne.vertices[index] + translation;
	for (int i = 0; i < testCount; ++i)
	{
		float normalVelocity = dot(suzanne.faces[testFaces[i]], normalize(velocity));
		if (normalVelocity < 0.0001f)continue;
		//Moller-Trumbore ray triangle intersection
		const vec3& A = suzanne.faceEdges[testFaces[i]].edge[0];
		const vec3& B = suzanne.faceEdges[testFaces[i]].edge[1];
		vec3 h = cross(velocity, B);
		float a = dot(A, h);
		if (a > -0.00001f && a < 0.00001f)//If a is essentially zero
			continue;
		float f = 1.0f / a;
		vec3 s = displaced - suzanne.vertices[suzanne.faceVerts[testFaces[i]].ind[0]];
		float u = f * dot(s, h);
		if (u < 0.0f || u > 1.0f) //Not within the first Barycentric coordinate
			continue;
		vec3 q = cross(s, A);
		float v = f * dot(velocity, q);
		if (v < 0.0f || (u + v) > 1.0f) //Not within Barycentric bounds
			continue;
		float t = f * dot(B, q);
		if (t > 0.000001f && t < 1.00001f)
		{
			//valid ToI
			if (t < min) min = t;
		}

	}
	projTimes[index] = min;
}

void hyperplaneReducedSphere(std::mutex& m, unsigned int index)
{
	float min = 1000.0f;
	vec3 displaced = sphere.vertices[index] + translation;
	for (int i = 0; i < testCount; ++i)
	{
		float normalVelocity = dot(sphere.faces[testFaces[i]], normalize(velocity));
		if (normalVelocity < 0.0001f)continue;
		//Moller-Trumbore ray triangle intersection
		const vec3& A = sphere.faceEdges[testFaces[i]].edge[0];
		const vec3& B = sphere.faceEdges[testFaces[i]].edge[1];
		vec3 h = cross(velocity, B);
		float a = dot(A, h);
		if (a > -0.00001f && a < 0.00001f)//If a is essentially zero
			continue;
		float f = 1.0f / a;
		vec3 s = displaced - sphere.vertices[sphere.faceVerts[testFaces[i]].ind[0]];
		float u = f * dot(s, h);
		if (u < 0.0f || u > 1.0f) //Not within the first Barycentric coordinate
			continue;
		vec3 q = cross(s, A);
		float v = f * dot(velocity, q);
		if (v < 0.0f || (u + v) > 1.0f) //Not within Barycentric bounds
			continue;
		float t = f * dot(B, q);
		if (t > 0.000001f && t < 1.00001f)
		{
			//valid ToI
			if (t < min) min = t;
		}

	}
	projTimes[index] = min;
}

int main()
{
	//Create an instance of a Thread Pool, give it the optimal number of available threads minus 1 (to avoid context switching with main thread)
	ThreadPool pool(std::thread::hardware_concurrency() - 1);
	//wait for thread pool to initialize
	pool.initialized();

	//Test the dot product result a few times...
	std::chrono::time_point<std::chrono::steady_clock> startTime;
	uint32_t delta = 0;
	float serialVectorResult = 0.0f;
	startTime = std::chrono::steady_clock::now();
	for (unsigned int i = 0; i < 14; ++i)
	{
		serialVectorResult += vectorA[i] * vectorB[i];
	}
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "serial 14 component dot product: " << delta << '\n';

	for (unsigned int i = 0; i < 10000; ++i)
	{
		//Give the task count and task function to the dispatch call
		if (i == 9999)
		{
			startTime = std::chrono::steady_clock::now();
			pool.dispatch(14, &fourteenDotProduct);
			delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
			std::cout << "last execution of 14 dot product: " << delta << '\n';
		}
		else pool.dispatch(14, &fourteenDotProduct);
		//compare to verify the result is the same to the serial method
		//std::cout << ((serialVectorResult - vectorResult < 0.0001f) ? std::to_string(i) + " is equal to serial\n" : std::to_string(i) + " is wrong!\n");
		//std::cout << "serial: " << serialVectorResult << " parallel: " << vectorResult << "\n";
		//reset the result from the dispatch to dispatch again for additional comparison tests
		vectorResult = 0.0f;
	}

	std::minstd_rand0 nextVal;

	unsigned int quantity;
	std::cin >> quantity;

	for (unsigned int i = 1; i < (quantity + 1); ++i)
	{
		std::cout << "nextVal: " << nextVal();
		nVectorA.push_back(static_cast<float>(nextVal() % 1000) * 0.001f);
		nVectorB.push_back(static_cast<float>(nextVal() % 1000) * 0.001f);
	}
	std::cout << "Finished building vectors\n";
	serialVectorResult = 0.0f;
	startTime = std::chrono::steady_clock::now();
	for (unsigned int i = 0; i < quantity; ++i)
	{
		float component = nVectorA[i] * nVectorB[i];
		serialVectorResult += component;
	}
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "serial 1000 component dot product: " << delta << '\n';

	for (unsigned int i = 0; i < 1000; ++i)
	{
		//Give the task count and task function to the dispatch call
		if (i == 999)
		{
			startTime = std::chrono::steady_clock::now();
			pool.dispatch(quantity, &nDotProduct);
			delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
			std::cout << "last execution 1000 dot product: " << delta << '\n';
			std::cout << "serial: " << serialVectorResult << " parallel: " << (0.00001f * static_cast<float>(atomicResult.load(std::memory_order_relaxed))) << '\n';
		}
		else pool.dispatch(quantity, &nDotProduct);
		//compare to verify the result is the same to the serial method
		//std::cout << ((serialVectorResult - vectorResult < 0.0001f) ? std::to_string(i) + " is equal to serial\n" : std::to_string(i) + " is wrong!\n");
		//std::cout << "serial: " << serialVectorResult << " parallel: " << vectorResult << "\n";
		//reset the result from the dispatch to dispatch again for additional comparison tests
		//vectorResult = 0.0f;
		atomicResult.store(0, std::memory_order_relaxed);
	}

	//Geometric testing begins
	initShapes();
	buildFaces();

	//GJK
	translation = vec3(5.0f, 0.0f, 0.0f);
	float distance = 0.0f;
	startTime = std::chrono::steady_clock::now();
	distance = gjkDistance(cube, cube, translation);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Cube to Cube distance is: " << distance << " and took: " << delta << " microseconds\n";
	startTime = std::chrono::steady_clock::now();
	distance = gjkDistance(suzanne, suzanne, translation);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Suzanne to Suzanne distance is: " << distance << " and took: " << delta << " microseconds\n";
	startTime = std::chrono::steady_clock::now();
	distance = gjkDistance(sphere, sphere, translation);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Sphere to Sphere distance is: " << distance << " and took: " << delta << " microseconds\n";
	
	//ToI GJK
	velocity = vec3(-5.0f, 0.0f, 0.0f);
	std::pair<float, float> timeDistance;
	startTime = std::chrono::steady_clock::now();
	timeDistance = gjkToI(cube, cube, translation, velocity);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Cube to Cube Time is: " << timeDistance.first << " and took: " << delta << " microseconds\n";
	startTime = std::chrono::steady_clock::now();
	timeDistance = gjkToI(suzanne, suzanne, translation, velocity);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Suzanne to Suzanne Time is: " << timeDistance.first << " and took: " << delta << " microseconds\n";
	startTime = std::chrono::steady_clock::now();
	timeDistance = gjkToI(sphere, sphere, translation, velocity);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Sphere to Sphere Time is: " << timeDistance.first << " and took: " << delta << " microseconds\n";

	//Projection Timing Tests
	//Box
	startTime = std::chrono::steady_clock::now();
	pool.dispatch(8, &hyperplaneCubeAllFacesToI);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Cube all faces tested: " << delta << " microseconds\n";
	//Suzanne
	startTime = std::chrono::steady_clock::now();
	pool.dispatch(66, &hyperplaneSuzanneAllFacesToI);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Suzanne all faces tested: " << delta << " microseconds\n";
	//Sphere
	startTime = std::chrono::steady_clock::now();
	pool.dispatch(242, &hyperplaneSphereAllFacesToI);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Sphere all faces tested: " << delta << " microseconds\n";
	//Culling Performance
	//	Suzanne
	startTime = std::chrono::steady_clock::now();
	pool.dispatch(128, &hyperplaneCullSuzanne);

	for (int i = 0; i < 128; ++i)
	{
		if (validFaces[i])
		{
			testFaces[testCount++] = i;
			validFaces[i] = false;
		}
	}
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Suzanne Culling time: " << delta << " microseconds\n";

	startTime = std::chrono::steady_clock::now();
	pool.dispatch(66, &hyperplaneReducedSuzanne);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Suzanne Reduced time: " << delta << " microseconds\n";
	//	Sphere
	testCount = 0;
	startTime = std::chrono::steady_clock::now();
	pool.dispatch(480, &hyperplaneCullSphere);

	for (int i = 0; i < 480; ++i)
	{
		if (validFaces[i])
		{
			testFaces[testCount++] = i;
			validFaces[i] = false;
		}
	}
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Sphere Culling time: " << delta << " microseconds\n";

	startTime = std::chrono::steady_clock::now();
	pool.dispatch(242, &hyperplaneReducedSphere);
	delta = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count();
	std::cout << "Sphere Reduced time: " << delta << " microseconds\n";

	char c;
	std::cin >> c;

	return 0;
}