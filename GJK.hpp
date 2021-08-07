/*
Author: Dan Rehberg
Date Modified: 8/5/2020

Purpose: This is a 3D implementation of the distance GJK algorithm.
	It does not include the case of orientation changes.
	Note, if anyone wishes to reused this, feel free to exchange the 
		translation only information for the Minkowski Difference to
		use 4x4 matrics for linear and rotational changes, or to combine
		quaternions and translations separately to limit floating point
		operations.
	Additionally, there might be redundancies in this implementation for
		testing out cases.
		The distance version of GJK is a little finnicky with exit conditions
			in 3D (compared to the boolean version of the algorithm).
		So, feel free to verify the conditions and alter them as best fits.
	Moreover, this build does not always function correctly.
		It might conservately assume a separation too often.
	Note, if this is reused, consider building a better Shape class as the 
		one for this program is capped in its mesh complexity.
*/

#ifndef __DISTANCE_GJK__
#define __DISTANCE_GJK__

#include <cstdint>
#include <utility>
#include "VectorMath.hpp"
#include "Shape.hpp"

struct simplex
{
	vec3 verts[4];//should initialze to zero vectors... see above
	int aID[4] = { 0, 0, 0, 0 };
	int bID[4] = { 0, 0, 0, 0 };
	//aternatively : a,b,c,d
	uint32_t count = 0;
};

int dirSign(float val)
{
	if (val < -0.0001)return -1;
	return 1;
}

void simplexMin(simplex& S, const vec3& P)
{
	//reduces 3d simplex to 2d
	if (S.count == 4)
	{
		vec3 ABC = closestPoint(P, S.verts[3], S.verts[2], S.verts[1]),
			ABD = closestPoint(P, S.verts[3], S.verts[2], S.verts[0]),
			BCD = closestPoint(P, S.verts[2], S.verts[1], S.verts[0]),
			ACD = closestPoint(P, S.verts[3], S.verts[1], S.verts[0]);
		if (dot(ABC, ABC) < dot(ABD, ABD))
		{
			if (dot(ABC, ABC) < dot(BCD, BCD))
			{
				S.verts[0] = S.verts[3];
				//S.inds[0] = S.inds[3];
				S.aID[0] = S.aID[3];
				S.bID[0] = S.bID[3];
			}
		}
		else
		{
			if (dot(ABD, ABD) < dot(BCD, BCD))
			{
				S.verts[1] = S.verts[0];
				S.verts[0] = S.verts[3];
				//S.inds[1] = S.inds[0];
				S.aID[1] = S.aID[0];
				S.bID[1] = S.bID[0];
				//S.inds[0] = S.inds[3];
				S.aID[0] = S.aID[3];
				S.bID[0] = S.bID[3];
			}
		}
		S.count = 3;
	}
}

float gjkDistance(const Shape& A, const Shape& B, const vec3& bOffset)
{
	vec3 D(1.0f, 0.25f, 0.5f);
	uint32_t supportA = A.supportPoint((-1.0f * D)),
		supportB = B.supportPoint(D);
	vec3 minkowskiDifference = (B.getVertex(supportB) + bOffset) -
		A.getVertex(supportA);
	simplex S;
	S.verts[0] = minkowskiDifference;
	//S.inds[0].first = supportA;
	//S.inds[0].second = supportB;
	S.aID[0] = supportA;
	S.bID[0] = supportB;
	D = (D * -1.0f);
	//get line segment -- I.E. 1D simplex
	supportA = A.supportPoint((-1.0f * D));
	supportB = B.supportPoint(D);
	minkowskiDifference = (B.getVertex(supportB) + bOffset) -
		A.getVertex(supportA);
	S.verts[1] = minkowskiDifference;
	//S.inds[1].first = supportA;
	//S.inds[1].second = supportB;
	S.aID[1] = supportA;
	S.bID[1] = supportB;
	S.count = 2;
	D = closestPoint(vec3(0.0f, 0.0f, 0.0f), S.verts[1], S.verts[0]);
	bool intersection = false;
	bool end = false;
	while (!end)
	{
		D = -1.0f * D;
		if (D.x == 0.0f && D.y == 0.0f && D.z == 0.0f)
		{
			intersection = true;
			end = true;
			break;
		}
		supportA = A.supportPoint((-1.0f * D));
		supportB = B.supportPoint(D);
		minkowskiDifference = (B.getVertex(supportB) + bOffset) -
			A.getVertex(supportA);
		switch (S.count)
		{
		case 2:
		{
			if (S.aID[0] == supportA && S.bID[0] == supportB)
			{
				intersection = false;
				end = true;
				break;
			}
			if (S.aID[1] == supportA && S.bID[1] == supportB)
			{
				intersection = false;
				end = true;
				break;
			}
			S.verts[2] = minkowskiDifference;
			//S.inds[2].first = supportA;
			//S.inds[2].second = supportB;
			S.aID[2] = supportA;
			S.bID[2] = supportB;
			S.count = 3;
			break;
		}
		case 3:
		{
			if (S.aID[0] == supportA && S.bID[0] == supportB)
			{
				intersection = false;
				end = true;
				break;
			}
			if (S.aID[1] == supportA && S.bID[1] == supportB)
			{
				intersection = false;
				end = true;
				break;
			}
			if (S.aID[2] == supportA && S.bID[2] == supportB)
			{
				intersection = false;
				end = true;
				break;
			}
			S.verts[3] = minkowskiDifference;
			//S.inds[3].first = supportA;
			//S.inds[3].second = supportB;
			S.aID[3] = supportA;
			S.bID[3] = supportB;
			S.count = 4;
			simplexMin(S, vec3(0.0f, 0.0f, 0.0f));
			break;
		}
		default:break;//bad dimensions
		}
		if (end)
		{
			//repeating indices, likely no intersection
			break;
		}
		vec3 tempD;
		switch (S.count)
		{
		case 2:
		{
			tempD = closestPoint(vec3(0.0f, 0.0f, 0.0f), S.verts[1], S.verts[0]);
			break;
		}
		case 3:
		{
			tempD = closestPoint(vec3(0.0f, 0.0f, 0.0f), S.verts[2], S.verts[1], S.verts[0]);
			break;
		}
		default: end = true; intersection = true; break;//intersection likely
		}
		if (dot(D, D) <= dot(tempD, tempD))
		{
			//shortest distance probably found, likely no intersection
			end = true;
			intersection = false;
			break;
		}
		else
		{
			D = tempD;
		}
	}
	if (intersection)
	{
		return -1.0f;
	}
	else
	{
		return std::sqrt(dot(D, D));
	}
}

std::pair<float, float> gjkToI(const Shape& A, const Shape& B, const vec3& bOffset, vec3 bVelocity)
{
	std::pair<float, float> result(std::pair<float, float>(-1.0f, -1.0f));
	float start = 0.0f, end = 1.0f, current;
	//Assuming 1 arbitrary time unit traveled
	float intersection = 0.01f;
	float distance = 1.0f;
	unsigned int maximumItr = 0; //In case a valid (i.e. will intersect) case is provided, stop after this many iterations
	while (maximumItr < 1000)
	{
		current = (start + end) * 0.5f;
		distance = gjkDistance(A, B, (bVelocity * current) + bOffset);
		if (distance < -0.5)
		{
			//intersection occurred
			end = current;
			start = 0.0f;
		}
		else if (distance < intersection)
		{
			//found the ToI
			result.first = current;
			result.second = distance;
			break;
		}
		else
		{
			//Not within an applicable distance yet
			start = current;
		}
		++maximumItr;
	}

	return result;
}

#endif