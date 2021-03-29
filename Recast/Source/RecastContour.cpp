//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "Recast.h"
#include "RecastAlloc.h"
#include "RecastAssert.h"


static int getCornerHeight(int x, int y, int i, int dir,
						   const rcCompactHeightfield& chf,
						   bool& isBorderVertex)
{
	const rcCompactSpan& s = chf.spans[i];
	int ch = (int)s.y; // 用于保存四个采样点里的最大高度
	int dirp = (dir+1) & 0x3;

	// 四个方向
	// 23
	// 10
	//
	// 12
	// 03
	//
	// 01
	// 32
	//
	// 30
	// 21

	unsigned int regs[4] = {0,0,0,0};
	
	// Combine region and area codes in order to prevent
	// border vertices which are in between two areas to be removed.
	regs[0] = chf.spans[i].reg | (chf.areas[i] << 16);
	
	if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
	{
		const int ax = x + rcGetDirOffsetX(dir);
		const int ay = y + rcGetDirOffsetY(dir);
		const int ai = (int)chf.cells[ax+ay*chf.width].index + rcGetCon(s, dir);
		const rcCompactSpan& as = chf.spans[ai];
		ch = rcMax(ch, (int)as.y);
		regs[1] = chf.spans[ai].reg | (chf.areas[ai] << 16);
		if (rcGetCon(as, dirp) != RC_NOT_CONNECTED)
		{
			const int ax2 = ax + rcGetDirOffsetX(dirp);
			const int ay2 = ay + rcGetDirOffsetY(dirp);
			const int ai2 = (int)chf.cells[ax2+ay2*chf.width].index + rcGetCon(as, dirp);
			const rcCompactSpan& as2 = chf.spans[ai2];
			ch = rcMax(ch, (int)as2.y);
			regs[2] = chf.spans[ai2].reg | (chf.areas[ai2] << 16);
		}
	}
	if (rcGetCon(s, dirp) != RC_NOT_CONNECTED)
	{
		const int ax = x + rcGetDirOffsetX(dirp);
		const int ay = y + rcGetDirOffsetY(dirp);
		const int ai = (int)chf.cells[ax+ay*chf.width].index + rcGetCon(s, dirp);
		const rcCompactSpan& as = chf.spans[ai];
		ch = rcMax(ch, (int)as.y);
		regs[3] = chf.spans[ai].reg | (chf.areas[ai] << 16);
		if (rcGetCon(as, dir) != RC_NOT_CONNECTED)
		{
			const int ax2 = ax + rcGetDirOffsetX(dir);
			const int ay2 = ay + rcGetDirOffsetY(dir);
			const int ai2 = (int)chf.cells[ax2+ay2*chf.width].index + rcGetCon(as, dir);
			const rcCompactSpan& as2 = chf.spans[ai2];
			ch = rcMax(ch, (int)as2.y);
			regs[2] = chf.spans[ai2].reg | (chf.areas[ai2] << 16);
		}
	}

	// Check if the vertex is special edge vertex, these vertices will be removed later.
	// 注意这里检查的四个点，构成了一个正方形
	for (int j = 0; j < 4; ++j)
	{
		const int a = j;
		const int b = (j+1) & 0x3;
		const int c = (j+2) & 0x3;
		const int d = (j+3) & 0x3;

		// The vertex is a border vertex there are two same exterior cells in a row,
		// followed by two interior cells and none of the regions are out of bounds.
		// 判断 a b 两个点，是否是 tile 边界上相同 region、相同 area 的区域的点
		// 见 rcBuildRegions->paintRectRegion
		const bool twoSameExts = (regs[a] & regs[b] & RC_BORDER_REG) != 0 && regs[a] == regs[b];
		// 判断对面的两个点 c d，是否均不是边界点
		const bool twoInts = ((regs[c] | regs[d]) & RC_BORDER_REG) == 0;
		// 判断对面的两个点 c d，是否属于相同的 area
		// area 是用于用户自定义地形类别的，寻路时也可能用于计算不同的花费
		// 所以这里是：即便处于相同的 region，但是如果 area 不同，依旧算是边缘
		// c 和 d 必须分属相同的 area，但是 region 可以不同。。。
		const bool intsSameArea = (regs[c]>>16) == (regs[d]>>16);
		// TODO 四个点均在合法 region 内
		const bool noZeros = regs[a] != 0 && regs[b] != 0 && regs[c] != 0 && regs[d] != 0;
		if (twoSameExts && twoInts && intsSameArea && noZeros)
		{
			isBorderVertex = true;
			break;
		}
	}
	
	return ch;
}

static void walkContour(int x, int y, int i,
						rcCompactHeightfield& chf,
						unsigned char* flags, rcIntArray& points)
{
	// Choose the first non-connected edge
	// 按左、上、右、下的顺序找到第一个未连通边的方向
	unsigned char dir = 0;
	while ((flags[i] & (1 << dir)) == 0)
		dir++;

	unsigned char startDir = dir;
	int starti = i;

	const unsigned char area = chf.areas[i];
	
	int iter = 0;
	while (++iter < 40000)
	{
		// 注意这里为真代表这个方向上没有邻接 span，或者邻接的 span 不属于同一个region
		if (flags[i] & (1 << dir))
		{
			// Choose the edge corner
			bool isBorderVertex = false;
			bool isAreaBorder = false;
			int px = x;

			// CornerHeight 是 当前span、左、左下、下 四个方向里的最大高度
			// ?为什么只计算左下方？
			// isBorderVertex，是否是边缘上的节点
			// 含义是对于组成一个正方形的相邻四个 span，如果横平竖直方向上，两两的连线正好分属两个不同的区域
			// 则可以认为这几个点正好在边缘上
			// A B
			// C D
			// AB|CD
			// BD|CA
			// DC|AB
			// CA|BD
			int py = getCornerHeight(x, y, i, dir, chf, isBorderVertex);
			int pz = y;
			switch(dir)
			{
				// 这里 px pz 是格子的 x y 坐标
				case 0: pz++; break; // 左方不连通，获取下方的格子坐标
				case 1: px++; pz++; break; // 上方不连通，获取右下方的格子坐标？
				case 2: px++; break; // 右方不连通，就获取右方格子的坐标？？？
			}

			int r = 0;
			const rcCompactSpan& s = chf.spans[i];
			// 一个 span 在某个方向上不连通，这不算是 Area Border 的性质之一？
			if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
			{
				const int ax = x + rcGetDirOffsetX(dir);
				const int ay = y + rcGetDirOffsetY(dir);
				const int ai = (int)chf.cells[ax+ay*chf.width].index + rcGetCon(s, dir);
				r = (int)chf.spans[ai].reg;

				// 如果当前 span 指定方向上的邻接 span 不属于同一个 region，则当前 span 是 region 边缘
				if (area != chf.areas[ai])
					isAreaBorder = true;
			}
			else
			{
				int abc = 100;
			}
			if (isBorderVertex)
				r |= RC_BORDER_VERTEX; // 65536 - 0x1 0000
			if (isAreaBorder)
				r |= RC_AREA_BORDER; // 131072 - 0x2 0000
			points.push(px);
			points.push(py);
			points.push(pz);
			points.push(r);

			flags[i] &= ~(1 << dir); // Remove visited edges
			// 顺时针转向
			dir = (dir+1) & 0x3;  // Rotate CW
		}
		else
		{
			// 如果对应方向有邻接 span
			// 则从下一个 span 开始，逆时针转向
			int ni = -1;
			const int nx = x + rcGetDirOffsetX(dir);
			const int ny = y + rcGetDirOffsetY(dir);
			const rcCompactSpan& s = chf.spans[i];
			if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
			{
				const rcCompactCell& nc = chf.cells[nx+ny*chf.width];
				ni = (int)nc.index + rcGetCon(s, dir);
			}
			if (ni == -1)
			{
				// Should not happen.
				return;
			}
			x = nx;
			y = ny;
			i = ni;
			// 逆时针转向
			dir = (dir+3) & 0x3;	// Rotate CCW
		}

		// 绕回到起始点，可以结束循环
		// 注意这里绕回来的方向和起始方向不一致的话，需要继续循环
		if (starti == i && startDir == dir)
		{
			break;
		}
	}
}

static float distancePtSeg(const int x, const int z,
						   const int px, const int pz,
						   const int qx, const int qz)
{
	float pqx = (float)(qx - px);
	float pqz = (float)(qz - pz);
	float dx = (float)(x - px);
	float dz = (float)(z - pz);
	float d = pqx*pqx + pqz*pqz;
	float t = pqx*dx + pqz*dz;
	if (d > 0)
		t /= d;
	if (t < 0)
		t = 0;
	else if (t > 1)
		t = 1;
	
	dx = px + t*pqx - x;
	dz = pz + t*pqz - z;
	
	return dx*dx + dz*dz;
}

static void simplifyContour(rcIntArray& points, rcIntArray& simplified,
							const float maxError, const int maxEdgeLen, const int buildFlags)
{
	// Add initial points.
	// 这里遍历轮廓边缘格子，从 r 中取出低 16 位里的 region id
	// 这里的 points[i+3]，低 16 位是 region_id，17、18 位是特殊标记位
	// 找到第一个 region id 不为 0 的点
	bool hasConnections = false;
	for (int i = 0; i < points.size(); i += 4)
	{
		if ((points[i+3] & RC_CONTOUR_REG_MASK) != 0)
		{
			hasConnections = true;
			break;
		}
		else
		{
			int b = 0;
		}
	}
	
	if (hasConnections)
	{
		// The contour has some portals to other regions.
		// Add a new point to every location where the region changes.
		for (int i = 0, ni = points.size()/4; i < ni; ++i)
		{
			// ii 是轮廓连线上，当前点 i 的下一个点的序号
			int ii = (i+1) % ni;
			// 判断当前点与下一个点是否属于同一个 region id
			const bool differentRegs = (points[i*4+3] & RC_CONTOUR_REG_MASK) != (points[ii*4+3] & RC_CONTOUR_REG_MASK);
			// 两个点一个是 area border 一个不是，证明当前是边界
			// ? 两个点都是边界的情况？
			const bool areaBorders = (points[i*4+3] & RC_AREA_BORDER) != (points[ii*4+3] & RC_AREA_BORDER);
			// 这里为了简化什么？
			// 条件取反一下，即为当 !differentRegs && !areaBorders 时，丢弃点
			// 所以如果两个点：属于同一个region && （都不是 border 点 || 都是 border 点）时，可以忽略掉当前点（可以视为和下一个点合并了）
			if (differentRegs || areaBorders)
			{
				simplified.push(points[i*4+0]);
				simplified.push(points[i*4+1]);
				simplified.push(points[i*4+2]);
				simplified.push(i);
			}
		}
	}
	
	if (simplified.size() == 0)
	{
		// If there is no connections at all,
		// create some initial points for the simplification process.
		// Find lower-left and upper-right vertices of the contour.
		int llx = points[0];
		int lly = points[1];
		int llz = points[2];
		int lli = 0;
		int urx = points[0];
		int ury = points[1];
		int urz = points[2];
		int uri = 0;
		for (int i = 0; i < points.size(); i += 4)
		{
			int x = points[i+0];
			int y = points[i+1];
			int z = points[i+2];
			if (x < llx || (x == llx && z < llz))
			{
				llx = x;
				lly = y;
				llz = z;
				lli = i/4;
			}
			if (x > urx || (x == urx && z > urz))
			{
				urx = x;
				ury = y;
				urz = z;
				uri = i/4;
			}
		}
		simplified.push(llx);
		simplified.push(lly);
		simplified.push(llz);
		simplified.push(lli);
		
		simplified.push(urx);
		simplified.push(ury);
		simplified.push(urz);
		simplified.push(uri);
	}
	
	// Add points until all raw points are within
	// error tolerance to the simplified shape.
	const int pn = points.size()/4;
	for (int i = 0; i < simplified.size()/4; )
	{
		int ii = (i+1) % (simplified.size()/4);
		
		int ax = simplified[i*4+0];
		int az = simplified[i*4+2];
		int ai = simplified[i*4+3];

		int bx = simplified[ii*4+0];
		int bz = simplified[ii*4+2];
		int bi = simplified[ii*4+3];

		// Find maximum deviation from the segment.
		float maxd = 0;
		int maxi = -1;
		// ci - current i
		// cinc - c inc
		// endi
		int ci, cinc, endi;

		// Traverse the segment in lexilogical order so that the
		// max deviation is calculated similarly when traversing
		// opposite segments.
		if (bx > ax || (bx == ax && bz > az))
		{
			cinc = 1;
			ci = (ai+cinc) % pn;
			endi = bi;
		}
		else
		{
			cinc = pn-1;
			ci = (bi+cinc) % pn;
			endi = ai;
			rcSwap(ax, bx);
			rcSwap(az, bz);
		}
		
		// Tessellate only outer edges or edges between areas.
		if ((points[ci*4+3] & RC_CONTOUR_REG_MASK) == 0 ||
			(points[ci*4+3] & RC_AREA_BORDER))
		{
			// 这里求出两个点连线上的点与连线的距离最大平方 maxd，以及距离最大的点 maxi
			while (ci != endi)
			{
				float d = distancePtSeg(points[ci*4+0], points[ci*4+2], ax, az, bx, bz);
				if (d > maxd)
				{
					maxd = d;
					maxi = ci;
				}
				ci = (ci+cinc) % pn;
			}
		}
		
		
		// If the max deviation is larger than accepted error,
		// add new point, else continue to next segment.
		// 如果距离超出最大误差距离，则在这个误差距离最大的地方插入一个中间点
		// 然后下一次循环以这个点为起点继续
		// 如果在误差距离以内，则直接从下一个点开始遍历
		if (maxi != -1 && maxd > (maxError*maxError))
		{
			// Add space for the new point.
			simplified.resize(simplified.size()+4);
			const int n = simplified.size()/4;
			// 将 i+1 后的点后移一个位置
			for (int j = n-1; j > i; --j)
			{
				simplified[j*4+0] = simplified[(j-1)*4+0];
				simplified[j*4+1] = simplified[(j-1)*4+1];
				simplified[j*4+2] = simplified[(j-1)*4+2];
				simplified[j*4+3] = simplified[(j-1)*4+3];
			}
			// Add the point.
			simplified[(i+1)*4+0] = points[maxi*4+0];
			simplified[(i+1)*4+1] = points[maxi*4+1];
			simplified[(i+1)*4+2] = points[maxi*4+2];
			simplified[(i+1)*4+3] = maxi;
		}
		else
		{
			++i;
		}
	}
	
	// Split too long edges.
	if (maxEdgeLen > 0 && (buildFlags & (RC_CONTOUR_TESS_WALL_EDGES|RC_CONTOUR_TESS_AREA_EDGES)) != 0)
	{
		for (int i = 0; i < simplified.size()/4; )
		{
			const int ii = (i+1) % (simplified.size()/4);
			
			const int ax = simplified[i*4+0];
			const int az = simplified[i*4+2];
			const int ai = simplified[i*4+3];
			
			const int bx = simplified[ii*4+0];
			const int bz = simplified[ii*4+2];
			const int bi = simplified[ii*4+3];
			
			// Find maximum deviation from the segment.
			int maxi = -1;
			int ci = (ai+1) % pn;
			
			// Tessellate only outer edges or edges between areas.
			bool tess = false;
			// Wall edges.
			if ((buildFlags & RC_CONTOUR_TESS_WALL_EDGES) && (points[ci*4+3] & RC_CONTOUR_REG_MASK) == 0)
				tess = true;
			// Edges between areas.
			if ((buildFlags & RC_CONTOUR_TESS_AREA_EDGES) && (points[ci*4+3] & RC_AREA_BORDER))
				tess = true;
			
			if (tess)
			{
				int dx = bx - ax;
				int dz = bz - az;
				if (dx*dx + dz*dz > maxEdgeLen*maxEdgeLen)
				{
					// Round based on the segments in lexilogical order so that the
					// max tesselation is consistent regardles in which direction
					// segments are traversed.
					const int n = bi < ai ? (bi+pn - ai) : (bi - ai);
					if (n > 1)
					{
						if (bx > ax || (bx == ax && bz > az))
							maxi = (ai + n/2) % pn;
						else
							maxi = (ai + (n+1)/2) % pn;
					}
				}
			}
			
			// If the max deviation is larger than accepted error,
			// add new point, else continue to next segment.
			if (maxi != -1)
			{
				// Add space for the new point.
				simplified.resize(simplified.size()+4);
				const int n = simplified.size()/4;
				for (int j = n-1; j > i; --j)
				{
					simplified[j*4+0] = simplified[(j-1)*4+0];
					simplified[j*4+1] = simplified[(j-1)*4+1];
					simplified[j*4+2] = simplified[(j-1)*4+2];
					simplified[j*4+3] = simplified[(j-1)*4+3];
				}
				// Add the point.
				simplified[(i+1)*4+0] = points[maxi*4+0];
				simplified[(i+1)*4+1] = points[maxi*4+1];
				simplified[(i+1)*4+2] = points[maxi*4+2];
				simplified[(i+1)*4+3] = maxi;
			}
			else
			{
				++i;
			}
		}
	}
	
	for (int i = 0; i < simplified.size()/4; ++i)
	{
		// The edge vertex flag is take from the current raw point,
		// and the neighbour region is take from the next raw point.
		const int ai = (simplified[i*4+3]+1) % pn;
		const int bi = simplified[i*4+3];
		simplified[i*4+3] = (points[ai*4+3] & (RC_CONTOUR_REG_MASK|RC_AREA_BORDER)) | (points[bi*4+3] & RC_BORDER_VERTEX);
	}
	
}

static int calcAreaOfPolygon2D(const int* verts, const int nverts)
{
	int area = 0;
	for (int i = 0, j = nverts-1; i < nverts; j=i++)
	{
		const int* vi = &verts[i*4];
		const int* vj = &verts[j*4];
		area += vi[0] * vj[2] - vj[0] * vi[2];
	}
	return (area+1) / 2;
}

// TODO: these are the same as in RecastMesh.cpp, consider using the same.
// Last time I checked the if version got compiled using cmov, which was a lot faster than module (with idiv).
inline int prev(int i, int n) { return i-1 >= 0 ? i-1 : n-1; }
inline int next(int i, int n) { return i+1 < n ? i+1 : 0; }

inline int area2(const int* a, const int* b, const int* c)
{
	return (b[0] - a[0]) * (c[2] - a[2]) - (c[0] - a[0]) * (b[2] - a[2]);
}

//	Exclusive or: true iff exactly one argument is true.
//	The arguments are negated to ensure that they are 0/1
//	values.  Then the bitwise Xor operator may apply.
//	(This idea is due to Michael Baldwin.)
inline bool xorb(bool x, bool y)
{
	return !x ^ !y;
}

// Returns true iff c is strictly to the left of the directed
// line through a to b.
inline bool left(const int* a, const int* b, const int* c)
{
	return area2(a, b, c) < 0;
}

inline bool leftOn(const int* a, const int* b, const int* c)
{
	return area2(a, b, c) <= 0;
}

inline bool collinear(const int* a, const int* b, const int* c)
{
	return area2(a, b, c) == 0;
}

//	Returns true iff ab properly intersects cd: they share
//	a point interior to both segments.  The properness of the
//	intersection is ensured by using strict leftness.
static bool intersectProp(const int* a, const int* b, const int* c, const int* d)
{
	// Eliminate improper cases.
	if (collinear(a,b,c) || collinear(a,b,d) ||
		collinear(c,d,a) || collinear(c,d,b))
		return false;
	
	return xorb(left(a,b,c), left(a,b,d)) && xorb(left(c,d,a), left(c,d,b));
}

// Returns T iff (a,b,c) are collinear and point c lies
// on the closed segement ab.
static bool between(const int* a, const int* b, const int* c)
{
	if (!collinear(a, b, c))
		return false;
	// If ab not vertical, check betweenness on x; else on y.
	if (a[0] != b[0])
		return	((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]));
	else
		return	((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]));
}

// Returns true iff segments ab and cd intersect, properly or improperly.
static bool intersect(const int* a, const int* b, const int* c, const int* d)
{
	if (intersectProp(a, b, c, d))
		return true;
	else if (between(a, b, c) || between(a, b, d) ||
			 between(c, d, a) || between(c, d, b))
		return true;
	else
		return false;
}

static bool vequal(const int* a, const int* b)
{
	return a[0] == b[0] && a[2] == b[2];
}

static bool intersectSegCountour(const int* d0, const int* d1, int i, int n, const int* verts)
{
	// For each edge (k,k+1) of P
	for (int k = 0; k < n; k++)
	{
		int k1 = next(k, n);
		// Skip edges incident to i.
		if (i == k || i == k1)
			continue;
		const int* p0 = &verts[k * 4];
		const int* p1 = &verts[k1 * 4];
		if (vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1))
			continue;
		
		if (intersect(d0, d1, p0, p1))
			return true;
	}
	return false;
}

static bool	inCone(int i, int n, const int* verts, const int* pj)
{
	const int* pi = &verts[i * 4];
	const int* pi1 = &verts[next(i, n) * 4];
	const int* pin1 = &verts[prev(i, n) * 4];
	
	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if (leftOn(pin1, pi, pi1))
		return left(pi, pj, pin1) && left(pj, pi, pi1);
	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1));
}


static void removeDegenerateSegments(rcIntArray& simplified)
{
	// Remove adjacent vertices which are equal on xz-plane,
	// or else the triangulator will get confused.
	int npts = simplified.size()/4;
	for (int i = 0; i < npts; ++i)
	{
		int ni = next(i, npts);
		
		if (vequal(&simplified[i*4], &simplified[ni*4]))
		{
			// Degenerate segment, remove.
			for (int j = i; j < simplified.size()/4-1; ++j)
			{
				simplified[j*4+0] = simplified[(j+1)*4+0];
				simplified[j*4+1] = simplified[(j+1)*4+1];
				simplified[j*4+2] = simplified[(j+1)*4+2];
				simplified[j*4+3] = simplified[(j+1)*4+3];
			}
			simplified.resize(simplified.size()-4);
			npts--;
		}
	}
}


static bool mergeContours(rcContour& ca, rcContour& cb, int ia, int ib)
{
	const int maxVerts = ca.nverts + cb.nverts + 2;
	int* verts = (int*)rcAlloc(sizeof(int)*maxVerts*4, RC_ALLOC_PERM);
	if (!verts)
		return false;
	
	int nv = 0;
	
	// Copy contour A.
	for (int i = 0; i <= ca.nverts; ++i)
	{
		int* dst = &verts[nv*4];
		const int* src = &ca.verts[((ia+i)%ca.nverts)*4];
		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
		dst[3] = src[3];
		nv++;
	}

	// Copy contour B
	for (int i = 0; i <= cb.nverts; ++i)
	{
		int* dst = &verts[nv*4];
		const int* src = &cb.verts[((ib+i)%cb.nverts)*4];
		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
		dst[3] = src[3];
		nv++;
	}
	
	rcFree(ca.verts);
	ca.verts = verts;
	ca.nverts = nv;
	
	rcFree(cb.verts);
	cb.verts = 0;
	cb.nverts = 0;
	
	return true;
}

struct rcContourHole
{
	rcContour* contour;
	int minx, minz, leftmost;
};

struct rcContourRegion
{
	rcContour* outline;
	rcContourHole* holes;
	int nholes;
};

struct rcPotentialDiagonal
{
	int vert;
	int dist;
};

// Finds the lowest leftmost vertex of a contour.
static void findLeftMostVertex(rcContour* contour, int* minx, int* minz, int* leftmost)
{
	*minx = contour->verts[0];
	*minz = contour->verts[2];
	*leftmost = 0;
	for (int i = 1; i < contour->nverts; i++)
	{
		const int x = contour->verts[i*4+0];
		const int z = contour->verts[i*4+2];
		if (x < *minx || (x == *minx && z < *minz))
		{
			*minx = x;
			*minz = z;
			*leftmost = i;
		}
	}
}

static int compareHoles(const void* va, const void* vb)
{
	const rcContourHole* a = (const rcContourHole*)va;
	const rcContourHole* b = (const rcContourHole*)vb;
	if (a->minx == b->minx)
	{
		if (a->minz < b->minz)
			return -1;
		if (a->minz > b->minz)
			return 1;
	}
	else
	{
		if (a->minx < b->minx)
			return -1;
		if (a->minx > b->minx)
			return 1;
	}
	return 0;
}


static int compareDiagDist(const void* va, const void* vb)
{
	const rcPotentialDiagonal* a = (const rcPotentialDiagonal*)va;
	const rcPotentialDiagonal* b = (const rcPotentialDiagonal*)vb;
	if (a->dist < b->dist)
		return -1;
	if (a->dist > b->dist)
		return 1;
	return 0;
}


static void mergeRegionHoles(rcContext* ctx, rcContourRegion& region)
{
	// Sort holes from left to right.
	for (int i = 0; i < region.nholes; i++)
		findLeftMostVertex(region.holes[i].contour, &region.holes[i].minx, &region.holes[i].minz, &region.holes[i].leftmost);
	
	qsort(region.holes, region.nholes, sizeof(rcContourHole), compareHoles);
	
	int maxVerts = region.outline->nverts;
	for (int i = 0; i < region.nholes; i++)
		maxVerts += region.holes[i].contour->nverts;
	
	rcScopedDelete<rcPotentialDiagonal> diags((rcPotentialDiagonal*)rcAlloc(sizeof(rcPotentialDiagonal)*maxVerts, RC_ALLOC_TEMP));
	if (!diags)
	{
		ctx->log(RC_LOG_WARNING, "mergeRegionHoles: Failed to allocated diags %d.", maxVerts);
		return;
	}
	
	rcContour* outline = region.outline;
	
	// Merge holes into the outline one by one.
	for (int i = 0; i < region.nholes; i++)
	{
		rcContour* hole = region.holes[i].contour;
		
		int index = -1;
		int bestVertex = region.holes[i].leftmost;
		for (int iter = 0; iter < hole->nverts; iter++)
		{
			// Find potential diagonals.
			// The 'best' vertex must be in the cone described by 3 cosequtive vertices of the outline.
			// ..o j-1
			//   |
			//   |   * best
			//   |
			// j o-----o j+1
			//         :
			int ndiags = 0;
			const int* corner = &hole->verts[bestVertex*4];
			for (int j = 0; j < outline->nverts; j++)
			{
				if (inCone(j, outline->nverts, outline->verts, corner))
				{
					int dx = outline->verts[j*4+0] - corner[0];
					int dz = outline->verts[j*4+2] - corner[2];
					diags[ndiags].vert = j;
					diags[ndiags].dist = dx*dx + dz*dz;
					ndiags++;
				}
			}
			// Sort potential diagonals by distance, we want to make the connection as short as possible.
			qsort(diags, ndiags, sizeof(rcPotentialDiagonal), compareDiagDist);
			
			// Find a diagonal that is not intersecting the outline not the remaining holes.
			index = -1;
			for (int j = 0; j < ndiags; j++)
			{
				const int* pt = &outline->verts[diags[j].vert*4];
				bool intersect = intersectSegCountour(pt, corner, diags[i].vert, outline->nverts, outline->verts);
				for (int k = i; k < region.nholes && !intersect; k++)
					intersect |= intersectSegCountour(pt, corner, -1, region.holes[k].contour->nverts, region.holes[k].contour->verts);
				if (!intersect)
				{
					index = diags[j].vert;
					break;
				}
			}
			// If found non-intersecting diagonal, stop looking.
			if (index != -1)
				break;
			// All the potential diagonals for the current vertex were intersecting, try next vertex.
			bestVertex = (bestVertex + 1) % hole->nverts;
		}
		
		if (index == -1)
		{
			ctx->log(RC_LOG_WARNING, "mergeHoles: Failed to find merge points for %p and %p.", region.outline, hole);
			continue;
		}
		if (!mergeContours(*region.outline, *hole, index, bestVertex))
		{
			ctx->log(RC_LOG_WARNING, "mergeHoles: Failed to merge contours %p and %p.", region.outline, hole);
			continue;
		}
	}
}


/// @par
///
/// The raw contours will match the region outlines exactly. The @p maxError and @p maxEdgeLen
/// parameters control how closely the simplified contours will match the raw contours.
///
/// Simplified contours are generated such that the vertices for portals between areas match up.
/// (They are considered mandatory vertices.)
///
/// Setting @p maxEdgeLength to zero will disabled the edge length feature.
///
/// See the #rcConfig documentation for more information on the configuration parameters.
///
/// @see rcAllocContourSet, rcCompactHeightfield, rcContourSet, rcConfig
bool rcBuildContours(rcContext* ctx, rcCompactHeightfield& chf,
					 const float maxError, const int maxEdgeLen,
					 rcContourSet& cset, const int buildFlags)
{
	rcAssert(ctx);
	
	const int w = chf.width;
	const int h = chf.height;
	const int borderSize = chf.borderSize;
	
	rcScopedTimer timer(ctx, RC_TIMER_BUILD_CONTOURS);
	
	rcVcopy(cset.bmin, chf.bmin);
	rcVcopy(cset.bmax, chf.bmax);
	if (borderSize > 0)
	{
		// borderSize 代表要在高度场包围盒边界处保留一定数量 cell 大小的边界区域
		// 构建轮廓时需要排除掉这个范围
		// If the heightfield was build with bordersize, remove the offset.
		const float pad = borderSize*chf.cs;
		cset.bmin[0] += pad;
		cset.bmin[2] += pad;
		cset.bmax[0] -= pad;
		cset.bmax[2] -= pad;
	}
	cset.cs = chf.cs;
	cset.ch = chf.ch;
	cset.width = chf.width - chf.borderSize*2;
	cset.height = chf.height - chf.borderSize*2;
	cset.borderSize = chf.borderSize;
	cset.maxError = maxError;
	
	int maxContours = rcMax((int)chf.maxRegions, 8);
	cset.conts = (rcContour*)rcAlloc(sizeof(rcContour)*maxContours, RC_ALLOC_PERM);
	if (!cset.conts)
		return false;
	cset.nconts = 0;
	
	rcScopedDelete<unsigned char> flags((unsigned char*)rcAlloc(sizeof(unsigned char)*chf.spanCount, RC_ALLOC_TEMP));
	if (!flags)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'flags' (%d).", chf.spanCount);
		return false;
	}
	
	ctx->startTimer(RC_TIMER_BUILD_CONTOURS_TRACE);
	
	// Mark boundaries.
	// 这里是将所有 span 四方向里不连通的边标记出来
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			const rcCompactCell& c = chf.cells[x+y*w];
			for (int i = (int)c.index, ni = (int)(c.index+c.count); i < ni; ++i)
			{
				unsigned char res = 0;
				const rcCompactSpan& s = chf.spans[i];
				if (!chf.spans[i].reg || (chf.spans[i].reg & RC_BORDER_REG))
				{
					// 如果不在任何区域内，或者是边界 region，则标记为 0，表示后面可以略过处理
					flags[i] = 0;
					continue;
				}

				for (int dir = 0; dir < 4; ++dir)
				{
					unsigned short r = 0;
					if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
					{
						const int ax = x + rcGetDirOffsetX(dir);
						const int ay = y + rcGetDirOffsetY(dir);
						const int ai = (int)chf.cells[ax+ay*w].index + rcGetCon(s, dir);
						r = chf.spans[ai].reg;
					}

					// 该方向上，当前 span 与邻接 span 的 region 相同，代表是连通的
					// 标记该方向为连通
					if (r == chf.spans[i].reg)
						res |= (1 << dir);
				}

				// 前面是连通标记，这里取反，得到不连通标记
				flags[i] = res ^ 0xf; // Inverse, mark non connected edges.
			}
		}
	}
	
	ctx->stopTimer(RC_TIMER_BUILD_CONTOURS_TRACE);
	
	rcIntArray verts(256); // 用于存放轮廓上顶点的临时数组
	rcIntArray simplified(64); // 用于存放平滑处理后的轮廓顶点的临时数组
	
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			const rcCompactCell& c = chf.cells[x+y*w];
			for (int i = (int)c.index, ni = (int)(c.index+c.count); i < ni; ++i)
			{
				// flags[i] == 0：无reg、边界reg、处于 region 内部位置（因为四方向邻接的都是相同 reg 的 span），或者几个方向均已处理完毕，这四种情况可以直接略过
				// flags[i] == 0xf：四方向均不连通或均为不同的 region，也可以直接略过
				// 因为这些情况里，这个 span 肯定不是轮廓线上的位置，对于构建轮廓线没有帮助
				if (flags[i] == 0 || flags[i] == 0xf)
				{
					flags[i] = 0;
					continue;
				}

				const unsigned short reg = chf.spans[i].reg;
				if (!reg || (reg & RC_BORDER_REG))
					continue;
				const unsigned char area = chf.areas[i];
				// 到这里， i 代表的一定是在 region 轮廓线上的 span

				verts.clear();
				simplified.clear();
				
				ctx->startTimer(RC_TIMER_BUILD_CONTOURS_TRACE);

				// walkContour: 绕着 region 边缘遍历一圈，把边缘的格子点坐标都保存到 verts 数组中
				// 数据格式为：x h y r
				// 分别对应格子 x 轴坐标、span 的边缘高度、格子的 y 坐标、region 标记（分段保存了 id、是否边缘节点、是否区域边缘）
				// ? RC_BORDER_VERTEX RC_AREA_BORDER RC_BORDER_REG 什么区别？
				// 猜想：（好像是错的）
				// RC_BORDER_VERTEX 是位于一条分界线上的点，另一边可能是 area 可能不是？
				// RC_AREA_BORDER 是明确的两个 area 之间的界限
				// RC_BORDER_REG 代表高度场里一个 span 属于 border，属于不可行走区域
				walkContour(x, y, i, chf, flags, verts);

				ctx->stopTimer(RC_TIMER_BUILD_CONTOURS_TRACE);

				ctx->startTimer(RC_TIMER_BUILD_CONTOURS_SIMPLIFY);
				// simplifyContour: 将区域轮廓进行简化处理
				// 主要是：
				// 1.在 maxError 范围内对折线做平滑处理
				// 2.对大于 maxEdgeLen 的线段做拆分处理
				simplifyContour(verts, simplified, maxError, maxEdgeLen, buildFlags);
				removeDegenerateSegments(simplified);

				ctx->stopTimer(RC_TIMER_BUILD_CONTOURS_SIMPLIFY);

				// Store region->contour remap info.
				// Create contour.
				if (simplified.size()/4 >= 3)
				{
					if (cset.nconts >= maxContours)
					{
						// Allocate more contours.
						// This happens when a region has holes.
						const int oldMax = maxContours;
						maxContours *= 2;
						rcContour* newConts = (rcContour*)rcAlloc(sizeof(rcContour)*maxContours, RC_ALLOC_PERM);
						for (int j = 0; j < cset.nconts; ++j)
						{
							newConts[j] = cset.conts[j];
							// Reset source pointers to prevent data deletion.
							cset.conts[j].verts = 0;
							cset.conts[j].rverts = 0;
						}
						rcFree(cset.conts);
						cset.conts = newConts;
						
						ctx->log(RC_LOG_WARNING, "rcBuildContours: Expanding max contours from %d to %d.", oldMax, maxContours);
					}
					
					rcContour* cont = &cset.conts[cset.nconts++];
					
					cont->nverts = simplified.size()/4;
					cont->verts = (int*)rcAlloc(sizeof(int)*cont->nverts*4, RC_ALLOC_PERM);
					if (!cont->verts)
					{
						ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'verts' (%d).", cont->nverts);
						return false;
					}
					memcpy(cont->verts, &simplified[0], sizeof(int)*cont->nverts*4);
					if (borderSize > 0)
					{
						// If the heightfield was build with bordersize, remove the offset.
						for (int j = 0; j < cont->nverts; ++j)
						{
							int* v = &cont->verts[j*4];
							v[0] -= borderSize;
							v[2] -= borderSize;
						}
					}
					
					cont->nrverts = verts.size()/4;
					cont->rverts = (int*)rcAlloc(sizeof(int)*cont->nrverts*4, RC_ALLOC_PERM);
					if (!cont->rverts)
					{
						ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'rverts' (%d).", cont->nrverts);
						return false;
					}
					memcpy(cont->rverts, &verts[0], sizeof(int)*cont->nrverts*4);
					if (borderSize > 0)
					{
						// If the heightfield was build with bordersize, remove the offset.
						for (int j = 0; j < cont->nrverts; ++j)
						{
							int* v = &cont->rverts[j*4];
							v[0] -= borderSize;
							v[2] -= borderSize;
						}
					}
					
					cont->reg = reg;
					cont->area = area;
				}
			}
		}
	}
	
	// Merge holes if needed.
	if (cset.nconts > 0)
	{
		// Calculate winding of all polygons.
		rcScopedDelete<signed char> winding((signed char*)rcAlloc(sizeof(signed char)*cset.nconts, RC_ALLOC_TEMP));
		if (!winding)
		{
			ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'hole' (%d).", cset.nconts);
			return false;
		}
		int nholes = 0;
		for (int i = 0; i < cset.nconts; ++i)
		{
			rcContour& cont = cset.conts[i];
			// If the contour is wound backwards, it is a hole.
			winding[i] = calcAreaOfPolygon2D(cont.verts, cont.nverts) < 0 ? -1 : 1;
			if (winding[i] < 0)
				nholes++;
		}
		
		if (nholes > 0)
		{
			// Collect outline contour and holes contours per region.
			// We assume that there is one outline and multiple holes.
			const int nregions = chf.maxRegions+1;
			rcScopedDelete<rcContourRegion> regions((rcContourRegion*)rcAlloc(sizeof(rcContourRegion)*nregions, RC_ALLOC_TEMP));
			if (!regions)
			{
				ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'regions' (%d).", nregions);
				return false;
			}
			memset(regions, 0, sizeof(rcContourRegion)*nregions);
			
			rcScopedDelete<rcContourHole> holes((rcContourHole*)rcAlloc(sizeof(rcContourHole)*cset.nconts, RC_ALLOC_TEMP));
			if (!holes)
			{
				ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'holes' (%d).", cset.nconts);
				return false;
			}
			memset(holes, 0, sizeof(rcContourHole)*cset.nconts);
			
			for (int i = 0; i < cset.nconts; ++i)
			{
				rcContour& cont = cset.conts[i];
				// Positively would contours are outlines, negative holes.
				if (winding[i] > 0)
				{
					if (regions[cont.reg].outline)
						ctx->log(RC_LOG_ERROR, "rcBuildContours: Multiple outlines for region %d.", cont.reg);
					regions[cont.reg].outline = &cont;
				}
				else
				{
					regions[cont.reg].nholes++;
				}
			}
			int index = 0;
			for (int i = 0; i < nregions; i++)
			{
				if (regions[i].nholes > 0)
				{
					regions[i].holes = &holes[index];
					index += regions[i].nholes;
					regions[i].nholes = 0;
				}
			}
			for (int i = 0; i < cset.nconts; ++i)
			{
				rcContour& cont = cset.conts[i];
				rcContourRegion& reg = regions[cont.reg];
				if (winding[i] < 0)
					reg.holes[reg.nholes++].contour = &cont;
			}
			
			// Finally merge each regions holes into the outline.
			for (int i = 0; i < nregions; i++)
			{
				rcContourRegion& reg = regions[i];
				if (!reg.nholes) continue;
				
				if (reg.outline)
				{
					mergeRegionHoles(ctx, reg);
				}
				else
				{
					// The region does not have an outline.
					// This can happen if the contour becaomes selfoverlapping because of
					// too aggressive simplification settings.
					ctx->log(RC_LOG_ERROR, "rcBuildContours: Bad outline for region %d, contour simplification is likely too aggressive.", i);
				}
			}
		}
		
	}
	
	return true;
}
